from oakvar import BasePostAggregator
from pathlib import Path
import sys
cur_path = str(Path(__file__).parent)
sys.path.append(cur_path)
import sqlite3
import thrombophilia_ref_homo


class CravatPostAggregator (BasePostAggregator):
    sql_insert = """ INSERT INTO thrombophilia (
                        rsid,
                        gene,
                        risk_allele,
                        genotype,
                        conclusion,
                        genotype_conclusion,
                        weight,
                        pmid,
                        population,
                        studydesign,
                        pvalue,
                        weightcolor
                    ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?) """
    ref_homo = thrombophilia_ref_homo.ThrombophiliaRefHomo()

    def check(self):
        return True

    def setup (self):
        self.ref_homo.init(self, self.sql_insert)
        modules_path = str(Path(__file__).parent)
        sql_file = modules_path + "/data/thrombophilia.sqlite"
        if Path(sql_file).exists():
            self.thrombophilia_conn = sqlite3.connect(sql_file)
            self.thrombophilia_cursor = self.thrombophilia_conn.cursor()

        self.result_path = Path(self.output_dir, self.run_name + "_longevity.sqlite")
        self.longevity_conn = sqlite3.connect(self.result_path)
        self.longevity_cursor = self.longevity_conn.cursor()
        sql_create = """ CREATE TABLE IF NOT EXISTS thrombophilia (
            id integer NOT NULL PRIMARY KEY,
            rsid text,
            gene text,
            risk_allele text,
            genotype text,
            conclusion text,
            genotype_conclusion text,
            weight float,
            pmid text,
            population text,
            studydesign text,
            pvalue text,
            weightcolor text
            )"""
        self.longevity_cursor.execute(sql_create)
        self.longevity_conn.commit()
        self.longevity_cursor.execute("DELETE FROM thrombophilia;")
        self.ref_homo.setup()


    def cleanup (self):
        if self.longevity_cursor is not None:
            self.longevity_cursor.close()
        if self.longevity_conn is not None:
            self.longevity_conn.commit()
            self.longevity_conn.close()
        if self.thrombophilia_cursor is not None:
            self.thrombophilia_cursor.close()
        if self.thrombophilia_conn is not None:
            self.thrombophilia_conn.close()
        return


    def get_color(self, w, scale = 1.5):
        w = float(w)
        if w < 0:
            w = w * -1
            w = 1 - w * scale
            if w < 0:
                w = 0
            color = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = "ff" + color + color
        else:
            w = 1 - w * scale
            if w < 0:
                w = 0
            color = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = color + "ff" + color

        return color

    def annotate (self, input_data):
        rsid = str(input_data['dbsnp__rsid'])
        if rsid == '':
            return

        self.ref_homo.process_row(input_data)

        if not rsid.startswith('rs'):
            rsid = "rs" + rsid

        alt = input_data['base__alt_base']
        ref = input_data['base__ref_base']

        zygot = input_data['vcfinfo__zygosity']
        genome = alt + ref
        gen_set = {alt, ref}
        if zygot == 'hom':
            genome = alt + alt
            gen_set = {alt, alt}

        zygot:str = input_data['vcfinfo__zygosity']
        if zygot is None or zygot == "":
            zygot = "het"

        query_for_pv:str = f"SELECT p_value FROM weight WHERE rsid = '{rsid}' AND weight.allele='{alt}' AND weight.zygosity='{zygot}'"

        self.thrombophilia_cursor.execute(query_for_pv)
        pvalue:tuple = self.thrombophilia_cursor.fetchone()
        if pvalue is None:
            substring:str = ''
            query:str = "SELECT rsids.risk_allele, gene, genotype, genotype_specific_conclusion, " \
            " rsid_conclusion, weight.weight, pmids, studies.population, weight.p_value, pubmed_id" \
            f" FROM rsids, weight WHERE rsids.rsid ='{rsid}' AND weight.rsid = '{rsid}' " \
            f" AND weight.allele='{alt}' AND weight.zygosity='{zygot}' "

        else:
            pv:str = ''
            pv += pvalue[0]
            strings = pv.split("[PMID: ")
            string = []
            for i in strings:
                if i != '':
                    string.append(i[0:8])
            substring = 'AND (studies.pubmed_id='

            if len(string) > 1:
                substring += f"'{string[0]}'"
                for i in range(1, len(string)):
                    if i == (len(string)-1):
                        substring += f"OR studies.pubmed_id='{string[i]}')"
                        break
                    substring += f"OR studies.pubmed_id='{string[i]}'"
            else:
                substring += f"'{string[0]}'"
        
            query:str = "SELECT rsids.risk_allele, gene, genotype, genotype_specific_conclusion, " \
            " rsid_conclusion, weight.weight, pmids, population, studies.populations, weight.p_value, pubmed_id" \
            f" FROM rsids, studies, weight WHERE rsids.rsid ='{rsid}' AND weight.rsid = '{rsid}' " \
            f" AND weight.allele='{alt}' AND weight.zygosity='{zygot}'" + substring

        self.thrombophilia_cursor.execute(query)
        rows:list[tuple] = self.thrombophilia_cursor.fetchall()

        if len(rows) == 0:
            return

        study_design=''
        if len(rows) > 1:
            for row in rows:
                for i in range(len(string)):
                    if string[i] == str(row[10]):
                        study_design += '[PMID: '+ string[i] + "] " + row[8] + "\n"
        else:
            study_design += rows[0][8]

        row_gen = {rows[0][2][0], rows[0][2][1]}
        if pvalue is None:
            task = (rsid, rows[0][1], rows[0][0], genome, rows[0][4], rows[0][3], float(rows[0][5]), rows[0][6], rows[0][7], '',
                    rows[0][8], self.get_color(rows[0][5], 0.6))
        else:
            task = (rsid, rows[0][1], rows[0][0], genome, rows[0][4], rows[0][3], float(row[5]), rows[0][6], rows[0][7], study_design,
                rows[0][9], self.get_color(rows[0][5], 0.6))

        if gen_set == row_gen:
            self.longevity_cursor.execute(self.sql_insert, task)

        return {"col1":""}


    def postprocess(self):
        self.ref_homo.end()
        pass

