import sqlite3
from sqlite3 import Error


class ThrombophiliaRefHomo:

    rsid_map:dict = {}

    def init(self, reporter, sql_insert:str):
        self.parent = reporter
        self.sql_insert:str = sql_insert


    def setup(self):
        sql:str = "SELECT rsid, risk_allele, weight FROM weight WHERE state = 'ref' AND zygosity = 'hom'"
        self.parent.thrombophilia_cursor.execute(sql)
        rows:list[tuple] = self.parent.thrombophilia_cursor.fetchall()
        for row in rows:
            self.rsid_map[row[0]] = {'exist':True, 'allele':row[1], 'weight':row[2]}

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

    def process_row(self, row):
        rsid:str = str(row['dbsnp__rsid'])
        if rsid == '':
            return

        if not rsid.startswith('rs'):
            rsid = "rs"+rsid

        item:dict = self.rsid_map.get(rsid)
        if item:
            self.rsid_map[rsid]['exist'] = False


    def end(self):
        for rsid in self.rsid_map:
            if self.rsid_map[rsid]['exist']:
                allele:str = self.rsid_map[rsid]['allele']
                genotype:str = allele+allele

                query_for_pv:str = f"SELECT p_value FROM weight WHERE rsid = '{rsid}' AND weight.allele='{allele}' AND weight.state = 'ref' AND weight.zygosity='hom'"
                self.parent.thrombophilia_cursor.execute(query_for_pv)

                pvalue:tuple = self.parent.thrombophilia_cursor.fetchone()
                if pvalue is None or "NoneType":
                    substring:str = ''
                    query:str = "SELECT rsids.risk_allele, gene, genotype, genotype_specific_conclusion, " \
                    " rsid_conclusion, weight.weight, pmids, population, weight.p_value" \
                    f" FROM rsids, weight WHERE rsids.rsid ='{rsid}' AND weight.rsid = '{rsid}' AND weight.state = 'ref' AND weight.zygosity='hom'"

                else:
                    pv:str = ''
                    pv += pvalue[0]
                    strings = pv.split("[PMID: ")
                    string:list = []
                    for i in strings:
                        if i != '':
                            string.append(i[0:8])
                    substring:str = 'AND (studies.pubmed_id='

                    if len(string) > 1:
                        substring += f"'{string[0]}'"
                        for i in range(1, len(string)):
                            if i == (len(string)-1):
                                substring += f"OR studies.pubmed_id='{string[i]}')"
                                break
                            substring += f"OR studies.pubmed_id='{string[i]}'"
                    else:
                        substring += f"'{string[0]}')"
                
                    query:str = "SELECT rsids.risk_allele, gene, genotype, genotype_specific_conclusion, " \
                    " rsid_conclusion, weight.weight, pmids, population, weight.p_value, studies.populations,  pubmed_id" \
                    f" FROM rsids, studies, weight WHERE rsids.rsid ='{rsid}' AND weight.rsid = '{rsid}' AND weight.state = 'ref'" \
                    f" weight.zygosity='hom'" + substring

                self.parent.thrombophilia_cursor.execute(query)
                rows:list[tuple]  = self.parent.thrombophilia_cursor.fetchall()
                study_design:str = ''
                if len(rows) != 0:
                    if pvalue is None or "NoneType":
                        task:tuple = (rsid, rows[0][1], allele, genotype, rows[0][4], rows[0][3], float(rows[0][5]), rows[0][6], rows[0][7], '',
                                rows[0][8], self.get_color(rows[0][5], 0.6))
                    else:
                        if len(rows) > 1:
                            for row in rows:
                                for i in range(len(string)):
                                    if string[i] == str(row[10]):
                                        study_design += '[PMID: '+ string[i] + "] " + row[8] + "\n"
                        else:
                            study_design += rows[0][9]
                        task:tuple = (rsid, rows[0][1], rows[0][0], genotype, rows[0][4], rows[0][3], float(rows[0][5]), rows[0][6], rows[0][7], study_design,
                            rows[0][8], self.get_color(rows[0][5], 0.6))
                    self.parent.longevity_cursor.execute(self.sql_insert, task)