# Sqlite exchange
import sqlite3
import ssl
from textwrap import indent
import json


class Databasevar:
    def __init__(self, df_variants, db, vcfname, tablename, config, header_explode):
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()
        self.vcfname = vcfname
        self.db = db
        self.df_variants = df_variants
        self.tablename = tablename
        self.config = config
        self.header = header_explode[0]
        print("#[INFO] SQLite " + db + " connected ")

    def request(self, request):
        try:
            print("#[INFO] SQL REQUEST " + request)
            self.c.execute(request)
            self.conn.commit()
        except sqlite3.Error as error:
            print("Failed request " + request, error)
        finally:
            if self.conn:
                self.conn.close()
                print("#[INFO] The SQLite connection is closed")

    def update(self):
        pass

    def create_table(self):
        self.df_variants.to_sql(
            self.tablename, self.conn, if_exists="replace", index=False
        )

    def chromosome_check(self):
        # request = "SELECT * FROM " + self.tablename + " WHERE #CHROM LIKE ", (chr,)
        # self.request(request)
        self.c.execute(
            "SELECT * from " + self.tablename + " WHERE [#CHROM] Like ?",
            ("chr%",),
        )
        res = self.c.fetchall()
        # print(len(res))
        # print(type(self.c.fetchall()))
        # print(len(self.df_variants.index))
        #If not all variants got chr in #CHROM column
        #print(len(res))
        #print(len(self.df_variants.index))
        if len(res) != len(self.df_variants.index) and len(res) != 0:
            self.c.execute(
                "SELECT [#CHROM], POS from "
                + self.tablename
                + " WHERE [#CHROM] NOT IN ({seq})".format(
                    seq=",".join(["?"] * len(self.config["variants"]["chr"]))
                ),
                self.config["variants"]["chr"],
            )
            res_chr = self.c.fetchall()
            if not res_chr:
                # print(res_chr)
                print("#[INFO] Chromosomes check OK")
                return False
            else:
                print("WARNING some variants do not have 'chr' in #CHROM col ", res_chr)
                return self.c.fetchall()
        #All variants do not have chr in #CHROM column
        elif len(res) == 0:
            print("WARNING any variants got 'chr' in #CHROM col")
            return self.c.execute(
                "SELECT [#CHROM], POS from "
                + self.tablename)





    def update_value(self, old):
        self.c.execute(
            "UPDATE "
            + self.tablename
            + " SET [#CHROM] = 'chr"
            + str(old[0])
            + "' WHERE [#CHROM] = '"
            + str(old[0])
            + "' AND POS = '"
            + str(old[1])
            + "'"
        )
        self.conn.commit()
        self.c.execute("SELECT * FROM " + self.tablename)

    def check_annotations(self):
        miss = []
        curs = self.c.execute('SELECT * from '+self.tablename)
        names = [desc[0] for desc in curs.description]
        for col in names:
            if not col in self.header.keys() and not col in self.header['FORMAT'] and not col in self.header['INFO'] \
                and not col in self.config['header']['all']:
                miss.append(col)   
        return miss

    def db_to_dataframe(self):
        return self.c.fetchall()
