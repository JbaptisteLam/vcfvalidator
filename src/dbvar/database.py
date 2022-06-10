# Sqlite exchange
import sqlite3
import ssl
from textwrap import indent
import json
from tqdm import tqdm
from utils.utils import (
    parse_sample_field,
    parse_info_field
)

class Databasevar:
    def __init__(self, df_variants, db, vcfname, tablename, config, header_explode):
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()
        self.vcfname = vcfname
        self.db = db
        self.df_variants = df_variants
        self.tablename = tablename
        self.config = config
        try:
            self.header = header_explode[0]
        except TypeError:
            self.header = False
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

    def table_exists(self):
        #print("SELECT name FROM sqlite_temp_master WHERE type='table' AND name='"+self.tablename+"'")
        #self.c.execute("SELECT name FROM sqlite_temp_master WHERE type='table' AND name='"+self.tablename+"'")
        self.c.execute('SELECT name from sqlite_master where type= "table"')
        for tables in self.c.fetchall():
            if tables[0] == self.tablename:
                print("#[INFO] Table '"+self.tablename+"' already exists, remove if the db need cleaning args TODO")
                return True
        return False
    
    def explode_annotations(self):
        df_variants, badanno, list_sample = parse_sample_field(
        parse_info_field(self.df_variants))
        return df_variants, badanno, list_sample


    def create_table(self):
        #self.df_variants.to_sql(
        #    self.tablename, self.conn, if_exists="replace", index=False
        #)
        print("#[INFO] Create Table "+self.tablename)
        self.insert_with_progress(self.explode_annotations()[0], self.db)

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
                return res_chr
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

    @staticmethod
    def chunker(seq, size):
        # from http://stackoverflow.com/a/434328
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    def insert_with_progress(self, df, dbfile):
        #from https://stackoverflow.com/a/39495229
        con = sqlite3.connect(dbfile)
        #case of less than 10 var in dataframe
        chunksize = int(len(df) / 10) # 10%
        if chunksize < 1:
            chunksize = 1
        with tqdm(total=len(df), desc="#[INFO] From dataframe to DB, chunksize: "+str(chunksize)) as pbar:
            for i, cdf in enumerate(self.chunker(df, chunksize)):
                replace = "replace" if i == 0 else "append"
                cdf.to_sql(con=con, name=self.tablename, if_exists=replace, index=False)
                pbar.update(chunksize)


    def get_col_name(self):
        self.c.execute('SELECT * from '+self.tablename)
        names = [description[0] for description in self.c.description]
        return names

    #Edit INFO fields annotations 
    def db_to_dataframe(self):
        return self.c.fetchall()
    
    def correct_variants(self):
        ## Act on header vcf
        # check if value are correct
        for actions, values in self.dico_args.items():
            ##if user need adding row
            #if actions == "add" and values:
            #    for rows in values:
            #        if not self.matching_line(rows[0], rows[1]):
            #            self.add_row(*rows)
            #if actions == "edit" and values:
            #    for key, rows in values.items():
            #        self.edit_row(key, rows)
            if actions == "remove" and values:
                for rows in values:
                    self.remove_row(*rows)
        return self.header

    def add_info(self):
        pass

    def rm_info(self, col):
        self.c.execute("ALTER TABLE "+self.tablename+" DROP COLUMN "+col)
        self.conn.commit()

    def edit_info(self):
        pass