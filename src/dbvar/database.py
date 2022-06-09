# Sqlite exchange
import sqlite3


class Databasevar:
    def __init__(self, df_variants, db, vcfname, tablename, config):
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()
        self.vcfname = vcfname
        self.db = db
        self.df_variants = df_variants
        self.tablename = tablename
        self.config = config
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

    def chromosome_check(self, mode):
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
        if len(res) != len(self.df_variants.index):
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
            else:
                print("WARNING some variants are not conform ", res_chr)
                print("#[INFO] correct CHR field")
                for values in res_chr:
                    self.update_value(values)

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

    def db_to_dataframe(self):
        return self.c.fetchall()
