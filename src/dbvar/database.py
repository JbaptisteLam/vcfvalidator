# Sqlite exchange
import sqlite3


class Databasevar:
    def __init__(self, df_variants, db, vcfname):
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()
        # self.db = df_variants.to_sql(
        #    "VCF_variants", self.conn, if_exists="replace", index=False
        # )
        self.vcfname = vcfname
        pass

    def request(self, request):
        print("#[INFO] SQL REQUEST " + request)
        self.c.execute(request)
        self.conn.commit()
