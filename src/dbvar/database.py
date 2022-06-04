# Sqlite exchange
import sqlite3


class Databasevar:
    def __init__(self, df_variants, db, vcfname):
        self.db = df_variants.to_sql(
            "VCF_variants", self.conn, if_exists="replace", index=False
        )
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()
        pass

    def intialize_dfb(self):
        self.execute("#[INFO] Create DB for " + self.vcfname)
        self.conn.commit()
