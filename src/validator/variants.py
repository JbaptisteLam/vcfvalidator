from utils.utils import (
    read_json,
    preprocess_vcf,
    var_to_dataframe
)

class VCFpreprocess:
    def __init__(self, vcf):
        self.vcf = vcf
        # self.config = config
        self.header, self.skiprows = preprocess_vcf(self.vcf)
        self.df = var_to_dataframe(self.vcf, self.skiprows, columns=False)

    def explode_columns(self):
        df = self.df.copy()
        # test = pd.DataFrame(df.row.str.split(";")[1].tolist())
        df.INFO.str.split(";").tolist()

    def get_values(self):
        # print(self.df)
        return self.df, self.header


class Checkvariants:
    def __init__(self, df, samplename):
        self.df = df
        self.samplename = samplename

    def check_col(self):
        miss = []
        colus = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ]
        for col in colus:
            if col not in self.df.columns:
                miss.append(col)
        if self.df.columns[-1] in colus and miss:
            print("WARNING missing SAMPLE " + " ".join(miss) + " columns")
            return True
        else:
            print("#[INFO] First 9 columns field are OK")
            return False

    def add_basic(self):
        if self.check_col:
            self.df["FORMAT"] = "GT:AD:DP"
            self.df[self.samplename] = "0/1:50,50:100"
            return self.df
        else:
            return False