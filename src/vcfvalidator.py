# aim: validator de vcf avec vcf en db sqlite pour les variants

from os import remove
import pandas as pd
import subprocess
import sys
import re
from os.path import join as osj
from validator.parseargs import Parseoptions

# from sqlalchemy import create_engine

# example config
# config = {
#    "header": {
#        "field" {
#            "fileformat": "add or edit or remove"
#            #if INFO or FORMAT
#            "INFO": {
#                    action: add edit or remove  value =list   for fileformat and other value
#                ID: value
#                "Number": value
#                "Type": Value
#                "Description": Value
#                opt "Source": value,
#                opt "Version"
#           }
#        }
#    },
#    "fields": {},
#    "variants": {},
# }

# utils
def preprocess_vcf(file):
    """
    from vcf file as input
    return dico with 2 lists, header with each row from header and field which contain the vcf header field,
            and number of row of header (without field)
    """
    data = {}
    skip = []
    with open(file, "r") as f:
        data["header"] = []
        for i, lines in enumerate(f):
            data["header"].append(lines.strip())
            if lines.split("\t")[0] == "#CHROM":
                print(lines)
                data["fields"] = lines.strip().split("\t")
                skip.append(i)
                break
    return data, skip[0]


def systemcall(command):
    """
    Passing command to the shell, return first item in list containing stdout lines
    """
    try:
        print("#[INFO] " + command)
        p = subprocess.check_output([command], shell=True)
    except subprocess.CalledProcessError as err:
        sys.exit("ERROR " + str(err.returncode))
    return (
        subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
        .stdout.read()
        .decode("utf8")
        .strip()
        .split("\n")[0]
    )


def win_to_unix(input, output):
    return systemcall('awk \'{ sub("\r$", ""); print }\' ' + input + " > " + output)


class Checkheader:
    def __init__(self, header_dict, config, trueconfig):
        self.header = header_dict["header"]
        self.fields = header_dict["fields"]
        self.config = config
        # list of field id to process for each
        # self.rm = self.config["remove"]
        self.add = self.config["add"]
        self.trueconfig = trueconfig
        # self.edit = self.config["edit"]
        # self.rmwhole = self.config["rmall"]

    def check_integrity(self):
        pass

    def add_assembly(self):
        # need install of gatk
        return

    def add_row(self, field, id, number, type, description, **kwargs):
        self.header.append(
            "##"
            + field
            + "=<"
            + ",".join(
                [
                    "ID=" + id,
                    "Number=" + number,
                    "Type=" + type,
                    "Description=" + description,
                ]
            )
            + ">"
        )

    def remove_row(self):
        for rows in self.rm:
            index = self.matchingline()
            if index:
                del self.header[index]

    def matching_line(self, field, id):
        """
        return index of row in global header
        """
        matching = []
        for i, rows in enumerate(self.header):
            x = re.search(field + "=<ID=" + id)
            if x:
                matching.append(x)
        if not matching:
            print("WARNING " + field + " ," + id + " does not match any row in header")
        return matching[0]

    def edit_row(self):
        """
        edit only one row loop is in class.process func
        """
        return

    def remove_whole(self):
        print("#[INFO] Clear whole header")
        self.header.clear()

    def process(self):
        # if self.rm:
        #    self.removerow()
        if self.add:
            for rows in self.add:
                print(rows)
                self.addrow(*rows)
                print(self.header)
        # if self.edit:
        #    for j, r in self.edit:
        #        self.editrow()
        # if self.rmwhole:
        #    self.removewhole()


class VCFpreprocess:
    def __init__(self, vcf, config):
        self.vcf = vcf
        self.config = config
        self.header, self.skiprows = preprocess_vcf(self.vcf)
        self.df = self.vartodataframe(self.skiprows, columns=False)

    def var_to_dataframe(self, skiprows, columns):
        """
        from tsv file to pandas dataframe, skiprows number of row to reach header, columns: col which need change (from , to .) to allowed excel filter
        """
        if skiprows:
            df = pd.read_csv(
                self.vcf,
                skiprows=skiprows,
                sep="\t",
                header=0,
                chunksize=10000,
                low_memory=False,
            )
        else:
            df = pd.read_csv(
                self.vcf,
                sep="\t",
                header=0,
                chunksize=10000,
                low_memory=False,
            )
        # if "index" in df:
        #    df = df.drop(columns="index")
        if columns:
            for col in columns:
                if col in df.columns:
                    df[col] = df[col].apply(lambda x: x.replace(",", "."))
                    df[col] = df[col].astype(float)
        return df

    def get_values(self):
        return self.df, self.header

        # def dfTosql(df, con, dico):
        #    return


def main():
    opt = Parseoptions()
    # return options parsed in dico format and argparse.Namespace respectively
    dico_args, args = opt.get_args()
    # print(args)
    # get options parser and process

    # Load vcf file
    vcf = VCFpreprocess(args.vcf, args.config)
    variants, header = vcf.get_values()

    # Act on header vcf
    for actions, values in dico_args.items():
        #    #if user need adding row
        if actions == "add" and values:
            pass
    #    vcf = VCFpreprocess(args.vcf, arg)
    #    df, header = vcf.getvalues()
    #    dico = {"add": [arg]}
    #    pheader = Checkheader(header, dico)
    #    print(pheader.process())


if __name__ == "__main__":
    main()
