# aim: validator de vcf avec vcf en db sqlite pour les variants

from multiprocessing.sharedctypes import Value
from os import remove
import argparse
import pandas as pd
import subprocess
import sys
import re
from os.path import join as osj

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
def preprocessvcf(file):
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


def wintounix(input, output):
    return systemcall('awk \'{ sub("\r$", ""); print }\' ' + input + " > " + output)


class Parseoptions:
    def __init__(self, config):
        self.config = config

    def extract(self):
        return self.config.split(",")


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

    def addassembly(self):
        # need install of gatk
        return

    def addrow(self, field, id, number, type, description, **kwargs):
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

    def removerow(self):
        for rows in self.rm:
            index = self.matchingline()
            if index:
                del self.header[index]

    def matchingline(self, field, id):
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

    def editrow(self):
        """
        edit only one row loop is in class.process func
        """
        return

    def removewhole(self):
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
        self.header, self.skiprows = preprocessvcf(self.vcf)
        self.df = self.vartodataframe(self.skiprows, columns=False)

    def vartodataframe(self, skiprows, columns):
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

    def getvalues(self):
        return self.df, self.header
        # def dfTosql(df, con, dico):
        #    return


def parseargs():
    parser = argparse.ArgumentParser(
        description="Filter tsv in DPNI and POOL context, basically tsv have to come from varank analysis "
    )
    parser.add_argument(
        "-i",
        "--vcf",
        type=str,
        help="Absolute path of input file, should be vcf but may miss columns like sample or format, script will be handle this",
    )
    subparsers = parser.add_subparsers(dest="command")
    parser_add = subparsers.add_parser(
        "add", help="a help", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_add.add_argument(
        "-v",
        "--values",
        type=str,
        help="value to add in header",
    )

    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    opt = Parseoptions(args.values)
    arg = opt.extract()
    if args.command == "add":
        vcf = VCFpreprocess(args.vcf, arg)
        df, header = vcf.getvalues()
        dico = {"add": [arg]}
        pheader = Checkheader(header, dico)
        print(pheader.process())


if __name__ == "__main__":
    main()
