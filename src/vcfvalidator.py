# aim: validator de vcf avec vcf en db sqlite pour les variants

from os import remove
import pandas as pd
import subprocess
import sys
import json
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
            if lines.split("\t")[0] != "#CHROM":
                data["header"].append(lines.strip())
            else:
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


def read_json(file):
    with open(file, "r") as jsonFile:
        jconf = json.load(jsonFile)
        return jconf


def win_to_unix(input, output):
    return systemcall('awk \'{ sub("\r$", ""); print }\' ' + input + " > " + output)


class Checkheader:
    def __init__(self, header_dict, dico_args, trueconfig):
        self.header = header_dict
        self.dico_args = dico_args
        # list of field id to process for each
        # self.rm = self.config["remove"]
        # self.add = self.config["add"]
        self.config = read_json(trueconfig)
        # self.edit = self.config["edit"]
        # self.rmwhole = self.config["rmall"]

    def check_integrity(self):
        print(self.dico_args)
        for key, val in self.dico_args.items():
            if val:
                if key == "edit":
                    for j, val in val.items():
                        self.raise_integrity(j.split(".")[0])
                else:
                    for rows in val:
                        self.raise_integrity(rows[0])

    def raise_integrity(self, elem):
        if (
            elem not in self.config["header"]["field"]
            and elem not in self.config["header"]["extrafield"]
        ):

            raise ValueError(
                "'"
                + elem
                + "' is not a correct value field \n\t Correct values are "
                + ",".join(self.config["header"]["field"])
                + ",".join(self.config["header"]["extrafield"])
            )

    def add_assembly(self):
        # need install of gatk
        return

    def add_row(self, field, id, number, type, description, **kwargs):
        row_add = (
            "##"
            + field
            + "=<"
            + ",".join(
                [
                    "ID=" + id,
                    "Number=" + number,
                    "Type=" + type,
                    'Description="' + description,
                ]
            )
            + '">'
        )
        print("#[INFO] Adding " + row_add + " in header")
        self.header["header"].append(row_add)

    def remove_row(self):
        index = self.matchingline()
        if index:
            print("#[INFO] Removing " + index + " from header")
            del self.header[index]

    def matching_line(self, field, id):
        """
        return index of row in global header
        """
        matching = []
        for i, rows in enumerate(self.header["header"]):
            # print("look ", ".*" + field + "=<ID=" + id + ".*")
            x = re.findall(".*" + field + "=<ID=" + id + ".*", rows)
            if x:
                matching.extend(x)
        if not matching:
            return []
        elif len(matching) > 1:
            print(
                "WARNING '"
                + field
                + "."
                + id
                + "' got more than one match, verify your header"
            )
            return []

        return matching

    def edit_row(self, id, fields):
        """
        edit only one row loop is in class.process func, id field.ID   fieldskey: fieldsnewvalue,
        """
        key = id.split(".")[0]
        val = id.split(".")[1]
        # if more than one line matching or no matching return empty list and do not process anything
        row_match = self.matching_line(key, val)
        # print(row_match)
        if row_match:
            # print(id, fields)
            dico_val = {}
            # explode old header row
            for f in re.search(r"<(.*?)>", row_match[0]).group(1).split(","):
                dico_val[f.split("=")[0]] = f.split("=")[1]
            # replace by new value
            for k, v in fields.items():
                dico_val[k] = v
                if k == "Description":
                    dico_val[k] = '"' + v + '"'
            # reconstruct row string from modify dict
            row_process = (
                "##"
                + key
                + "=<"
                + ",".join(["=".join([s, r]) for s, r in dico_val.items()])
                + ">"
            )
            # print(row_process)
            print("#[INFO] Edit index " + key + "." + val + " return " + row_process)
            # print(dico_val)
            self.header["header"] = list(
                map(
                    lambda x: x.replace(row_match[0], row_process),
                    self.header["header"],
                )
            )
        else:
            print(
                "WARNING '"
                + key
                + "."
                + val
                + "' does not match any row in header PASS"
            )

    def remove_whole(self):
        print("#[INFO] Clear whole header")
        self.header.clear()

    def process(self):
        ## Act on header vcf
        # check if value are correct
        self.check_integrity()
        for actions, values in self.dico_args.items():
            ##if user need adding row
            if actions == "add" and values:
                for rows in values:
                    if not self.matching_line(rows[0], rows[1]):
                        self.add_row(*rows)
            if actions == "edit" and values:
                for key, rows in values.items():
                    self.edit_row(key, rows)
            if actions == "remove" and values:
                for rows in values:
                    self.remove_row(*rows)

        print(self.header)
        # if self.rmwhole:
        #    self.removewhole()


class VCFpreprocess:
    def __init__(self, vcf):
        self.vcf = vcf
        # self.config = config
        self.header, self.skiprows = preprocess_vcf(self.vcf)
        self.df = self.var_to_dataframe(self.skiprows, columns=False)

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
    vcf = VCFpreprocess(args.vcf)
    variants, header = vcf.get_values()

    headercheck = Checkheader(header, dico_args, args.config)
    headercheck.process()

    #    vcf = VCFpreprocess(args.vcf, arg)
    #    df, header = vcf.getvalues()
    #    dico = {"add": [arg]}
    #    pheader = Checkheader(header, dico)
    #    print(pheader.process())


if __name__ == "__main__":
    main()
