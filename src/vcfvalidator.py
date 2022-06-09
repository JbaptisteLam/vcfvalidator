# aim: validator de vcf avec vcf en db sqlite pour les variants


import os
import re
import sys
from os import remove
from os.path import join as osj
from utils.utils import (
    read_json,
    preprocess_vcf,
    var_to_dataframe,
    parse_sample_field,
    parse_info_field,
    systemcall,
    df_to_vcflike,
    create_vcf,
    explode_header,
    is_utf8
)
import pandas as pd
from dbvar.database import Databasevar as dbv
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


class Checkheader:
    def __init__(self, header_dict, dico_args, trueconfig):
        self.header = header_dict
        self.dico_args = dico_args
        # list of field id to process for each
        # self.rm = self.config["remove"]
        # self.add = self.config["add"]
        self.config = trueconfig
        # self.edit = self.config["edit"]
        # self.rmwhole = self.config["rmall"]

    def check_integrity(self):
        #print(self.dico_args)
        for key, val in self.dico_args.items():
            if val:
                if key == "edit":
                    for j, val in val.items():
                        self.raise_integrity(j.split(".")[0], False,  False, None)
                else:
                    for rows in val:
                        self.raise_integrity(rows[0], False,  False, None)

    def raise_integrity(self, elem, rows, check, type):
        if (
            elem not in self.config["header"]["field"]
            and elem not in self.config["header"]["extrafield"]
        ):
            if not check:
                raise ValueError(
                    "'"
                    + elem
                    + "' is not a correct value field \n\t Correct values are "
                    + ",".join(self.config["header"]["field"])
                    + ",".join(self.config["header"]["extrafield"])
                )
            else:
                if type == 'correctvalue':
                    print("WARNING "+elem+" not a correct value field", self.id_issues(rows, self.header['header']))
                elif type == 'malformation':
                    print("WARNING "+elem+" field is not correctly formated EOL", self.id_issues(rows, self.header['header']))

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

    def correct_header(self):
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

        return self.header
    
    def id_issues(self, lines, iterable):
        return "rows "+str(iterable.index(lines)+1)+" "+ lines
    
    def header_check(self):
        match = []
        #If json is malforammted it header got a problem
        header_explode, error = explode_header(self.header)
        for lines in error:
            self.raise_integrity('Description', lines, True, 'malformation')
        #check if field are allowed
        for lines in self.header['header']:
            self._fields(lines)
            try:
                str(lines)
                #print('ok')
            except SyntaxError:
                match.append(self.id_issues(lines))
        return match
    
    def _fields(self, rows):
        #dico_values = get_header_id(self.header, self.config)
        #dico_values = explode_header(self.header)
        field = re.findall(r'(?<=##)[^=]+', rows)[0]
        #header values not in config json allowed
        if field not in self.config["header"]["field"] and field not in self.config["header"]["extrafield"]:
            self.raise_integrity(field, rows, True, 'correctvalue')
            #print("WARNING "+field+" is not an allowed value ", self.id_issues(field, self.header['header']))


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
            return True
        else:
            return False

    #def process(self):
    #    if self.check_col():
    #        self.add_basic()
    #        return self.df
    #    else:
    #        return pd.DataFrame({"emp1": [], "empt2": []})


def main():
    opt = Parseoptions()
    # return options parsed in dico format and argparse.Namespace respectively
    dico_args, args = opt.get_args()

    #read config
    config = read_json(args.config)

    # if no output specified generate vcf as input file name with correct
    if args.output == "basic":
        output = osj(
            os.path.dirname(args.vcf),
            os.path.basename(args.vcf).split(".")[0] + "_correct.vcf",
        )
    if os.path.exists(output):
        print("ERROR " + output + " exist remove or choose an other output filename")
        exit()

    # Load vcf file
    variants, header = VCFpreprocess(args.vcf).get_values()

    # choice main func
    if args.command == "Correct":
        main_correct(variants, header, args, dico_args, output, config)
    elif args.command == "Scan":
        main_scan(variants, header, args, config)
    elif args.command == "Annotate":
        main_annotate(variants, header, args, output, config)
    else:
        main_test(header, config)


def main_annotate(variants, header, args, dico_args, output, config):
    if dico_args:
        headercheck = Checkheader(header, dico_args, args.config)
        hprocess = headercheck.correct_header()

    else:
        hprocess = header

    create_vcf(variants, hprocess, output)

def main_scan(variants, header, args, config):
    print("Check integrity of "+os.path.abspath(args.vcf)+" ...")
    errors = is_utf8(args.vcf)
    if errors:
        for err in errors:
            print(err)
    else:
        print("#[INFO] ASCII Check OK")
    cvar = Checkvariants(variants, "sample")
    variants_new = cvar.check_col()

    headercheck = Checkheader(
            header, {}, config
        ).header_check()
    explode_header(header)



def main_correct(variants, header, args, output, config):
    #only if correct mode activate #TODO
    #print("Check integrity of "+args.vcf+" ...")
    # add correct column
    # header with basic GT AD and DP add
    cvar = Checkvariants(variants, "sample")
    #return True if basic columns where add false otherwise
    variants_new = cvar.add_basic()

    if not variants_new:
        headercheck = Checkheader(
            header, config["variants"]["basic"], args.config
        ).correct_header()
    #explode INFO and SAMPLE field
    variants_explode, badannno, list_sample = parse_sample_field(
        parse_info_field(variants)
    )
    # worked
    # tablename = "variants"
    if args.dbname:
        dbname = args.dbname
    else:
        dbname = os.path.join("dbvar", os.path.basename(args.vcf).split(".")[0] + ".db")
    db = dbv(
        variants_explode,
        dbname,
        os.path.basename(args.vcf),
        args.tablename,
        config,
    )
    #FROM dataframe to sql db
    db.create_table()

    #chromosome check only if correct mode activate
    db.chromosome_check()

    #Generate vcf only if correction mode is enable
    #from db sql to dataframe after managing
    dfwrite_tmp = pd.DataFrame(db.db_to_dataframe(), columns=variants_explode.columns)
    dfwrite = df_to_vcflike(dfwrite_tmp, "Likely_benin")

    #Regenerate vcf final 
    create_vcf(dfwrite, header, output)

def main_analyze(variants, header, args, output):
    '''
    check many possible malformation in vcf
    '''
    pass

def main_test(header, config):
    #print(get_header_id(header, config))
    pass


if __name__ == "__main__":
    main()
