import argparse

# Usage: python vcfvalidator.py -c ../config/trueconfig.json -i ../share/vcftest.vcf Scan


class Parseoptions:
    def __init__(self):
        self.args = self.parseargs()

    def extract(self):
        print(self.args)
        if self.args.command == "Annotate":
            # for items
            dico = {"add": [], "remove": [], "edit": {}}
            if self.args.add:
                # add new row could be pass more than one row to add format 6fields(field,ID,number,type,description,valuestoadd), 6 fields and so on
                for values in self.args.add.split(";"):
                    dico["add"].append(values.split(",", 5))
            if self.args.remove:
                # remove row could be pass more than one row to remove format 2fields(field, ID)
                for values in self.args.remove.split(";"):
                    dico["remove"].append(values.split(","))
            if self.args.edit:
                # remove row could be pass more than one row to edit format 2fields(field.ID) to find which row to edit
                # then add field needed to be replace and the new entry like so INFO.gene|Type:Integer;   (change old type to integer)
                # or for multiple field INFO.gene|Type:Integer,Descritpion:newdescription;....
                edit = dico["edit"]
                # split if there are more than one edit
                for n_rows in self.args.edit.split(";"):
                    rows = n_rows.split("|")[0]
                    values = n_rows.split("|")[1]
                    edit[rows] = {}
                    for field in values.split(","):
                        key = field.split(":")[0]
                        val = field.split(":")[1]
                        edit[rows][key] = val

            # dico["vcf"] = self.args.vcf
            # dico["config"] = self.args.config
            return dico

    def parseargs(self):
        parser = argparse.ArgumentParser(
            description="VCFvalidator, manipulate variant call format file, basic: fix vcf both header and variant regaridng samtools spcifications version 4.2",
        )
        parser.add_argument(
            "-i",
            "--vcf",
            type=str,
            default="",
            help="Absolute path of input file, should be vcf but may miss columns like sample or format, script will be handle this, REQUIRED",
        )
        parser.add_argument(
            "-c",
            "--config",
            type=str,
            default="/home1/BAS/lamouchj/scripts/vcfvalidator/config/trueconfig.json",
            help="Absolute path of config file containing value allowed in vcfvalidator for each field both header and variant part",
        )
        parser.add_argument(
            "-o",
            "--output",
            type=str,
            default="basic",
            help="Output vcf with absolute path",
        )
        parser.add_argument(
            "-db",
            "--database",
            action="store_true",
            help="If table with same name are already present a new one will be create, set arg if you want to use an old one",
        )
        subparsers = parser.add_subparsers(dest="command")

        # VCF Annotate -------------------
        parser_annotate = subparsers.add_parser(
            "Annotate",
            help="Add, edit or remove row in VCF header"
            # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser_annotate.add_argument(
            "-a",
            "--add",
            default="",
            type=str,
            help="value to add in header",
        ),
        parser_annotate.add_argument(
            "-e",
            "--edit",
            default="",
            type=str,
            help="edit value in header",
        ),
        parser_annotate.add_argument(
            "-rm",
            "--remove",
            default="",
            type=str,
            help="remove row in header by field, ID pair",
        )
        parser_annotate.add_argument(
            "-db",
            "--dbname",
            type=str,
            help="Name of variant database, if not provided it will be set at input vcf name",
        )
        parser_annotate.add_argument(
            "-t",
            "--tablename",
            type=str,
            help="name of table database",
        )

        # VCF Scan ----------------
        parser_scan = subparsers.add_parser(
            "Scan",
            help="Identify anomalies in vcf both header and variants"
            # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser_scan.add_argument(
            "-db",
            "--dbname",
            type=str,
            help="Name of variant database, if not provided it will be set at input vcf name",
        )
        parser_scan.add_argument(
            "-t",
            "--tablename",
            type=str,
            help="name of table database",
        )

        # VCF Correct ----------------
        parser_correct = subparsers.add_parser(
            "Correct",
            help="Correct VCF header regarding variants annotations",
            # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser_correct.add_argument(
            "-db",
            "--dbname",
            type=str,
            help="Name of variant database, if not provided it will be set at input vcf name",
        )
        parser_correct.add_argument(
            "-t",
            "--tablename",
            type=str,
            help="name of table database",
        )

        # VCF Database and keep header to manipulate with other soft ------------

        parser_db = subparsers.add_parser(
            "Database",
            help="Create DB file regarding vcf and keep header in file, <vcf>_header.txt",
        )
        parser_db.add_argument(
            "-s",
            "--save",
            help="save header in input older extension _header.txt",
            type=str,
        )
        parser_db.add_argument(
            "-m",
            "--merge",
            help="header and db path to merge ex :<absoluteheader>,<absolutepathdb>",
            type=str,
        )
        parser_db.add_argument(
            "-db",
            "--dbname",
            type=str,
            help="Name of variant database, if not provided it will be set at input vcf name",
        )
        parser_db.add_argument(
            "-t",
            "--tablename",
            type=str,
            help="name of table database",
        )
        parser_db.add_argument(
            "-sn",
            "--samplename",
            type=str,
            help="sample name in final header",
        )
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # parser_db.add_argument(        )

        args = parser.parse_args()
        return args

    def get_args(self):
        # print("#[INFO] args ", self.parseargs())
        # print(type(self.parseargs()))
        # print(self.parseargs())
        return self.extract(), self.parseargs()
