import argparse


class Parseoptions:
    def __init__(self):
        self.args = self.parseargs()

    def extract(self):
        # for items
        dico = {"add": [], "remove": [], "edit": {}}
        if self.args.add:
            # add new row could be pass more than one row to add format 5fields(field,ID,number,type,description), 5 fields and so on
            dico["add"].append(self.args.add.split(";"))
        if self.args.remove:
            # remove row could be pass more than one row to remove format 2fields(field, ID)
            dico["remove"].append(self.args.add.split(";"))
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
            description="Filter tsv in DPNI and POOL context, basically tsv have to come from varank analysis ",
        )
        parser.add_argument(
            "-i",
            "--vcf",
            type=str,
            required=True,
            help="Absolute path of input file, should be vcf but may miss columns like sample or format, script will be handle this",
        )
        parser.add_argument(
            "-c",
            "--config",
            type=str,
            default="/home1/BAS/lamouchj/scripts/vcfvalidator/config/trueconfig.json",
            help="Absolute path of config file containing value allowed in vcfvalidator for each field both header and variant part",
        )
        subparsers = parser.add_subparsers(dest="command")

        # VCF Annotate -------------------
        parser_annotate = subparsers.add_parser(
            "Annotate",
            help="Add, edit or remove row in VCF header",
            # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser_annotate.add_argument(
            "-a",
            "--add",
            type=str,
            help="value to add in header",
        ),
        parser_annotate.add_argument(
            "-e",
            "--edit",
            type=str,
            help="edit value in header",
        ),
        parser_annotate.add_argument(
            "-rm",
            "--remove",
            type=str,
            help="remove row in header by field, ID pair",
        )

        # VCF Correct ----------------
        parser_correct = subparsers.add_parser(
            "Correct",
            help="Correct VCF header regarding variants annotations",
            # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

        args = parser.parse_args()
        return args

    def get_args(self):
        # print("#[INFO] args ", self.parseargs())
        # print(type(self.parseargs()))
        return self.extract(), self.parseargs()
