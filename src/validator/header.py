import re
from utils.utils import explode_header

class Checkheader:
    def __init__(self, header_dict, dico_args, trueconfig):
        self.header = header_dict
        self.dico_args = dico_args
        # list of field id to process for each
        # self.rm = self.config["remove"]
        # self.add = self.config["add"]
        self.config = trueconfig
        print(self.config)
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
                    #print("WARNING "+elem+" not a correct value field", self.id_issues(rows, self.header['header']))
                    pass
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
                print(values)
                for rows in values:
                    self.remove_row(rows)

        return self.header
    
    def id_issues(self, lines, iterable):
        return "rows "+str(iterable.index(lines)+1)+" "+ lines
    
    def header_check(self):
        match = []
        #If json is malforammted it header got a problem
        header_explode, error = explode_header(self.header)
        if error:
            for lines in error:
                self.raise_integrity('Description', lines, True, 'malformation')
        else:
            print("#[INFO] Header properly formated")
        #check if field are allowed
        for lines in self.header['header']:
            self._fields(lines)
            try:
                str(lines)
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

