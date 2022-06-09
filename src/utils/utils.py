from posixpath import sep
import pandas as pd
import json
import subprocess
import numpy as np
import os
import sys
import re
from itertools import zip_longest

# git combo FF, dossier TEST TODO
def fancystdout(style, text):
    subprocess.call("pyfiglet -f " + style + " -w 100 " + text)

def launch():
    message = '''\n
    <> Author: Jean-Baptiste Lamouche
    <> Mail: Jbaptiste.lamouche@gmail.com
    <> Github: https://github.com/JbaptisteLam
    <> Version: 0.1\n
    '''
    return message

def is_utf8(vcf):
    error = []
    with open(vcf, 'r') as f:
        for i, line in enumerate(f):
            line = line.encode()
            try:
                line.decode('ascii')
                #print("#[INFO] Valid utf-8")
                #print(line.decode('utf-8'))
            except UnicodeDecodeError:
                error.append('WARNING row '+str(i+1)+' invalid utf-8')
    return error

#DEPRECIATED TODO
def get_header_id(header, config):
    header_infos = {}
    for rows in header['header']:
        field = re.findall(r'(?<=##)[^=]+', rows)[0]
        if field and field in config['header']['field'] and field != 'ID':
            if not field in header_infos:
                header_infos[field] = []
            id = re.findall(r'(?<=ID=)[^,]+', rows)[0]
            #print(id)
            header_infos[field].append(id)
    #print(header_infos)
    return header_infos

def explode_header(header):
    '''
    from header of vcf file return huge dico separate first field second fiel list of values
    '''
    dico = {}
    error = []
    #for each row in header
    for lines in header['header']:
        effective = lines.split('##')[1]
        #if row is fill and need to be translate in dict
        if effective.endswith('>') and not effective.startswith("contig"):
            #print(effective)
            desc = effective.split("=", 1)
            #if key id is not already in dict create key
            if not desc[0] in dico.keys():
                dico[desc[0]] = {}
            #for each value in field ID, maximum split 3 to avoid split in description field (description should be the last) TODO
            tmp = {}
            for it in re.search('<(.*)>', desc[1]).group(1).split(',', 3):
                field = it.split('=')
                key = field[0]
                value = field[1]
                tmp[key] = value
                #Check if Description field is correct
                if key == 'Description':
                    if not value.startswith('"') or not value.endswith('"'):
                        error.append(lines)
            #STACK ID value as dict name and other values key pair in level -1
            if 'ID' in tmp.keys():
                value = tmp['ID']
                tmp.pop('ID')
                dico[desc[0]][value] = tmp
            else:
                dico[desc[0]] = tmp
        #extra field values 
        elif not effective.startswith("contig"):
            val = effective.split('=')
            dico[val[0]] = val[1]
    #print(json.dumps(dico, indent=4))
    return dico, error

def parse_sample_field(dfVar):
    #############
    ### Parsing Sample Field in VCF
    #############

    dico = []
    # dfTest = pd.Series(dfVar.TEM195660.values,index=dfVar.FORMAT).to_dict()
    bad_annotation = []
    sample_list = []

    # Parsing FORMAT field in VCF
    #print("[#INFO] Parsing FORMAT field")
    isample = list(dfVar.columns).index("FORMAT") + 1
    # index: line where the caller identify an event somethings
    for col in dfVar.columns[isample:]:
        #print("#[INFO] " + col + "\n")
        sample_list.append(col)
        for i, row in dfVar.iterrows():
            # print_progress_bar(i, len(dfVar.index)-1)
            # print('\n')
            # print(row['FORMAT'].split(':'), row['bwamem.VarScan_HUSTUMSOL.howard'].split(':'))
            if len(row["FORMAT"].split(":")) != len(row[col].split(":")):
                bad_annotation.append(pd.Series(row[:], index=dfVar.columns))
                continue
            else:
                toadd = pd.Series(
                    row[col].split(":"), index=row["FORMAT"].split(":")
                ).to_dict()
                # toadd.update({"caller":col})
                dico.append(toadd)

    dfSample = pd.DataFrame(dico)
    df_bad_anno = pd.DataFrame(bad_annotation)
    df_final = dfVar.join(dfSample, how="inner")
    df_final.drop(columns=sample_list, inplace=True)
    return df_final, df_bad_anno, sample_list


def parse_info_field(dfVar):
    """
    input: take a dataframe (from vcf)

    output: return a dataframe of the vcf when the info field is parsed
    """

    ############
    # Parsing INFO field from dfVar dataframe containing all informations from vcf
    ############

    # print(dfVar.head())
    infoList = []
    dicoInfo = []
    headers = []

    #print("#[INFO] Parsing INFO field")
    for i, elems in dfVar.iterrows():
        # print_progress_bar(i, len(dfVar.index)-1)
        infoList.append([x.split("=") for x in elems["INFO"].split(";")])

    #print("\n")

    [headers.append(elems[0]) for ite in infoList for elems in ite]
    dfInfo = pd.DataFrame(columns=np.unique(np.array(headers)))
    #print(np.unique(np.array(headers)))

    #print("#[INFO] From INFO field to Dataframe")
    for j, elems in enumerate(infoList):
        # print_progress_bar(j, len(infoList)-1)
        add = {}
        for fields in elems:
            if len(fields) <= 1:
                f = {fields[0]: "TRUE"}
                add.update(f)
            else:
                f = dict([fields])
                add.update(f)

        dicoInfo.append(add)

    #print("\n")
    # print(dicoInfo.keys())
    # print(dict(list(dicoInfo.items())[0:2]))

    df_final = pd.DataFrame(dicoInfo, columns=np.unique(np.array(headers)))

    dfInfo = dfVar.iloc[:, :7].join(df_final, how="inner")
    # Drop old columns info
    # dfInfo.drop(columns="INFO", inplace=True)
    df = dfInfo.join(dfVar.iloc[:, 8:], how="inner")
    return df


def var_to_dataframe(vcf, skiprows, columns):
    """
    from tsv file to pandas dataframe, skiprows number of row to reach header, columns: col which need change (from , to .) to allowed excel filter
    """
    if skiprows:
        df = pd.read_csv(
            vcf,
            skiprows=skiprows,
            sep="\t",
            header=0,
            # chunksize=10000,
            # low_memory=False,
        )
    else:
        df = pd.read_csv(
            vcf,
            sep="\t",
            header=0,
            # chunksize=10000,
            # low_memory=False,
        )
    # if "index" in df:
    #    df = df.drop(columns="index")
    if columns:
        for col in columns:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: x.replace(",", "."))
                df[col] = df[col].astype(float)
    return df


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
                # data["fields"] = lines.strip().split("\t")
                data["fields"] = lines.strip()
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


def cast_config(config):
    conf = read_json(config)
    result = []
    for key, values in conf["variants"]["cast"].items():
        result.append(key + " " + values)
    return ", ".join(result)


def win_to_unix(input, output):
    return systemcall('awk \'{ sub("\r$", ""); print }\' ' + input + " > " + output)


def df_to_vcflike(df, samplename):
    filter_col_pos = df.columns.get_loc("FILTER")
    format_col_pos = df.columns.get_loc("FORMAT")
    dfdone = df.iloc[:, 0 : filter_col_pos + 1]
    # Recreate INFO field from dataframe resulting from db query
    final = []
    l = []
    df_tmp = df.iloc[:, filter_col_pos + 1 : format_col_pos]
    for col in df_tmp.columns:
        tmp = []
        for row in df[col]:
            # if variant does not carry the annotation stack it at '.'
            if row == None:
                row = "."
            tmp.append(str(col) + "=" + str(row))
        l.append(tmp)
    for t in zip_longest(*l):
        final.append(";".join(t))
    dfdone["INFO"] = final
    dfdone["FORMAT"] = df["FORMAT"]
    dfdone[samplename] = df.iloc[:, format_col_pos + 1 :].apply(
        lambda x: ":".join(x.astype(str)), axis=1
    )
    return dfdone


def create_vcf(df, header_dict, output):
    print("#[INFO] Generate " + output)
    with open(output, "w+") as f:
        for row in header_dict["header"]:
            f.write(row + "\n")
        # f.write(header_dict["fields"] + "\n")
    print(df)
    df.to_csv(
        output,
        mode="a",
        doublequote=False,
        escapechar="\\",
        index=False,
        header=True,
        sep="\t",
    )
    return output
