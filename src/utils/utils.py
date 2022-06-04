import pandas as pd
import json
import subprocess
import numpy as np


def parse_sample_field(dfVar):
    #############
    ### Parsing Sample Field in VCF
    #############

    dico = []
    # dfTest = pd.Series(dfVar.TEM195660.values,index=dfVar.FORMAT).to_dict()
    bad_annotation = []
    sample_list = []

    # Parsing FORMAT field in VCF
    print("[#INFO] Parsing FORMAT field")
    isample = list(dfVar.columns).index("FORMAT") + 1
    # index: line where the caller identify an event somethings
    for col in dfVar.columns[isample:]:
        print("#[INFO] " + col + "\n")
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

    print("#[INFO] Parsing INFO field")
    for i, elems in dfVar.iterrows():
        # print_progress_bar(i, len(dfVar.index)-1)
        infoList.append([x.split("=") for x in elems["INFO"].split(";")])

    print("\n")

    [headers.append(elems[0]) for ite in infoList for elems in ite]
    dfInfo = pd.DataFrame(columns=np.unique(np.array(headers)))
    print(np.unique(np.array(headers)))

    print("#[INFO] From INFO field to Dataframe")
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

    print("\n")
    # print(dicoInfo.keys())
    # print(dict(list(dicoInfo.items())[0:2]))

    df_final = pd.DataFrame(dicoInfo, columns=np.unique(np.array(headers)))

    dfInfo = dfVar.iloc[:, 0:6].join(df_final, how="inner")
    print(dfVar.iloc[:, 0:6].columns)
    print(dfInfo.columns)
    # Drop old columns info
    # dfInfo.drop(columns="INFO", inplace=True)
    df = pd.concat([dfInfo, dfVar.iloc[8:]])
    print(df.columns)
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
    print(df)
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
