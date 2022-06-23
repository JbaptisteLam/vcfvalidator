# aim: validator de vcf avec vcf en db sqlite pour les variants


import os
from os.path import join as osj
from utils.utils import (
    read_json,
    df_to_vcflike,
    create_vcf,
    explode_header,
    is_utf8,
    fancystdout,
    Launch,
    uncompress_file,
)
from validator.header import Checkheader
import json
import pandas as pd
import gzip
from dbvar.database import Databasevar as dbv
from dbvar.database import ReadDB as rdb
from validator.parseargs import Parseoptions
from validator.variants import VCFpreprocess, Checkvariants


def main():
    opt = Parseoptions()
    # return options parsed in dico format and argparse.Namespace respectively
    dico_args, args = opt.get_args()

    # If vcf as input
    if args.vcf.endswith(".gz"):
        uncompress_file(args.vcf, ".".join(args.vcf.split(".")[:-1]))
        args.vcf = ".".join(args.vcf.split(".")[:-1])

    # read config
    # print(args.vcf)
    config = read_json(args.config)

    # if no output specified generate vcf as input file name with correct
    if args.output == "basic":
        output = osj(
            os.path.dirname(args.vcf),
            os.path.basename(args.vcf).split(".")[0] + "_correct.vcf",
        )
    else:
        output = args.output
    if os.path.exists(output):
        print("ERROR " + output + " exist remove or choose an other output filename")
        exit()

    # Load vcf file
    if args.command != "Database":
        variants, header = VCFpreprocess(args.vcf).get_values()
    else:
        variants = None
        header = None

    # choice main func
    if args.command == "Correct":
        main_correct(variants, header, args, output, config, dico_args)
    elif args.command == "Scan":
        main_scan(variants, header, args, config)
    elif args.command == "Annotate":
        main_annotate(variants, header, args, dico_args, output, config)
        # No need inpu
    elif args.command == "Database":
        main_database(variants, header, args, config)
    else:
        main_test(header, config)


def main_annotate(variants, header, args, dico_args, output, config):
    # could return df if change have been done
    # TODO test
    variants = variants.iloc[:5000]

    # header modified
    if dico_args:
        print("#[INFO] Actions ", dico_args)
        headercheck = Checkheader(header, dico_args, config)
        hprocess = headercheck.correct_header()

    else:
        hprocess = header

    if args.dbname:
        dbname = args.dbname
    else:
        dbname = os.path.join("dbvar", os.path.basename(args.vcf).split(".")[0] + ".db")

    if args.tablename:
        tablename = args.tablename
    else:
        tablename = os.path.basename(args.vcf).split(".", 1)[0]

    xh, _chrom, _errs = explode_header(header)
    db = dbv(
        variants,
        dbname,
        os.path.basename(args.vcf),
        tablename,
        config,
        [xh, _chrom],
        args.database,
        dico_args,
    )

    if not db.table_exists():
        # FROM dataframe to sql db
        db.create_table()

    # Act on variant, add, edit or remove
    db.correct_variants()
    dfwrite_tmp = db.db_to_dataframe()
    dfwrite = df_to_vcflike(dfwrite_tmp, "Likely_benin_Clinvar", db.get_col_name())
    create_vcf(dfwrite, hprocess, output)


def main_scan(variants, header, args, config):
    print("Check integrity of " + os.path.abspath(args.vcf) + " ...\n")
    errors = is_utf8(args.vcf)
    if errors:
        for err in errors:
            print(err)
    else:
        print("#[INFO] ASCII Check OK")
    Checkvariants(variants, "sample").check_col()
    Checkheader(header, {}, config).header_check()
    xh, _, _errs = explode_header(header)

    if args.dbname:
        dbname = args.dbname
    else:
        dbname = os.path.join("dbvar", os.path.basename(args.vcf).split(".")[0] + ".db")

    if args.tablename:
        tablename = args.tablename
    else:
        tablename = os.path.basename(args.vcf).split(".", 1)[0]

    db = dbv(
        variants,
        dbname,
        os.path.basename(args.vcf),
        tablename,
        config,
        xh,
        args.database,
        None,
    )
    if not db.table_exists():
        # FROM dataframe to sql db
        db.create_table()

    db.chromosome_check()
    chann = db.check_annotations()
    if chann:
        for items in chann:
            print("WARNING missing " + items + " in VCF header")
    elif not chann and "INFO" in xh:
        print("#[INFO] all variants annotations are linked with header annotation")

    with open("json_data.json", "w+") as outfile:
        json.dump(xh, outfile, indent=4)


def main_correct(variants, header, args, output, config, dico_args):
    print("Correct " + os.path.abspath(args.vcf) + " ...\n")
    # only if correct mode activate #TODO
    # print("Check integrity of "+args.vcf+" ...")
    # add correct column
    # header with basic GT AD and DP add
    cvar = Checkvariants(variants, "sample")

    # return True if basic columns where add false otherwise
    variants_new = cvar.add_basic()

    if variants_new:
        Checkheader(header, config["variants"]["basic"], config).correct_header()

    xh, _, _errs = explode_header(header)
    # worked
    if args.dbname:
        dbname = args.dbname
    else:
        dbname = os.path.join("dbvar", os.path.basename(args.vcf).split(".")[0] + ".db")

    if args.tablename:
        tablename = args.tablename
    else:
        tablename = os.path.basename(args.vcf).split(".", 1)[0]

    db = dbv(
        variants,
        dbname,
        os.path.basename(args.vcf),
        tablename,
        config,
        xh,
        args.database,
        dico_args,
    )

    if not db.table_exists():
        # FROM dataframe to sql db
        db.create_table()

    # chromosome check only if correct mode activate
    res_chr = db.chromosome_check()
    if res_chr:
        print("#[INFO] correct CHR field")
        for values in res_chr:
            db.update_value(values)

    # Generate vcf only if correction mode is enable
    # from db sql to dataframe after managing  columns contains INFO and format split so more than 10 for multisample
    dfwrite_tmp = db.db_to_dataframe()
    dfwrite = df_to_vcflike(dfwrite_tmp, "Likely_benin", db.get_col_name())

    print(dfwrite)
    # Regenerate vcf final
    create_vcf(dfwrite, header, output)


def main_database(variants, header, args, config):
    if args.output == "basic":
        print("ERROR specify output file for vcf")
        exit()
    if args.dbname:
        dbname = args.dbname
    else:
        dbname = os.path.join("dbvar", os.path.basename(args.vcf).split(".")[0] + ".db")

    if args.tablename:
        tablename = args.tablename
    else:
        tablename = os.path.basename(args.vcf).split(".", 1)[0]

    if not args.merge:
        xh, _, _errs = explode_header(header)
        db = dbv(
            variants,
            dbname,
            os.path.basename(args.vcf),
            tablename,
            config,
            xh,
            args.database,
            None,
        )

        if not db.table_exists():
            # FROM dataframe to sql db
            db.create_table()

        header_file = osj(
            os.path.dirname(args.output),
            os.path.basename(args.output).split(".")[0] + "_header.txt",
        )
        if args.save:
            print("#[INFO] Write header in " + header_file)
            with open(header_file, "w+") as f:
                for row in header["header"]:
                    f.write(row + "\n")

    else:
        head_load = {}
        head_tmp = args.merge.split(",", 1)[0]
        with open(head_tmp, "r") as f:
            print(head_tmp)
            head_load["header"] = []
            for row in f:
                row = row.strip()
                head_load["header"].append(row)

        db_load = args.merge.split(",", 1)[1]
        db = rdb(
            db_load,
            tablename,
        )
        df_tmp = df_to_vcflike(db.db_to_dataframe(), args.samplename, db.get_col_name())
        create_vcf(df_tmp, head_load, args.output)


# def main_analyze(variants, header, args, output):
#    '''
#    check many possible malformation in vcf
#    '''
#    pass


def main_test(header, config):
    # print(explode_header(header))
    pass


if __name__ == "__main__":
    fancystdout("speed", "VCFvalidator")
    Launch().__str__()
    main()
