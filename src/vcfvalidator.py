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
    Launch,uncompress_file
)
from validator.header import Checkheader
import json
import pandas as pd
import gzip
from dbvar.database import Databasevar as dbv
from validator.parseargs import Parseoptions
from validator.variants import VCFpreprocess, Checkvariants


def main():
    opt = Parseoptions()
    # return options parsed in dico format and argparse.Namespace respectively
    dico_args, args = opt.get_args()
    
    if args.vcf.endswith('.gz'):
        uncompress_file(args.vcf, '.'.join(args.vcf.split('.')[:-1]))
        args.vcf = '.'.join(args.vcf.split('.')[:-1])
        
    #read config
    #print(args.vcf)
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
        main_correct(variants, header, args, output, config)
    elif args.command == "Scan":
        main_scan(variants, header, args, config)
    elif args.command == "Annotate":
        main_annotate(variants, header, args, output, config)
    else:
        main_test(header, config)


def main_annotate(variants, header, args, dico_args, output, config):
    #could return df if change have been done
 
    #header modified
    if dico_args:
        headercheck = Checkheader(header, dico_args, args.config)
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
        tablename = os.path.basename(args.vcf).split('.', 1)[0]

    db = dbv(
        variants,
        dbname,
        os.path.basename(args.vcf),
        tablename,
        config,
        None,
        args.database
    )

    if not db.table_exists():
        #FROM dataframe to sql db
        db.create_table()

    create_vcf(variants, hprocess, output)

def main_scan(variants, header, args, config):
    print("Check integrity of "+os.path.abspath(args.vcf)+" ...\n")
    errors = is_utf8(args.vcf)
    if errors:
        for err in errors:
            print(err)
    else:
        print("#[INFO] ASCII Check OK")
    Checkvariants(variants, "sample").check_col()
    Checkheader(
            header, {}, config
        ).header_check()
    xh = explode_header(header)

    if args.dbname:
        dbname = args.dbname
    else:
        dbname = os.path.join("dbvar", os.path.basename(args.vcf).split(".")[0] + ".db")
    
    if args.tablename:
        tablename = args.tablename
    else:
        tablename = os.path.basename(args.vcf).split('.', 1)[0]

    db = dbv(
        variants,
        dbname,
        os.path.basename(args.vcf),
        tablename,
        config,
        xh,
        args.database
    )
    if not db.table_exists():
        #FROM dataframe to sql db
        db.create_table()

    db.chromosome_check()
    chann = db.check_annotations()
    if chann:
        for items in chann:
            print("WARNING missing "+items+" in VCF header")
    elif not chann and "INFO" in xh:
        print("#[INFO] all variants annotations are linked with header annotation")
    
    with open('json_data.json', 'w+') as outfile:
        json.dump(xh, outfile, indent=4)


def main_correct(variants, header, args, output, config):
    print("Correct "+os.path.abspath(args.vcf)+" ...\n")
    #only if correct mode activate #TODO
    #print("Check integrity of "+args.vcf+" ...")
    # add correct column
    # header with basic GT AD and DP add
    cvar = Checkvariants(variants, "sample")
    
    #return True if basic columns where add false otherwise
    variants_new = cvar.add_basic()

    if variants_new:
        Checkheader(
            header, config["variants"]["basic"], config
        ).correct_header()

    # worked
    if args.dbname:
        dbname = args.dbname
    else:
        dbname = os.path.join("dbvar", os.path.basename(args.vcf).split(".")[0] + ".db")
    
    if args.tablename:
        tablename = args.tablename
    else:
        tablename = os.path.basename(args.vcf).split('.', 1)[0]

    db = dbv(
        variants,
        dbname,
        os.path.basename(args.vcf),
        tablename,
        config,
        None,
        args.database
    )

    if not db.table_exists():
        #FROM dataframe to sql db
        db.create_table()

    #chromosome check only if correct mode activate
    res_chr = db.chromosome_check()
    if res_chr:
        print("#[INFO] correct CHR field")
        for values in res_chr:
            db.update_value(values)

    #Generate vcf only if correction mode is enable
    #from db sql to dataframe after managing  columns contains INFO and format split so more than 10 for multisample
    dfwrite_tmp = pd.DataFrame(db.db_to_dataframe(), columns=db.get_col_name())
    dfwrite = df_to_vcflike(dfwrite_tmp, "Likely_benin")

    print(dfwrite)
    #Regenerate vcf final 
    #create_vcf(dfwrite, header, output)

#def main_analyze(variants, header, args, output):
#    '''
#    check many possible malformation in vcf
#    '''
#    pass

def main_test(header, config):
    #print(explode_header(header))
    pass


if __name__ == "__main__":
    fancystdout("speed", "VCFvalidator")
    Launch().__str__()
    main()
