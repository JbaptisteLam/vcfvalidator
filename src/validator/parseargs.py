import argparse


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
    parser_annotate = subparsers.add_parser(
        "Annotate",
        help="Add, edit or remove row in VCF header",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
    parser_correct = subparsers.add_parser(
        "Correct",
        help="Correct VCF header regarding variants annotations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    args = parser.parse_args()
    return args
