import sys
from collections import defaultdict

def main(args):
    og_tsv = args.og_tsv
    genomeid_file = args.genomeid_file
    offsprings_lvl_file = args.offsprings_lvl_file
    outfile = args.outfile
    cnv = args.cnv
    bi = args.bi
    suffix = "PA.table"
    format = "Presence absence"
    if cnv:
        suffix = "CN.table"
        format = "Copy number"
    if bi:
        cnv = False # suppress cnv format
        suffix = "BiPA.table"
        format = "0/1 Presence absence"

    print("[STEP.1] Read genome accession list.")
    genomeid_list = []
    with open(genomeid_file) as f:
        for line in f:
            genomeid_list.append(line.strip())

    print("[STEP.2] Read offsprings level list.")
    offsprings_lvl = []
    if offsprings_lvl_file:
        with open(offsprings_lvl_file) as f:
            for line in f:
                offsprings_lvl.append(line.strip())
    else:
        print("No offsprings level file. Skip.")

    print("[STEP.3] Get og dict.")
    og_dict = defaultdict(dict)
    with open(og_tsv) as f:
        for line in f:
            # uid, geneId, og_path, status = line.strip().split("\t")
            # geneId, og_path = line.strip().split("\t")
            uid, geneId, og_path = line.strip().split("\t")
            genomeId = geneId.rsplit("_", 1)[0]
            # genomeId must be in genomeid_list
            if genomeId not in genomeid_list:
                continue
            if offsprings_lvl:
                # the first valid orthoTag equal to or lower than this node
                og_list = [ogtag for ogtag in og_path.split("/") if (not ogtag.startswith("-at")) and (ogtag.split("at")[-1] in offsprings_lvl)]
            else:
                # the first valid orthoTag
                og_list = [ogtag for ogtag in og_path.split("/") if not ogtag.startswith("-at")]
            if og_list:
                og_tag = og_list[0] # first orthoTag
                if genomeId in og_dict[og_tag]:
                    og_dict[og_tag][genomeId].append(geneId)
                else:
                    og_dict[og_tag][genomeId] = [geneId]

    print("[STEP.4] Printout", format, "format.")
    g = open(outfile+"."+suffix, "w")
    if cnv:
        g.write("Orthogroup\t"+"\t".join(genomeid_list)+"\n")
        for og_tag in og_dict:
            gene_list = []
            for acc in genomeid_list:
                if acc in og_dict[og_tag]:
                    gene_list.append(str(len(og_dict[og_tag][acc])))
                else:
                    gene_list.append("0")
            g.write("\t".join([og_tag]+gene_list)+"\n")
    elif bi:
        g.write("Orthogroup\t"+"\t".join(genomeid_list)+"\n")
        for og_tag in og_dict:
            gene_list = []
            for acc in genomeid_list:
                if acc in og_dict[og_tag]:
                    gene_list.append("1")
                else:
                    gene_list.append("0")
            g.write("\t".join([og_tag]+gene_list)+"\n")
    else:
        g.write("Orthogroup\t"+"\t".join(genomeid_list)+"\n")
        for og_tag in og_dict:
            gene_list = []
            for acc in genomeid_list:
                if acc in og_dict[og_tag]:
                    gene_list.append(",".join(og_dict[og_tag][acc]))
                else:
                    gene_list.append("")
            g.write("\t".join([og_tag]+gene_list)+"\n")
    g.close()

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Create Presence-Absence/Copy Number table for given genome list and orthologer tsv file.'
    )
    parser.add_argument('-i',
                        help='input orthopath file (including "-atINT123")',
                        required=True,
                        dest='og_tsv',
                        metavar='file_path')
    parser.add_argument('-g',
                        help='genome id file',
                        type=str,
                        required=True,
                        dest='genomeid_file',
                        metavar='file_path')
    parser.add_argument('-n',
                        help='file stored the ranks equal to and lower than internal node of interest',
                        type=str,
                        required=False,
                        dest='offsprings_lvl_file',
                        metavar='file_path')
    parser.add_argument('-o',
                        help='output file prefix',
                        type=str,
                        required=True,
                        dest='outfile',
                        metavar='file_path')
    parser.add_argument('-c',
                        help='output copy number instead of gene list in the table.',
                        action='store_true',
                        dest='cnv')
    parser.add_argument('-b',
                        help='output binary presence/absence instead of gene list in the table.\nIt will suppress -c option if given both.',
                        action='store_true',
                        dest='bi')

    args = parser.parse_args()

    main(args)
