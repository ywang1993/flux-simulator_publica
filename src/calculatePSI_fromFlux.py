"""
 # @author [Yuanyuan Wang]
 # @email [wyynju1993@gmail.com]
 # @create date 2020-05-28 04:19:59
 # @modify date 2020-05-28 04:19:59
 # @desc [description]
"""

from __future__ import print_function
import os, sys
import argparse
from datetime import datetime
import numpy as np

from detect_AS_from_GTF import *

# GLOBAL VARIABLE
USAGE = """python %(prog)s --gtf path/to/gtf.txt --pro path/to/XX.PRO [options]"""


def get_args():
    parser = argparse.ArgumentParser(usage=USAGE)
    parser.add_argument(
        "--od", action="store", help="output folder, default: working folder", default="./", dest="outpath"
    )
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument(
        "--gtf", action="store", help="file path to GTF file used for flux simulator", dest="fn_gtf", required=True
    )
    requiredNamed.add_argument(
        "--pro",
        action="store",
        help="file path of expression profile used in flux simulator",
        dest="fn_pro",
        required=True,
    )

    args = parser.parse_args()
    return args


def parse_pro(pro, tx2gene):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start: parsing exp profile.")

    count_gene_tx = {}
    with open(pro, "r") as fin:
        for line in fin:
            loc, tx, ctg, tlen, i, count = line.strip().split("\t")[0:6]

            if int(count) == 0:
                continue

            gene = tx2gene[tx]

            if gene not in count_gene_tx:
                count_gene_tx[gene] = {"TX": {}}
            count_gene_tx[gene]["TX"][tx] = int(count)

    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "end  : parsing exp profile.")
    return count_gene_tx


def SE_byTX(SE, gene_struct):
    count = {}

    for key in SE:
        x = SE_info(SE[key])

        if x.GeneID not in gene_struct:
            continue

        count[x.uniqID] = {"inc": [], "skp": [], "gene_name": x.geneSymbol, "gene_id": x.GeneID}

        for tx in gene_struct[x.GeneID]["TX"]:
            # remove one-exon tx
            if len(gene_struct[x.GeneID]["TX"][tx]) == 1:
                continue
            # initiation
            left = -1
            middle = -1
            right = -1
            # locate SE exons
            for idx in range(len(gene_struct[x.GeneID]["TX"][tx])):
                exon = gene_struct[x.GeneID]["TX"][tx][idx]
                if exon[1] == x.upstreamEE:
                    left = idx
                elif exon[0] == (x.exonStart_0base + 1) and exon[1] == x.exonEnd:
                    middle = idx
                elif exon[0] == (x.downstreamES + 1):
                    right = idx
            # determine inc or skp
            if not (left == -1 or right == -1):
                if middle == -1 and (right - left) == 1:
                    count[x.uniqID]["skp"].append(tx)
                elif (right - middle) == 1 and (middle - left) == 1:
                    count[x.uniqID]["inc"].append(tx)
    return count


def RI_byTX(RI, gene_struct):
    count = {}
    for key in RI:
        x = RI_info(RI[key])

        if x.GeneID not in gene_struct:
            continue

        count[x.uniqID] = {"inc": [], "skp": [], "gene_name": x.geneSymbol, "gene_id": x.GeneID}

        for tx in gene_struct[x.GeneID]["TX"]:
            # initiation
            left = -1
            right = -1
            include = False
            # locate RI exons
            for idx in range(len(gene_struct[x.GeneID]["TX"][tx])):
                exon = gene_struct[x.GeneID]["TX"][tx][idx]
                if exon[1] == x.upstreamEE:
                    left = idx
                elif exon[0] == (x.downstreamES + 1):
                    right = idx
                elif exon[0] <= x.upstreamEE and exon[1] >= (x.downstreamES + 1):
                    include = True
            # determine inc or skp
            if (left > -1) and (right > -1) and (right - left) == 1:
                count[x.uniqID]["skp"].append(tx)
            if include:
                count[x.uniqID]["inc"].append(tx)
    return count


def MXE_byTX(MXE, gene_struct):
    count = {}

    for key in MXE:
        x = MXE_info(MXE[key])

        if x.GeneID not in gene_struct:
            continue

        count[x.uniqID] = {"inc": [], "skp": [], "gene_name": x.geneSymbol, "gene_id": x.GeneID}

        for tx in gene_struct[x.GeneID]["TX"]:
            # initiation
            left = -1
            texon = -1
            aexon = -1
            right = -1
            # locate MXE
            for idx in range(len(gene_struct[x.GeneID]["TX"][tx])):
                exon = gene_struct[x.GeneID]["TX"][tx][idx]
                if exon[1] == x.upstreamEE:
                    left = idx
                elif exon[0] == (x.downstreamES + 1):
                    right = idx
                elif (exon[0] == x.stExonStart_0base + 1) and (exon[1] == x.stExonEnd):
                    texon = idx
                elif (exon[0] == x.ndExonStart_0base + 1) and (exon[1] == x.ndExonEnd):
                    aexon = idx
            # determine inc or skp
            if (left == -1) or (right == -1):
                continue
            if (texon > -1) and (right - texon == 1) and (texon - left == 1):
                count[x.uniqID]["inc"].append(tx)
            elif (aexon > -1) and (right - aexon == 1) and (aexon - left == 1):
                count[x.uniqID]["skp"].append(tx)

        if x.strand == "-":
            count[x.uniqID]["inc"], count[x.uniqID]["skp"] = count[x.uniqID]["skp"], count[x.uniqID]["inc"]

    return count


def AXSS_byTX(AXSS, gene_struct):
    count = {}

    for key in AXSS:
        x = AXSS_info(AXSS[key])

        if x.GeneID not in gene_struct:
            continue

        count[x.uniqID] = {"inc": [], "skp": [], "gene_name": x.geneSymbol, "gene_id": x.GeneID}

        for tx in gene_struct[x.GeneID]["TX"]:
            # remove one-exon tx
            if len(gene_struct[x.GeneID]["TX"][tx]) == 1:
                continue
            # initiation
            flank = -1
            elong = -1
            eshort = -1
            # locate AXSS exons
            for idx in range(len(gene_struct[x.GeneID]["TX"][tx])):
                exon = gene_struct[x.GeneID]["TX"][tx][idx]
                if x.flankingEE <= x.longExonStart_0base:
                    if exon[1] == x.flankingEE:
                        flank = idx
                    elif exon[0] == (x.longExonStart_0base + 1):
                        elong = idx
                    elif exon[0] == (x.shortES + 1):
                        eshort = idx
                elif x.flankingES >= x.longExonEnd:
                    if exon[0] == (x.flankingES + 1):
                        flank = idx
                    elif exon[1] == x.longExonEnd:
                        elong = idx
                    elif exon[1] == x.shortES:
                        eshort = idx
                else:
                    sys.exit("Error: AXSS event error")

            # determine inc or skp
            if x.flankingEE <= x.longExonStart_0base:
                if (flank > -1) and (elong > -1) and (elong - flank == 1):
                    count[x.uniqID]["inc"].append(tx)
                elif (flank > -1) and (eshort > -1) and (eshort - flank == 1):
                    count[x.uniqID]["skp"].append(tx)
            elif x.flankingES >= x.longExonEnd:
                if (flank > -1) and (elong > -1) and (elong - flank == -1):
                    count[x.uniqID]["inc"].append(tx)
                elif (flank > -1) and (eshort > -1) and (eshort - flank == -11):
                    count[x.uniqID]["skp"].append(tx)
    return count


def AS_byTX(gene_struct, SE, MXE, RI, A5SS, A3SS):
    se_tx = SE_byTX(SE, gene_struct)
    mxe_tx = MXE_byTX(MXE, gene_struct)
    ri_tx = RI_byTX(RI, gene_struct)
    a5ss_tx = AXSS_byTX(A5SS, gene_struct)
    a3ss_tx = AXSS_byTX(A3SS, gene_struct)
    return se_tx, mxe_tx, ri_tx, a5ss_tx, a3ss_tx


def CalculatePSI(AS, count_gene_tx):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start: calculating PSI")
    PSI = {}
    total_count = {}
    for key in AS:

        inc_count = 0
        skp_count = 0
        gene = AS[key]["gene_id"]

        if gene not in count_gene_tx:
            PSI[key] = np.nan
            total_count[key] = 0
            continue

        if len(AS[key]["inc"]) == 0 and len(AS[key]["skp"]) == 0:
            PSI[key] = np.nan
            total_count[key] = 0
            continue

        if len(AS[key]["inc"]) != 0:
            for tx in AS[key]["inc"]:
                if tx in count_gene_tx[gene]["TX"]:
                    inc_count += count_gene_tx[gene]["TX"][tx]
        if len(AS[key]["skp"]) != 0:
            for tx in AS[key]["skp"]:
                if tx in count_gene_tx[gene]["TX"]:
                    skp_count += count_gene_tx[gene]["TX"][tx]
        if inc_count + skp_count == 0:
            PSI[key] = np.nan
        else:
            PSI[key] = float(inc_count) / (inc_count + skp_count)
        total_count[key] = inc_count + skp_count
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "end  : calculating PSI")
    return total_count, PSI


def save(event, as_tx, basic_AS_keys, count_gene_tx, od):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start  : saving for ", event)

    line = "{}\t{}\t{}\t{}\t{}\n"
    total_count, PSI = CalculatePSI(as_tx, count_gene_tx)
    with open(os.path.join(os.path.abspath(od), "PSI_" + event + ".txt"), "w") as fo:
        fo.write(line.format("ID", "PSI", "count", "GeneID", "geneSymbol"))
        for key in PSI:
            fo.write(
                line.format(key, str(PSI[key]), str(total_count[key]), as_tx[key]["gene_id"], as_tx[key]["gene_name"])
            )
    with open(os.path.join(os.path.abspath(od), "PSI_basic_" + event + ".txt"), "w") as fo:
        fo.write(line.format("ID", "PSI", "count", "GeneID", "geneSymbol"))
        for key in basic_AS_keys:
            fo.write(
                line.format(key, str(PSI[key]), str(total_count[key]), as_tx[key]["gene_id"], as_tx[key]["gene_name"])
            )


def main():
    # fn_pro = "../profiles/count_modified.PRO"
    # fn_gtf = "../reference/gencode.v31lift37.annotation.gtf"
    args = get_args()

    gene_struct, tx2gene = parse_gtf(args.fn_gtf)
    gene_struct = build_graph(gene_struct)
    SE, MXE, RI, A5SS, A3SS = detect_AS(gene_struct)
    basic_SE_keys, basic_MXE_keys, basic_RI_keys, basic_A5SS_keys, basic_A3SS_keys = get_basic_AS(
        gene_struct, SE, MXE, RI, A5SS, A3SS
    )

    count_gene_tx = parse_pro(args.fn_pro, tx2gene)
    se_tx, mxe_tx, ri_tx, a5ss_tx, a3ss_tx = AS_byTX(gene_struct, SE, MXE, RI, A5SS, A3SS)
    # total_count, PSI = CalculatePSI(se_tx, count_gene_tx)

    save("SE", se_tx, basic_SE_keys, count_gene_tx, args.outpath)
    save("MXE", mxe_tx, basic_MXE_keys, count_gene_tx, args.outpath)
    save("RI", ri_tx, basic_RI_keys, count_gene_tx, args.outpath)
    save("A5SS", a5ss_tx, basic_A5SS_keys, count_gene_tx, args.outpath)
    save("A3SS", a3ss_tx, basic_A3SS_keys, count_gene_tx, args.outpath)


if __name__ == "__main__":
    main()
