"""
 # @author [Yuanyuan Wang]
 # @email [wyynju1993@gmail.com]
 # @create date 2020-05-28 01:22:52
 # @modify date 2020-05-28 04:19:50
 # @desc [description]
"""
from __future__ import print_function
import os, sys
from datetime import datetime


class GTF(object):
    def __init__(self, line):
        self.line = line
        (
            self.seqname,
            self.source,
            self.feature,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attribute,
        ) = line.strip().split("\t")
        self.start = int(self.start)
        self.end = int(self.end)
        self.chrom = self.seqname if "chr" in self.seqname else "chr" + self.seqname
        self.gene_id = self.attribute.split("gene_id")[1].split('"')[1]
        self.gene_name = self.attribute.split("gene_name")[1].split('"')[1]
        self.gene_type = self.attribute.split("gene_type")[1].split('"')[1]
        self.transcript_id = (
            self.attribute.split("transcript_id")[1].split('"')[1] if "transcript_id" in self.attribute else None
        )
        self.exon_number = (
            int(self.attribute.split("exon_number")[1].split(";")[0].strip(""))
            if "exon_number" in self.attribute
            else None
        )

    def __str__(self):
        return self.line


class SE_info(object):
    def __init__(self, info):
        self.info = info
        (
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.exonStart_0base,
            self.exonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
        ) = self.info
        self.uniqID = "|".join(
            [
                self.chrom + ":" + str(self.exonStart_0base) + "-" + str(self.exonEnd),
                self.strand,
                str(self.upstreamEE),
                str(self.downstreamES),
            ]
        )
        return

    def __str__(self):
        return str(self.info)


class RI_info(object):
    def __init__(self, info):
        self.info = info
        (
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.riExonStart_0base,
            self.riExonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
        ) = self.info
        self.uniqID = "|".join(
            [
                self.chrom + ":" + str(self.upstreamEE) + "-" + str(self.downstreamES),
                self.strand,
                str(self.upstreamES),
                str(self.downstreamEE),
            ]
        )
        return

    def __str__(self):
        return str(self.info)


class AXSS_info(object):
    def __init__(self, info):
        self.info = info
        (
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.longExonStart_0base,
            self.longExonEnd,
            self.shortES,
            self.shortEE,
            self.flankingES,
            self.flankingEE,
        ) = self.info

        if int(self.flankingEE) <= int(self.longExonStart_0base):
            self.uniqID = "|".join(
                [
                    self.chrom + ":" + str(self.flankingES) + "-" + str(self.flankingEE),
                    self.strand,
                    str(self.longExonStart_0base),
                    str(self.shortES),
                ]
            )
        elif int(self.flankingES) >= int(self.longExonEnd):
            self.uniqID = "|".join(
                [
                    self.chrom + ":" + str(self.flankingES) + "-" + str(self.flankingEE),
                    self.strand,
                    str(self.longExonEnd),
                    str(self.shortEE),
                ]
            )
        else:
            sys.exit("Error: check A5SS and A3SS file")
        return

    def __str__(self):
        return str(self.info)


class MXE_info(object):
    def __init__(self, info):
        self.info = info
        (
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.stExonStart_0base,
            self.stExonEnd,
            self.ndExonStart_0base,
            self.ndExonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
        ) = self.info

        self.uniqID = "|".join(
            [
                self.chrom
                + ":"
                + str(self.stExonStart_0base)
                + "-"
                + str(self.stExonEnd)
                + ":"
                + str(self.ndExonStart_0base)
                + "-"
                + str(self.ndExonEnd),
                self.strand,
                str(self.upstreamEE),
                str(self.downstreamES),
            ]
        )
        return

    def __str__(self):
        return str(self.info)


def parse_gtf(fn_gtf):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start: parsing GTF profile.")

    gene_struct = {}
    tx2gene = {}
    with open(fn_gtf, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            x = GTF(line)

            if x.feature == "exon":
                if x.gene_id not in gene_struct:
                    gene_struct[x.gene_id] = {
                        "TX": {x.transcript_id: [(x.start, x.end)]},
                        "exon": {(x.start, x.end)},
                        "chrom": x.chrom,
                        "strand": x.strand,
                        "gene_name": x.gene_name,
                    }
                elif x.transcript_id not in gene_struct[x.gene_id]["TX"]:
                    gene_struct[x.gene_id]["TX"][x.transcript_id] = [(x.start, x.end)]
                    gene_struct[x.gene_id]["exon"].add((x.start, x.end))
                else:
                    gene_struct[x.gene_id]["TX"][x.transcript_id].append((x.start, x.end))
                    gene_struct[x.gene_id]["exon"].add((x.start, x.end))
            elif x.feature == "transcript" and x.transcript_id not in tx2gene:
                tx2gene[x.transcript_id] = x.gene_id

    for gene in gene_struct:
        gene_struct[gene]["idx_exon"] = sorted(list(gene_struct[gene]["exon"]))
        gene_struct[gene]["exon_idx"] = {j: i for i, j in enumerate(gene_struct[gene]["idx_exon"])}
        del gene_struct[gene]["exon"]
        for tx in gene_struct[gene]["TX"]:
            gene_struct[gene]["TX"][tx].sort()

    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "end  : parsing GTF profile.")
    return gene_struct, tx2gene


def build_graph(gene_struct):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start: building splicing graph.")

    for gene in gene_struct:
        # initiation of sg
        num_exon = len(gene_struct[gene]["exon_idx"])
        gene_struct[gene]["sg"] = [[set() for j in range(num_exon)] for i in range(num_exon)]

        # build sg
        for tx in gene_struct[gene]["TX"]:
            num_exon = len(gene_struct[gene]["TX"][tx])
            if num_exon >= 2:
                left = gene_struct[gene]["exon_idx"][gene_struct[gene]["TX"][tx][0]]
                right = gene_struct[gene]["exon_idx"][gene_struct[gene]["TX"][tx][1]]
                gene_struct[gene]["sg"][left][right].add(left)
                gene_struct[gene]["sg"][right][left].add(right)
                if num_exon > 2:
                    for i in range(num_exon - 2):
                        left = gene_struct[gene]["exon_idx"][gene_struct[gene]["TX"][tx][i]]
                        mid = gene_struct[gene]["exon_idx"][gene_struct[gene]["TX"][tx][i + 1]]
                        right = gene_struct[gene]["exon_idx"][gene_struct[gene]["TX"][tx][i + 2]]
                        gene_struct[gene]["sg"][mid][right].add(left)
                        gene_struct[gene]["sg"][mid][left].add(right)

    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "end  : building splicing graph.")
    return gene_struct


def check_edge(gene_struct, gene, left, right):
    num_edge = len(gene_struct[gene]["sg"][left][right])
    if num_edge == 0:
        return False
    else:
        return num_edge


def detect_SE_MXE(gene_struct):
    # SE
    SE = {}
    MXE = {}
    for gene in gene_struct:
        num_exon = len(gene_struct[gene]["exon_idx"])
        for idx in range(num_exon):
            chrom = gene_struct[gene]["chrom"]
            strand = gene_struct[gene]["strand"]
            gene_name = gene_struct[gene]["gene_name"]
            exon = gene_struct[gene]["idx_exon"][idx]
            # if check_exon_overlap(gene_struct, gene, exon):
            #     continue
            for left in range(idx):
                if not check_edge(gene_struct, gene, left, idx):
                    continue
                for right in range(idx + 1, num_exon):
                    if not check_edge(gene_struct, gene, idx, right):
                        continue

                    se_key = (
                        chrom,
                        gene_struct[gene]["idx_exon"][left][1],
                        gene_struct[gene]["idx_exon"][right][0],
                        gene_struct[gene]["idx_exon"][idx][0],
                        gene_struct[gene]["idx_exon"][idx][1],
                    )

                    mxe_key = (
                        chrom,
                        gene_struct[gene]["idx_exon"][left][1],
                        gene_struct[gene]["idx_exon"][right][0],
                        gene_struct[gene]["idx_exon"][idx][0],
                        gene_struct[gene]["idx_exon"][idx][1],
                    )

                    for skip_left in range(idx):
                        if gene_struct[gene]["idx_exon"][skip_left][1] != gene_struct[gene]["idx_exon"][left][1]:
                            continue
                        for skip_right in range(idx + 1, num_exon):
                            if gene_struct[gene]["idx_exon"][skip_right][0] != gene_struct[gene]["idx_exon"][right][0]:
                                continue
                            if not check_edge(gene_struct, gene, skip_left, skip_right):
                                continue

                            if se_key not in SE:
                                SE[se_key] = (
                                    gene,
                                    gene_name,
                                    chrom,
                                    strand,
                                    gene_struct[gene]["idx_exon"][idx][0] - 1,
                                    gene_struct[gene]["idx_exon"][idx][1],
                                    gene_struct[gene]["idx_exon"][left][0] - 1,
                                    gene_struct[gene]["idx_exon"][left][1],
                                    gene_struct[gene]["idx_exon"][right][0] - 1,
                                    gene_struct[gene]["idx_exon"][right][1],
                                )
                            break
                        else:
                            continue
                        break

                    # MXE
                    for mid in range(idx + 1, num_exon):
                        if gene_struct[gene]["idx_exon"][mid][0] <= exon[1]:
                            continue
                        if not check_edge(gene_struct, gene, left, mid):
                            continue
                        if not check_edge(gene_struct, gene, mid, right):
                            continue

                        mxe_key = tuple(
                            list(mxe_key)
                            + [gene_struct[gene]["idx_exon"][mid][0], gene_struct[gene]["idx_exon"][mid][1],]
                        )

                        if mxe_key not in MXE:
                            MXE[mxe_key] = (
                                gene,
                                gene_name,
                                chrom,
                                strand,
                                gene_struct[gene]["idx_exon"][idx][0] - 1,
                                gene_struct[gene]["idx_exon"][idx][1],
                                gene_struct[gene]["idx_exon"][mid][0] - 1,
                                gene_struct[gene]["idx_exon"][mid][1],
                                gene_struct[gene]["idx_exon"][left][0] - 1,
                                gene_struct[gene]["idx_exon"][left][1],
                                gene_struct[gene]["idx_exon"][right][0] - 1,
                                gene_struct[gene]["idx_exon"][right][1],
                            )

    return SE, MXE


def detect_RI(gene_struct):
    RI = {}
    for gene in gene_struct:
        num_exon = len(gene_struct[gene]["exon_idx"])
        for idx in range(num_exon):
            chrom = gene_struct[gene]["chrom"]
            strand = gene_struct[gene]["strand"]
            gene_name = gene_struct[gene]["gene_name"]
            exon = gene_struct[gene]["idx_exon"][idx]

            for i in range(idx):
                if not check_edge(gene_struct, gene, i, idx):
                    continue
                tmp_pair = (gene_struct[gene]["idx_exon"][i][0], gene_struct[gene]["idx_exon"][idx][1])
                if tmp_pair not in gene_struct[gene]["exon_idx"]:
                    continue

                ri_key = (chrom, gene_struct[gene]["idx_exon"][i][1], gene_struct[gene]["idx_exon"][idx][0])

                if ri_key not in RI:
                    RI[ri_key] = (
                        gene,
                        gene_name,
                        chrom,
                        strand,
                        gene_struct[gene]["idx_exon"][i][0] - 1,
                        gene_struct[gene]["idx_exon"][idx][1],
                        gene_struct[gene]["idx_exon"][i][0] - 1,
                        gene_struct[gene]["idx_exon"][i][1],
                        gene_struct[gene]["idx_exon"][idx][0] - 1,
                        gene_struct[gene]["idx_exon"][idx][1],
                    )
    return RI


def detect_AXSS(gene_struct):
    A5SS = {}
    A3SS = {}
    for gene in gene_struct:
        num_exon = len(gene_struct[gene]["exon_idx"])

        for idx in range(num_exon):
            chrom = gene_struct[gene]["chrom"]
            strand = gene_struct[gene]["strand"]
            gene_name = gene_struct[gene]["gene_name"]
            exon = gene_struct[gene]["idx_exon"][idx]

            for i in range(idx):
                if not check_edge(gene_struct, gene, i, idx):
                    continue
                for j in range(i + 1, idx):
                    if not check_edge(gene_struct, gene, j, idx):
                        continue
                    if gene_struct[gene]["idx_exon"][i][0] != gene_struct[gene]["idx_exon"][j][0]:
                        continue

                    alt_key = (
                        chrom,
                        gene_struct[gene]["idx_exon"][i][1],
                        gene_struct[gene]["idx_exon"][j][1],
                        gene_struct[gene]["idx_exon"][idx][0] - 1,
                    )

                    if gene_struct[gene]["strand"] == "+":
                        if alt_key not in A5SS:
                            A5SS[alt_key] = (
                                gene,
                                gene_name,
                                chrom,
                                strand,
                                gene_struct[gene]["idx_exon"][i][0] - 1,
                                gene_struct[gene]["idx_exon"][i][1],
                                gene_struct[gene]["idx_exon"][j][0] - 1,
                                gene_struct[gene]["idx_exon"][j][1],
                                gene_struct[gene]["idx_exon"][idx][0] - 1,
                                gene_struct[gene]["idx_exon"][idx][1],
                            )
                    elif gene_struct[gene]["strand"] == "-":
                        if alt_key not in A3SS:
                            A3SS[alt_key] = (
                                gene,
                                gene_name,
                                chrom,
                                strand,
                                gene_struct[gene]["idx_exon"][i][0] - 1,
                                gene_struct[gene]["idx_exon"][i][1],
                                gene_struct[gene]["idx_exon"][j][0] - 1,
                                gene_struct[gene]["idx_exon"][j][1],
                                gene_struct[gene]["idx_exon"][idx][0] - 1,
                                gene_struct[gene]["idx_exon"][idx][1],
                            )

            for i in range(idx + 1, num_exon):
                if not check_edge(gene_struct, gene, idx, i):
                    continue
                for j in range(i + 1, num_exon):
                    if not check_edge(gene_struct, gene, idx, j):
                        continue
                    if gene_struct[gene]["idx_exon"][i][1] != gene_struct[gene]["idx_exon"][j][1]:
                        continue

                    alt_key = (
                        chrom,
                        gene_struct[gene]["idx_exon"][idx][1],
                        gene_struct[gene]["idx_exon"][i][0] - 1,
                        gene_struct[gene]["idx_exon"][j][0] - 1,
                    )

                    if gene_struct[gene]["strand"] == "-":
                        if alt_key not in A5SS:
                            A5SS[alt_key] = (
                                gene,
                                gene_name,
                                chrom,
                                strand,
                                gene_struct[gene]["idx_exon"][i][0] - 1,
                                gene_struct[gene]["idx_exon"][i][1],
                                gene_struct[gene]["idx_exon"][j][0] - 1,
                                gene_struct[gene]["idx_exon"][j][1],
                                gene_struct[gene]["idx_exon"][idx][0] - 1,
                                gene_struct[gene]["idx_exon"][idx][1],
                            )
                    elif gene_struct[gene]["strand"] == "+":
                        if alt_key not in A3SS:
                            A3SS[alt_key] = (
                                gene,
                                gene_name,
                                chrom,
                                strand,
                                gene_struct[gene]["idx_exon"][i][0] - 1,
                                gene_struct[gene]["idx_exon"][i][1],
                                gene_struct[gene]["idx_exon"][j][0] - 1,
                                gene_struct[gene]["idx_exon"][j][1],
                                gene_struct[gene]["idx_exon"][idx][0] - 1,
                                gene_struct[gene]["idx_exon"][idx][1],
                            )

    return A5SS, A3SS


def change_key(AS, AS_info):
    keys = AS.keys()

    for key in keys:
        x = AS_info(AS[key])
        AS[x.uniqID] = AS.pop(key)

    return AS


def detect_AS(gene_struct):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start: detecting AS event.")

    SE, MXE = detect_SE_MXE(gene_struct)
    RI = detect_RI(gene_struct)
    A5SS, A3SS = detect_AXSS(gene_struct)

    SE = change_key(SE, SE_info)
    MXE = change_key(MXE, MXE_info)
    RI = change_key(RI, RI_info)
    A5SS = change_key(A5SS, AXSS_info)
    A3SS = change_key(A3SS, AXSS_info)

    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "end  : detecting AS event.")
    return SE, MXE, RI, A5SS, A3SS


def check_exon_overlap(gene_struct, gene, exon, ignoreExon=""):
    for aexon in gene_struct[gene]["idx_exon"]:
        if aexon[1] < exon[0]:
            continue
        if aexon == exon:
            continue
        if aexon == ignoreExon:
            continue
        if aexon[0] > exon[1]:
            break
        else:
            return True
    return False


def update_idx_from_overlapping_exon(gene_struct, gene, uexon, dexon):
    idx_left = gene_struct[gene]["exon_idx"][uexon]
    idx_right = gene_struct[gene]["exon_idx"][dexon]
    new_left = idx_left
    new_right = idx_right

    for idx in range(idx_left + 1, idx_right):
        aexon = gene_struct[gene]["idx_exon"][idx]
        if (
            aexon[1] == gene_struct[gene]["idx_exon"][idx_left][1]
            and aexon[0] >= gene_struct[gene]["idx_exon"][idx_left][0]
        ):
            new_left = idx
        elif (
            aexon[0] == gene_struct[gene]["idx_exon"][idx_right][0]
            and aexon[1] <= gene_struct[gene]["idx_exon"][idx_right][1]
        ):
            new_right = idx
    return new_left, new_right


def get_basic_SE(SE, gene_struct):
    basic_SE_keys = set()

    for key in SE:
        x = SE_info(SE[key])
        texon = (x.exonStart_0base + 1, x.exonEnd)
        uexon = (x.upstreamES + 1, x.upstreamEE)
        dexon = (x.downstreamES + 1, x.downstreamEE)

        if check_exon_overlap(gene_struct, x.GeneID, texon):
            continue

        idx = gene_struct[x.GeneID]["exon_idx"][texon]
        uidx, didx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, uexon, dexon)

        if idx - uidx == 1 and didx - idx == 1:
            basic_SE_keys.add(key)

    return basic_SE_keys


def get_basic_MXE(MXE, gene_struct):
    basic_MXE_keys = set()

    for key in MXE:
        x = MXE_info(MXE[key])

        texon = (x.stExonStart_0base + 1, x.stExonEnd)
        aexon = (x.ndExonStart_0base + 1, x.ndExonEnd)
        uexon = (x.upstreamES + 1, x.upstreamEE)
        dexon = (x.downstreamES + 1, x.downstreamEE)

        if check_exon_overlap(gene_struct, x.GeneID, texon):
            continue
        if check_exon_overlap(gene_struct, x.GeneID, aexon):
            continue

        tidx = gene_struct[x.GeneID]["exon_idx"][texon]
        aidx = gene_struct[x.GeneID]["exon_idx"][aexon]
        uidx, didx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, uexon, dexon)

        if tidx - uidx == 1 and aidx - tidx == 1 and didx - aidx == 1:
            basic_MXE_keys.add(key)

    return basic_MXE_keys


def get_basic_RI(RI, gene_struct):
    basic_RI_keys = set()

    for key in RI:
        x = RI_info(RI[key])

        texon = (x.riExonStart_0base + 1, x.riExonEnd)
        uexon = (x.upstreamES + 1, x.upstreamEE)
        dexon = (x.downstreamES + 1, x.downstreamEE)
        uidx, didx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, uexon, dexon)
        tidx, didx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, texon, dexon)

        if tidx - uidx == 1 and didx - tidx == 1:
            basic_RI_keys.add(key)

    return basic_RI_keys


def get_basic_AXSS(AXSS, gene_struct):
    basic_AXSS_keys = set()

    for key in AXSS:
        x = AXSS_info(AXSS[key])

        lexon = (x.longExonStart_0base + 1, x.longExonEnd)
        sexon = (x.shortES + 1, x.shortEE)
        fexon = (x.flankingES + 1, x.flankingEE)

        if check_exon_overlap(gene_struct, x.GeneID, lexon, ignoreExon=sexon):
            continue

        if x.flankingEE <= x.longExonStart_0base:
            uidx, lidx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, fexon, lexon)
            uidx, sidx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, fexon, sexon)
            if lidx - uidx == 1 and sidx - lidx == 1:
                basic_AXSS_keys.add(key)
        elif x.flankingES >= x.longExonEnd:
            lidx, didx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, lexon, fexon)
            sidx, didx = update_idx_from_overlapping_exon(gene_struct, x.GeneID, sexon, fexon)
            if lidx - sidx == 1 and didx - lidx == 1:
                basic_AXSS_keys.add(key)

    return basic_AXSS_keys


def get_basic_AS(gene_struct, SE, MXE, RI, A5SS, A3SS):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start: filtering keys for basic AS event.")

    basic_SE_keys = get_basic_SE(SE, gene_struct)
    basic_MXE_keys = get_basic_MXE(MXE, gene_struct)
    basic_RI_keys = get_basic_RI(RI, gene_struct)
    basic_A5SS_keys = get_basic_AXSS(A5SS, gene_struct)
    basic_A3SS_keys = get_basic_AXSS(A3SS, gene_struct)

    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "end  : filtering keys for basic AS event.")
    return basic_SE_keys, basic_MXE_keys, basic_RI_keys, basic_A5SS_keys, basic_A3SS_keys


# def main():
#     # fn_gtf = "/mnt/isilon/xing_lab/wangy14/data/reference_genome/gencode.v31lift37.annotation.gtf"
#     fn_gtf = sys.argv[1]
#     gene_struct, tx2gene = parse_gtf(fn_gtf)
#     gene_struct = build_graph(gene_struct)
#     SE, MXE, RI, A5SS, A3SS = detect_AS(gene_struct)
#     basic_SE_keys, basic_MXE_keys, basic_RI_keys, basic_A5SS_keys, basic_A3SS_keys = get_basic_AS(
#         gene_struct, SE, MXE, RI, A5SS, A3SS
#     )


# if __name__ == "__main__":
#     main()
