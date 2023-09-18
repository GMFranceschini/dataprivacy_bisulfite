# BAMboozle.py - de-identify sequencing data
# Author: Christoph Ziegenhain / christoph.ziegenhain@ki.se
# Last update: 11-01-2021

import argparse
import itertools
import multiprocessing as mp
import os

import pysam


def makeBAMheader(args, v):
    bam = pysam.AlignmentFile(args.bam, "rb")
    hdr = bam.header.to_dict()
    bam.close()

    cmdlinecall = (
        "BAMboozle --bam "
        + args.bam
        + " --out "
        + args.out
        + " --fa "
        + args.fa
        + " --p "
        + str(args.p)
    )

    if args.strict:
        cmdlinecall = cmdlinecall + " --strict"
    if args.keepunmapped:
        cmdlinecall = cmdlinecall + " --keepunmapped"
    if args.keepsecondary:
        cmdlinecall = cmdlinecall + " --keepsecondary"

    pg = {"ID": "BAMboozle", "PN": "BAMboozle", "CL": cmdlinecall, "VN": v}

    if "PG" in hdr:
        pglines = hdr["PG"]
        pglines.append(pg)
    else:
        pglines = [pg]
    hdr["PG"] = pglines

    return hdr


def idx_bam(bam, threads):
    threads = str(threads)
    try:
        pysam.index(f"-@{threads}", bam)
    except:
        outcome = "idxerror"
    else:
        outcome = "idxsuccess"
    if outcome == "idxerror":
        print("indexing failed, trying to sort bam file...")
        inbam = bam
        bam = f"{bam}.sorted.bam"
        pysam.sort(f"-@{threads}", "-o", bam, inbam)
        print("indexing bam file...")
        pysam.index(bam)
    return bam


def collect_bam_chunks(inpath, chrs, outpath, unmapped):
    allpaths = [f"{inpath}.tmp.{c}.bam" for c in chrs]
    if unmapped:
        allpaths = [f"{inpath}.tmp.{c}.bam" for c in chrs[:-1]]
        allpaths.append(f"{inpath}.tmp.unmapped.bam")
    cat_args = ["-o", outpath] + allpaths
    pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]


def remove_tag(read, rtag):
    all_tags = read.get_tags()
    to_keep = [t[0] != rtag for t in all_tags]
    kept_tags = [tag for tag, keep in zip(all_tags, to_keep) if keep]
    read.set_tags(kept_tags)
    return read


def count_ref_consuming_bases(cigartuples):
    bases = 0
    for cig in cigartuples:
        if cig[0] in [0, 2, 7, 8]:
            bases = bases + cig[1]
    return bases

def tag_fix(incigar, bis_tag):
    new_tag = ""
    curr_tag = 0
    
    for cig, length in incigar:
        if cig == 0:
            new_tag += bis_tag[curr_tag:curr_tag+length]
            curr_tag = curr_tag + length
        elif cig == 1: # Insertion
            curr_tag = curr_tag + length
        elif cig == 2: # Deletion
            new_tag = new_tag + "Z"*length
    new_tag = "".join(new_tag)
    return(new_tag)



def clean_bam(
    inpath, threads, fastapath, chr, strict, keepunmapped, keepsecondary, anonheader
):
    fa = pysam.FastaFile(fastapath)

    if chr == "*":
        chrlabel = "unmapped"
    else:
        chrlabel = chr

    # open in/out files
    outpath = f"{inpath}.tmp.{chrlabel}.bam"
    tmppath = f"{inpath}.tmp.{chrlabel}.tmp.bam"
    inp = pysam.AlignmentFile(inpath, "rb", threads=threads)
    out = pysam.AlignmentFile(tmppath, "wb", header=anonheader, threads=threads)
    for read in inp.fetch(chr):
        # deal with unmapped reads
        if chrlabel == "unmapped":
            trim_tags = ["uT", "nM", "NM", "XN", "XM", "XO", "XG"]
            if strict:
                trim_tags += [
                    "NH",
                    "HI",
                    "IH",
                    "AS",
                    "MQ",
                    "H1",
                    "H2",
                    "OA",
                    "OC",
                    "OP",
                    "OQ",
                    "SA",
                    "SM",
                    "XA",
                    "XS",
                ]
            for t in trim_tags:
                if read.has_tag(t):
                    read = remove_tag(read, t)
            out.write(read)
            continue

        # only use primary alignments
        if not keepsecondary and read.is_secondary:
            continue

        # determine some basics
        readlen = read.query_length
        qual = read.query_qualities

        if read.is_paired:
            readtype = "PE"
            readtype_int = 2
        else:
            readtype = "SE"
            readtype_int = 1

        # modify tags
        trim_tags = ["MC", "XN", "XO", "XG"]

        


        for t in trim_tags:
            if read.has_tag(t):
                read = remove_tag(read, t)

        # some aligners like to keep the unmapped reads at the same posiiton as their mate, deal with them
        if read.is_unmapped:
            if keepunmapped:
                out.write(read)
            continue

        # look at cigar value
        incigar = read.cigartuples
        present_cigar_types = [x[0] for x in incigar]

        if present_cigar_types == [0]:
            fa_ref = fa.fetch(chr, read.reference_start, read.reference_start+readlen)
            seq = read.query_sequence
            btag = read.get_tags()[2][1]
            read.query_sequence = "".join([s if b.upper() == "Z" else r for s, r, b in zip(seq, fa_ref, btag)])
            read.set_tag(tag="MD", value_type="Z", value=read.query_sequence)
            final_outseq = read.query_sequence
            final_cigar = [(0, readlen)]
        else:
            fa_ref = fa.fetch(chr, read.reference_start, read.reference_start+readlen)
            seq = read.query_sequence
            btag = read.get_tags()[2][1]
            btag = tag_fix(incigar=incigar, bis_tag=btag)
            read.query_sequence = "".join([s if b.upper() == "Z" else r for s, r, b in zip(seq, fa_ref, btag)])
            read.set_tag(tag="MD", value_type="Z", value=read.query_sequence)
            final_outseq = []
            final_qual = []
            curr_query = 0
            curr_ref = 0

            for op, length in incigar:
                if op == 0:
                    final_outseq.extend(read.query_sequence[curr_ref:curr_ref+length])
                    final_qual.extend(qual[curr_query:curr_query+length])
                    curr_ref = curr_ref + length
                    curr_query = curr_query + length
                elif op == 1: # Insertion
                    curr_query = curr_query + length
                elif op == 2: # Deletion
                    final_outseq.extend("N" * length)
                    final_qual.extend([38] * length)
                    curr_ref = curr_ref + length
                    
            final_outseq = "".join(final_outseq)
            final_outseq = final_outseq.upper()
            
            qual = final_qual
            final_cigar = [(0, len(qual))]


        if len(final_outseq) != len(
            qual
        ):  # the sanitized output sequence cannot be longer than a contig (reason:deletions)
            len_diff = len(qual) - len(final_outseq)
            old_len = final_cigar[-1][1]
            new_len = old_len - len_diff
            final_cigar[-1] = (0, new_len)
            qual = qual[: len(final_outseq)]

        # set alignment record and write
        read.query_sequence = final_outseq
        read.query_qualities = qual
        read.cigartuples = final_cigar
        out.write(read)

    inp.close()
    out.close()

    # resorting SE reads is necessary since we may mess with read start positions
    try:
        test = readtype
    except NameError:
        readtype = "NA"
    if readtype == "SE":
        pysam.sort("-o", outpath, tmppath)
        if os.path.exists(tmppath):
            os.remove(tmppath)
    else:
        os.rename(tmppath, outpath)


def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        "--bam",
        type=str,
        metavar="FILENAME",
        help="Path to input BAM file",
        required=True,
    )
    parser.add_argument(
        "--out",
        type=str,
        metavar="FILENAME",
        help="Path to output bam file",
        required=True,
    )
    parser.add_argument(
        "--fa",
        type=str,
        metavar="FILENAME",
        help="Path to genome reference fasta",
        required=True,
    )
    parser.add_argument("--p", type=int, default=10, help="Number of processes to use")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Strict: also sanitize mapping score & auxiliary tags (eg. AS / NH).",
    )
    parser.add_argument(
        "--keepsecondary",
        action="store_true",
        help="Keep secondary alignments in output bam file.",
    )
    parser.add_argument(
        "--keepunmapped",
        action="store_true",
        help="Keep ummapped reads in output bam file.",
    )

    args = parser.parse_args()
    v = "0.5.0"
    print("BAMboozle.py v" + v)
    bampath = args.bam
    try:
        fa = pysam.FastaFile(args.fa)
    except ValueError:
        print("Error: Reference fasta file is not indexed!")
        print("Please run: samtools faidx " + args.fa)
        quit()

    if not os.path.exists(bampath + ".bai"):
        print("input bam index not found, indexing...")
        bampath = idx_bam(bampath, args.p)

    # Construct the new bam header to work with
    bamheader = makeBAMheader(args, v)
    print("Working...")

    chrs = pysam.idxstats(bampath).split("\n")
    chrs = [c.split("\t")[0] for c in chrs[:-1]]
    chrs = [
        c for c in chrs if c in fa.references
    ]  # we can only sanitize contigs that have a reference seq
    if args.keepunmapped:
        chrs.append("*")
    fa.close()

    if args.p > 20:
        pysam_workers = 2
        n_jobs = int(args.p / 2)
    else:
        pysam_workers = 1
        n_jobs = args.p

    pool = mp.Pool(n_jobs)
    results = [
        pool.apply_async(
            clean_bam,
            (
                args.bam,
                pysam_workers,
                args.fa,
                chr,
                args.strict,
                args.keepunmapped,
                args.keepsecondary,
                bamheader,
            ),
        )
        for chr in chrs
    ]
    x = [r.get() for r in results]
    # single threaded below:
    # [clean_bam(bampath,pysam_workers,args.fa,chr,args.strict) for chr in chrs]

    print("Creating final output .bam file...")
    collect_bam_chunks(
        inpath=bampath, chrs=chrs, outpath=args.out, unmapped=args.keepunmapped
    )
    print("Indexing final output .bam file...")
    y = idx_bam(args.out, args.p)

    print("Done!")


if __name__ == "__main__":
    main()
