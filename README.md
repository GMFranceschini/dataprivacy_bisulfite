

## anonymizeBAM.py: de-identification of sequencing reads

**anonymizeBAM.py** is a tool that can sanitize sequencing reads stored in BAM file format to protect the privacy and genetic information of donor individuals.

## Installation

anonymizeBAM.py is available through PyPI. To install, type the following command line, and add `-U` for upgrading:
`pip install anonymizeBAM`

Alternatively, you can install from this GitHub repository for the latest version:
`pip install git+https://github.com/sandberg-lab/dataprivacy`

Add `--user` if you don't have write permissions in your default python folder.

## Usage

anonymizeBAM.py requires only an aligned .bam file and the reference genome in fasta format.
Your fasta file should be indexed (`samtools faidx`).
The .bam file should be coordinate sorted and indexed, however `anonymizeBAM.py` will try to do this for you if not.


    usage: anonymizeBAM.py [-h] [--bam FILENAME] [--out FILENAME] [--fa FILENAME]
                            [--p P] [--strict] [--keepsecondary] [--keepunmapped]

    optional arguments:
      -h, --help      show this help message and exit
      --bam FILENAME  Path to input BAM file
      --out FILENAME  Path to output bam file
      --fa FILENAME   Path to genome reference fasta
      --p P           Number of processes to use
      --strict        Strict: also sanitize mapping score & auxiliary tags (eg. AS / NH).
      --keepsecondary  Keep secondary alignments in output bam file.
      --keepunmapped   Keep ummapped reads in output bam file.

## Description

**anonymizeBAM** sanitizes sequence reads to provide privacy protection and facilitate data sharing.

The anonymization procedure involves modification of the observed read sequence to the reference genome sequence and sanitation of auxiliary tags.

Here is an overview of the sequence correction strategy:

 1. SNPs: Mismatches to the reference (either explicitly *X* coded in the CIGAR value or within *M* matched segments) are replaced by the reference base.
 2. Insertions: The read sequence is extended by the length equal to the insertion while keeping the 5' mapping position constant.
 3. Deletions: The missing reference sequence is inserted into the read while removing an equal numbers of bases from the 3’ end.
 4. Clipping: soft or hard clipped bases (CIGAR: *S* / *H*) are replaced by matching reference sequence. If reads start with clipped bases in single-end data, the reference position of the read start is adjusted, however this is not possible for paired-end reads because it would invalidate the mate-pair information (TLEN and PNEXT fields). Instead for paired-end reads, the clipped sequence portion is added to the end of the read.
  5. Splicing: Splicing is observed and splice-sites are conserved even in the case of deletions and insertions.
  6. Multimapping: In the default behavior, only primary alignments are emitted. The user can choose to keep secondary but note that anonymization can not be guaranteed.
  7. Unmapped reads: Unmapped reads cannot be sanitized and are discarded in default settings.

Donor-related information could also be inferred from standard bam fields and auxiliary tags:

 1. CIGAR value is matched to the anonymized sequence (eg. 100M).
 2. MD are matched to the anonymized sequence, if present (eg. 100) .
 3. NM and nM tags are sanitized by replacement with 1.
 4. Tags containing information on the alignment are discarded (MC, XN, XM, XO, XG)

In `--strict` mode, the following tags are also changed:

 5. Mapping quality set to max/unavailable (255)
 6. AS and MQ are set to read length
 7. NH is set to 1
 8. Discarding of the following tags: HI, IH, H1, H2, OA, OC, OP, OQ, SA, SM, XA, XS



## Reference

[link to follow]
