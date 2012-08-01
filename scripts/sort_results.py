#!/usr/bin/python

from Bio import Seq
import os
import sys

def read_fasta(f):
    seqs = dict()

    seqid = None
    seq = ""
    header = None

    for line in open(f):
        line = line.strip()
        if line[0] == ">":
            if seqid != None:
                seqs[seqid] = (header, seq)
            i = line.find(" ")
            if i != -1:
                seqid = line[1:i]
                header = line[i:]
            else:
                seqid = line[1:]
                
            seq = ""
        else:
            seq += line.replace(" ", "")

    if seqid != None:
        seqs[seqid] = (header, seq)

    return seqs

def write_seq(s, seqid, dir, seqs):
    header, seq = seqs[seqid]

    if dir == "-":
        seq = Seq.Seq(seq).reverse_complement()

    if header:
        l = ">%s %s" % (seqid, header)
    else:
        l = ">%s" % seqid

    s.write("%s\n%s\n" % (l, seq))

def main(in_fasta, in_score):
    stem = in_fasta[:in_fasta.rfind(".")]
    stem = os.path.split(stem)[1]

    out_streams = dict()

    seqs = read_fasta(in_fasta)

    for line in open(in_score):
        if line[0] == "#":
            continue
        lexemes = line.strip().split("\t")
        if len(lexemes) < 5:
            continue

        seqid = lexemes[0]
        dir = lexemes[2]
        c = lexemes[-1]

        if c not in out_streams:
            out_streams[c] = open(stem + "_" + c + ".fasta", "w")

        write_seq(out_streams[c], seqid, dir, seqs)

    for c in out_streams:
        out_streams[c].close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print >>sys.stderr, "USAGE: sort_results.py <query_seqs> <result_file>"
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
