#!/usr/bin/python

import utils
import sys

seqs_out = open("test_seqs.fasta", "w")
labels_out = open("labels.txt", "w")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print >>sys.stderr, "USAGE: make_test.py <fasta1>,lab1 [fasta2,lab2 ...]"
        sys.exit(1)

    for test in sys.argv[1:]:
        lexemes = test.split(",")
        if len(lexemes) != 2:
            print >>sys.stderr, "Expected fasta file and label seperated by ,"
            sys.exit(1)

        fasta = lexemes[0]
        lab = lexemes[1]

        for seq in utils.read_fasta(fasta):
            labels_out.write("%s\t%s\n" % (seq["seqid"], lab))
            seqs_out.write(">%s\n%s\n" % (seq["seqid"], seq["seq"]))
