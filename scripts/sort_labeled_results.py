#!/usr/bin/python

import utils
import sys

if len(sys.argv) != 4:
    print >>sys.stderr, "USAGE: sort_results.py <seqs> <labels> <results>"
    sys.exit(1)

labels = dict()
for line in open(sys.argv[2]):
    lexemes = line.strip().split()
    if len(lexemes) == 2:
        labels[lexemes[0]] = lexemes[1]

rej = open("rejected.txt", "w");
misclassified = open("misclassified.txt", "w");
cor = open("correct.txt", "w");
rej_seq = open("rejected.fna", "w");
misclassified_seq = open("misclassified.fna", "w");
cor_seq = open("correct.fna", "w");

seqs = dict()

for seq in utils.read_fasta(sys.argv[1]):
    seqs[seq["seqid"]] = seq

r = 0
c = 0
m = 0
for line in open(sys.argv[3]):
    if line[0] == "#":
        continue

    lexemes = line.strip().split("\t")

    if len(lexemes) < 6:
        continue

    query = lexemes[0]
    pred = lexemes[-1]
    real = labels[query]

    if pred == "rejected":
        rej.write(line)
        rej_seq.write(">%s\n%s\n" % (seqs[query]["seqid"], seqs[query]["seq"]))
        r += 1
    elif pred == real:
        cor.write(line)
        cor_seq.write(">%s\n%s\n" % (seqs[query]["seqid"], seqs[query]["seq"]))
        c += 1
    else:
        misclassified.write(line)
        misclassified_seq.write(">%s\n%s\n" % (seqs[query]["seqid"], seqs[query]["seq"]))
        m += 1

tot = float(r + c + m)
print "Correctly classified %d (%f), misclassified %d (%f), rejected %d (%f)" % (c, c/tot, m, m/tot, r, r/tot)
