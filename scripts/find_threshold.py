#!/usr/bin/python

import numpy
import utils
import sys

if len(sys.argv) != 3:
    print >>sys.stderr, "USAGE: find_threshold.py <labels> <results>"
    sys.exit(1)

labels = dict()
for line in open(sys.argv[1]):
    lexemes = line.strip().split()
    if len(lexemes) == 2:
        labels[lexemes[0]] = lexemes[1]

c = []
m = []

for line in open(sys.argv[2]):
    if line[0] == "#":
        continue

    lexemes = line.strip().split("\t")

    if len(lexemes) < 5:
        continue

    query = lexemes[0]
    scores = [float(x) for x in lexemes[3:-1]]
    pred = lexemes[-1]
    real = labels[query]

    best = max(scores)
    second_best = -10000000
    for score in scores:
        if score == best:
            continue
        if score > second_best:
            second_best = score

    ll = best - second_best

    if pred == "rejected":
        raise ValueError("Cannot predict threshold from already thresholded predictions")
    elif pred == real:
        c.append(ll)
    else:
        m.append(ll)

print min(c), numpy.mean(c), max(c), numpy.std(c)
print min(m), numpy.mean(m), max(m), numpy.std(m)

