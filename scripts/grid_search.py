#!/usr/bin/python

import utils
import os
import sys
import random
import argparse

def get_range(a):
    lexemes = a.split(",")

    if len(lexemes) == 1:
        return [int(a)]
    elif len(lexemes) == 2:
        return range(int(lexemes[0]), int(lexemes[1]) + 1)
    elif len(lexemes) == 3:
        return range(int(lexemes[0]), int(lexemes[1]) + 1, int(lexemes[2]))        

    raise Exception("Cannot parse %s in to a range" % a)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform a grid search to find good parameters for build-icm.")
    parser.add_argument('--window', default=12, help="min,max")
    parser.add_argument('--depth', default=7, help="min,max")
    parser.add_argument('--period', default=1, help="min, max")
    parser.add_argument('-k', default=1, help="k-fold validation")
    parser.add_argument('model', nargs='+')

    args = parser.parse_args()

    seqs = []
    labels = dict()
    all_labels = set()

    for f in args.model:
        label = ".".join(os.path.split(f)[-1].split(".")[:-1])
        all_labels.add(label)

        for seq in utils.read_fasta(f):
            seqs.append(seq)
            labels[seq["seqid"]] = label

    best_args = [12, 7, 1, 100]

    print "width\tdepth\tperiod\terror_rate"
    for w in get_range(args.window):
        for d in get_range(args.depth):
            if w < d:
                continue

            for p in get_range(args.period):
                errors = utils.cross_validate(seqs, labels, int(args.k), window=w, depth=d, period=p)
                avg_error = sum(errors) / len(errors)

                print "%s\t%s\t%s\t%0.3f" % (w, d, p, avg_error)
                sys.stdout.flush()
                if avg_error < best_args[-1]:
                    best_args = [w, d, p, avg_error]

    print "Best paramters: w=%s, d=%s, p=%s with error of %0.3f" % (best_args[0], best_args[1], best_args[2], best_args[3])


