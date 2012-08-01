#!/usr/bin/python

import utils
import os
import sys
import random
import subprocess

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "USAGE: cross_validate.py <n-fold> <class1_samples> <class2_samples>..."
        sys.exit(1)

    seqs, labels, all_labels = utils.read_all_labeled(sys.argv[2:])
    folds = int(sys.argv[1])
    
    random.shuffle(seqs)
    errors = utils.cross_validate(seqs, labels, folds)

    print "Fold errors:", " ".join([str(x) for x in errors])
    print "Average error: %3.2f%%" % (sum(errors) * 100.0 / len(errors))


