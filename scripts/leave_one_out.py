#!/usr/bin/python

import utils
import os
import sys
import random
import subprocess

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "USAGE: leave_one_out.py <class1_samples> <class2_samples>..."
        sys.exit(1)

    seqs, labels, all_labels = utils.read_all_labeled(sys.argv[1:])

    errors = utils.leave_one_out(seqs, labels)

    print "Misclassified:", " ".join(errors)
    print "Average error: %2.2f%%" % (len(errors) / len(seqs))
