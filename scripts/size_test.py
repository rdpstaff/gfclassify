#!/usr/bin/python

import utils
import os
import sys
import random
import subprocess

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "USAGE: size_test.py <class1_samples> <class2_samples> ..."
        sys.exit(1)

    seqs = dict()
    labels = dict()
    all_labels = set()

    for f in sys.argv[1:]:
        label = ".".join(os.path.split(f)[-1].split(".")[:-1])
        all_labels.add(label)
        seqs[label] = []

        for seq in utils.read_fasta(f):
            seqs[label].append(seq)
            labels[seq["seqid"]] = label

    for label in seqs.keys():
        random.shuffle(seqs[label])

    errors = []

    sizes = []
    ratios = [.01, .05, .1, .25, .5, .75, .9];
    for ratio in ratios:
        training = []
        testing = []
        for label in seqs.keys():
            #random.shuffle(seqs[label]) #So our training sets are not proper subsets of eachother
            k = int(ratio * len(seqs[label]))
            training.extend(seqs[label][:k])
            testing.extend(seqs[label][k:])

        print >>sys.stderr, "Ratio %s, training: %s, testing: %s" % (ratio, len(training), len(testing))
        sizes.append(len(training))

        models = utils.build_models(training, labels)
	error, log_odds = utils.test_model(testing, labels, models)
        errors.append(error)

    #print errors

    print sizes
    print errors

