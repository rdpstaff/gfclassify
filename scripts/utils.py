#!/usr/bin/python

import os
import sys
import random
import subprocess
import math

def read_all_labeled(seq_files):
    seqs = []
    labels = dict()
    all_labels = set()

    for f in seq_files:
        label = ".".join(os.path.split(f)[-1].split(".")[:-1])
        all_labels.add(label)

        for seq in read_fasta(f):
            seqs.append(seq)
            labels[seq["seqid"]] = label

    return seqs, labels, all_labels

def read_fasta(seq_file):
    ret = []

    seq = None
    seqid = None
    for line in open(seq_file):
        if line[0] == ">":
            if seq != None and seqid != None:
                ret.append({"seqid" : seqid, "seq" : seq})
            seqid = line.split()[0][1:]
            seq = ""
        else:
            if seq == None:
                raise IOError("First life of file %s wasn't a fasta header" % (seq_file))
            seq += line.strip()

    if seq != None and seqid != None:
        ret.append({"seqid" : seqid, "seq" : seq.lower()})

    return ret

def cross_validate(seqs, labels, k, period = 1, window = 12, depth = 7):
    fold_errors = []
    for i in range(k):
        training, testing = get_cv_split(k, i, seqs)

#	hmms = build_hmms(training, labels)
#	errors = test_hmm(testing, labels, hmms)
#        print errors

        models = build_models(training, labels, period, window, depth)
        errors, log_odds = test_model(testing, labels, models)

        fold_errors.append(errors)

    return fold_errors

def leave_one_out(seqs, labels, period = 1, window = 12, depth = 7):
    errors = []
    for i in range(len(seqs)):
        testing = [seqs[i]]
        training = seqs[:i] + seqs[i+1:]

        models = build_models(training, labels, period, window, depth)
        e, log_odds = test_model(testing, labels, models)

        if e == 1:
            errors.append(seqs[i]["seqid"])

    return errors
    
def write_seqs(seqs, out):
    o = open(out, "w")
    for seq in seqs:
        o.write(">%s\n%s\n" % (seq["seqid"], seq["seq"]))
    o.close()

def get_cv_split(k, fold, seqs):
    if fold > k:
        raise ValueError("Fold cannot be higher than total fold cv")

    if k == 1:
        return seqs, seqs

    s = len(seqs) / float(k)
    idx = s * fold
    l = int(round(idx))
    u = int(round(idx + s))

    testing = seqs[l:u]
    training = seqs[:l] + seqs[u:]

    return training, testing

def build_models(training_seqs, labels, period = 1, window = 12, depth = 7):
    training_streams = dict()
    file_pattern = "/tmp/%s_training.fasta"
    icm_pattern = "/tmp/%s.icm"

    for seq in training_seqs:
        label = labels[seq["seqid"]]
        if label not in training_streams:
            training_streams[label] = open(file_pattern % label, "w")

        training_streams[label].write(">%s\n%s\n" % (seq["seqid"], seq["seq"]))

    ret = []
    for label in training_streams.keys():
        training_streams[label].close()
        f = icm_pattern % label

        if os.path.exists(f):
            os.remove(f)

        subprocess.check_call(["build-icm", "-p", str(period), "-w", str(window), "-d", str(depth), f], stdin=open(file_pattern % label))
        ret.append(f)

    return ret

def build_hmms(training_seqs, labels):
    training_streams = dict()
    file_pattern = "/tmp/%s_training.fasta"
    icm_pattern = "/tmp/%s.hmm"

    for seq in training_seqs:
        label = labels[seq["seqid"]]
        if label not in training_streams:
            training_streams[label] = open(file_pattern % label, "w")

        training_streams[label].write(">%s\n%s\n" % (seq["seqid"], seq["seq"]))

    ret = []
    for label in training_streams.keys():
        training_streams[label].close()
        f = icm_pattern % label

        if os.path.exists(f):
            os.remove(f)

        subprocess.check_call(["muscle", "-clw", "-in", file_pattern % label, "-out", "/tmp/aln.clw"], stdout=open("/dev/null"), stderr=open("/dev/null"))
        s = open("/tmp/aln.clw")
        clw = s.read().replace("*", "")
        s.close()

        out = open("/tmp/aln.sto", "w")
        out.write("# STOCKHOLM 1.0\n#")
        out.write(clw)
        out.write("\n")
        out.write("//\n")
        out.flush()
        out.close()        

        subprocess.check_call(["hmmbuild", "-n", label, f, "/tmp/aln.sto"], stdout=open("/dev/null"))
        ret.append(f)

    return ret

def test_hmm(testing, labels, models):
    testing_file = "/tmp/testing.fasta"
    write_seqs(testing, testing_file)

    best_labels = dict()

    for model in models:
        cmd = ["hmmsearch", "--tblout", "/tmp/output.txt", "-o", "/dev/null", model, testing_file]
        subprocess.check_call(cmd)


        for line in open("/tmp/output.txt"):
            if line[0] == "#":
                continue
            lexemes = line.strip().split()

            if len(lexemes) == 0:
                continue

            seqid = lexemes[0]
            pred = lexemes[2]
            score = float(lexemes[5])

            if seqid not in best_labels:
                best_labels[seqid] = [0, ""]

            if score > best_labels[seqid][0]:
                best_labels[seqid] = [score, pred]

    errors = 0
    for seqid in best_labels:
        pred = best_labels[seqid][1]

        if pred != labels[seqid]:
            errors += 1.0

    print "HMMER3 ERROR: %s, %s, %s" % (errors, (len(testing) - len(best_labels)), (errors + (len(testing) - len(best_labels))) / len(testing))
    errors += (len(testing) - len(best_labels))

    os.remove(testing_file)
    return errors / len(testing)

def test_model(testing, labels, models):
    testing_file = "/tmp/testing.fasta"
    write_seqs(testing, testing_file)

    out = open("/tmp/output.txt", "w")
    cmd = ["gf_classify"]
    cmd.extend(models)
    p = subprocess.check_call(cmd, stdin=open(testing_file), stdout=out)
    out.close()

    log_odds = dict()
    errors = 0

    for line in open("/tmp/output.txt"):
        if line[0] == "#":
            continue
        lexemes = line.strip().split("\t")

        if len(lexemes) == 0:
            continue

        seqid = lexemes[0]
        desc = lexemes[1]
        direction = lexemes[2]
        scores = [float(x) for x in lexemes[3:-1]]
        pred = lexemes[-1]

        log_odds[seqid] = scores

        if pred != labels[seqid]:
            errors += 1.0

    os.remove(testing_file)
    return errors / len(testing), log_odds
