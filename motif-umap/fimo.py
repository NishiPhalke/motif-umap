#!/usr/bin/env python

import os
import tempfile

class FIMO:

    def __init__(self, fa, motif, maxScores = 1000000):
        self.fa = fa
        self.motif = motif
        self.maxScores = maxScores
    
    def __enter__(self):
        self.tempdir = tempfile.TemporaryDirectory()
        os.system("fimo --parse-genomic-coord --max-stored-scores {maxScores} --oc {output} {motif} {fa}".format(
            maxScores = self.maxScores,
            output = self.tempdir.name,
            motif = self.motif,
            fa = self.fa
        ))
        self.fimo = open(os.path.join(self.tempdir.name, "fimo.tsv"), 'rt')
        return self.fimo

    def __exit__(self, *args):
        self.fimo.close()
        self.tempdir.cleanup()
