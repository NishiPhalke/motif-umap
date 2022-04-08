#!/usr/bin/env python3

import os
import sys
import argparse
import glob
import tempfile

import umap
import numpy
import random

from fa import FASTAFile
from fimo import FIMO
from metadata import reference_epigenome_signal_files, reference_epigenome
from matrix import signalmatrix

SCALE_FACTORS = { "H3K9me3": 10, "H3K27me3": 10, "H3K27ac": 1, "H3K4me3": 1, "H3K4me1": 3, "WGBS": 2, "CTCF": 4 }

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-epigenomes", nargs = '+', help = "ENCODE accessions for reference epigenomes to use", required = True)
    parser.add_argument("--motif-files", nargs = '+', help = "path to motif files in MEME format", required = True)
    parser.add_argument("--marks", nargs = '+', help = "epigenomic marks to use", default = [ "H3K4me3", "H3K27ac", "H3K4me1", "WGBS", "H3K27me3", "H3K9me3", "CTCF" ])
    parser.add_argument("--encode-mirror-format", help = "format string for paths to local ENCODE files", default = "/data/pdisk/signal/%s/%s.bigWig")
    parser.add_argument("--encode-Z-score-directory", help = "path to Z-score directory", default = "/data/pdisk/Signal-Files")
    parser.add_argument("--rDHS-coordinates", help = "path to rDHS coordinates", default = "/data/common/genome/GRCh38-rDHSs.bed")
    parser.add_argument("--output-directory", help = "path to output directory", required = True)
    parser.add_argument("--extra-signal-files", help = "path to extra signal files", nargs = '+', default = [])
    parser.add_argument("--extra-motif-sites", help = "optional path to extra motif sites", nargs = '+', default = [])
    parser.add_argument("--skip-active-filtering", help = "if set, cluser all motif sites, not just chromatin accessible ones", action = "store_true", default = False)
    return parser.parse_args()

def generate_if_necessary(v, p):
    if os.path.exists(p): return numpy.load(p)
    numpy.save(p, v())
    return numpy.load(p)

def map_local(reference_epigenome, mirror_format):
    def ensure_local(path, a):
        if not os.path.exists(path):
            os.system("wget https://www.encodeproject.org/files/%s/@@download/%s.bigWig -O %s" % (a, a, path))
        return path
    files = reference_epigenome_signal_files(reference_epigenome)
    return { k: ensure_local(mirror_format % v, v[1]) for k, v in files.items() }

def main():

    args = parse_args()
    args.marks = set(args.marks)

    print("ensuring local signal files are present for %d marks in %d epigenomes..." % (len(args.marks), len(args.reference_epigenomes)), file = sys.stderr)
    paths = { e: map_local(e, args.encode_mirror_format) for e in args.reference_epigenomes }

    print("getting active rDHSs for %d reference epigenomes..." % len(args.reference_epigenomes), file = sys.stderr)
    with open(args.rDHS_coordinates, 'r') as f:
        coordinates = { line.strip().split()[3]: tuple(line.strip().split()[:3]) for line in f }
    d = { e: reference_epigenome(e)["DNase-seq"] for e in args.reference_epigenomes }
    active = set()
    if not args.skip_active_filtering:
        for k, v in d.items():
            with open(glob.glob(os.path.join(args.encode_Z_score_directory, "%s-*" % v["accession"]))[0], 'rt') as f:
                active = active.union({ x.split()[0] for x in f if float(x.split()[1]) > 1.64 })
    else:
        with open(args.rDHS_coordinates, 'r') as f:
            active = { line.strip().split()[3] for line in f }
    print("found %d active rDHSs..." % len(active), file = sys.stderr)

    print("identifying motif instances...", file = sys.stderr)
    all_motifs = []
    with tempfile.NamedTemporaryFile('wt') as o:
        for x in active:
            o.write("%s\t%s\n" % ('\t'.join(coordinates[x]), x))
        o.flush()
        if not args.skip_active_filtering:
            with FASTAFile(o.name, "/data/common/genome/hg38.2bit") as fa:
                for motif in args.motif_files:
                    with FIMO(fa.name, motif) as fimo:
                        for line in fimo:
                            try:
                                all_motifs.append(( line.split()[1], int(line.split()[2]), int(line.split()[3]), line.split()[4] ))
                            except:
                                pass
        for x in args.extra_motif_sites:
            with tempfile.NamedTemporaryFile('rt') as f:
                os.system("bedtools intersect -a %s -b %s -wa | sort | uniq > %s" % (x, o.name, f.name))
                for line in f:
                    all_motifs.append(( line.split()[0], int(line.split()[1]), int(line.split()[2]),line.split()[3], line.split()[5] ))
    with open(os.path.join(args.output_directory, "motifs.bed"), 'wt') as o:
        o.write('\n'.join([ '\t'.join([ str(xx) for xx in x ]) for x in all_motifs ]) + '\n')
    random.shuffle(all_motifs)
    print("found %d total motif instances..." % len(all_motifs), file = sys.stderr)

    print("generating signal matrices...", file = sys.stderr)
    if not os.path.exists(os.path.join(args.output_directory, "motifs.150k.bed")):
        with open(os.path.join(args.output_directory, "motifs.150k.bed"), 'wt') as o:
            o.write('\n'.join([ '\t'.join([ str(xx) for xx in x ]) for x in all_motifs[:150000] ]) + '\n')
    minus_indexes = numpy.array([ i for i, x in enumerate(all_motifs[:150000]) if x[-1] == '-' ])
    def generate_signal(b):
        m = numpy.array( signalmatrix(os.path.join(args.output_directory, "motifs.bed"), b) )
        print(m.shape, file = sys.stderr)
        m[minus_indexes] = numpy.flip(m[minus_indexes], axis = 1)
        return m
    values = []
    for e, vv in paths.items():
        for k, v in vv.items():
            print(v + "...", file = sys.stderr)
            values.append(( e + "_" + k, generate_if_necessary(lambda: generate_signal(v), os.path.join(args.output_directory, os.path.basename(v).replace(".bigWig", ".npy"))) * SCALE_FACTORS[k] ))
    for sf in args.extra_signal_files:
        values.append(( sf, generate_if_necessary(lambda: generate_signal(sf), os.path.join(args.output_directory, os.path.basename(sf) + ".npy")) ))

    allvalues = generate_if_necessary(lambda: numpy.concatenate([ x[1] for x in values ], axis = 1), os.path.join(args.output_directory, "all-values.npy"))
    u = umap.UMAP(n_neighbors = 12, min_dist = 0.01)
    with open(os.path.join(args.output_directory, "mark-order.txt"), 'w') as o:
        o.write('\n'.join([ x[0] for x in values ]) + '\n' + '\n'.join(args.extra_signal_files) + '\n')
    generate_if_necessary(lambda: u.fit_transform(allvalues), os.path.join(args.output_directory, "umap.npy"))
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
