#!/usr/bin/env python

import numpy
import pyBigWig

from joblib import Parallel, delayed

def lsum(lists):
  r = [ [ xx for xx in x ] for x in lists[0] ]
  for x in lists[1:]:
    for i, xx in enumerate(x):
      r[i] += xx
  return numpy.array(r)

def flatten(l):
  r = []
  for x in l:
    r += x
  return r

def recenter(region, expsize):
  m = int((region[0] + region[1]) / 2)
  return ( m - expsize, m + expsize )

def chunk(lst, n):
  for i in range(0, len(lst), n):
    yield lst[i:i + n]

def bwreadg(bigwig, rr, s = 1):
  rx = []
  with pyBigWig.open(bigwig) as bw:
    for i, ra in enumerate(rr):
      r = ra[1]
      c = ra[0]
      if r[0] < 1:
        rx.append([ 0 for _ in range(r[0], r[1]) ])
        continue
      try:
        x = [ x * s if not numpy.isnan(x) else 0 for x in bw.values(c, r[0], r[1]) ]
        if ra[-1] == '-': x = reversed(x)
        rx.append(x)
      except:
        rx.append([ 0 for _ in range(r[0], r[1]) ])
  return rx

def averagem(s, r):
  rr = []
  for i in range(int(len(s) / r)):
    x = s[i * r : (i + 1) * r]
    rr.append(float(sum(x)) / float(len(x)))
  return rr

def signalmatrix(bed, bigwig, expsize = 2000, resolution = 20, s = 1):
  r = [ 0 for _ in range(-expsize, expsize) ]
  regions = []
  with open(bed, 'r') as f:
    for cp, line in enumerate(f):
      l = line.strip().split()
      regions.append( (l[0], recenter((int(l[1]), int(l[2])), expsize)) )
    chunks = [ x for x in chunk(regions, 16) ]
    rr = Parallel(n_jobs = 16)(delayed(bwreadg)(bigwig, c, s) for c in chunks)
    r = flatten(rr)
    return r if resolution == 1 else [ averagem(x, resolution) for x in r ]
