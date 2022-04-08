#!/usr/bin/env python3

import os
import tempfile

class FASTAFile:

    def __init__(self, bed, twoBit):
        self.regions = bed
        self.twoBit = twoBit
    
    def __enter__(self):
        self.tempfile = tempfile.NamedTemporaryFile('rt')
        with tempfile.NamedTemporaryFile('wt') as o:
            with open(self.regions, 'r') as ff:
                for line in ff:
                    o.write('\t'.join([
                        "%s\t%s\t%s\t%s:%s-%s\n" % (line.split()[0], line.split()[1], line.split()[2], line.split()[0], line.split()[1], line.split()[2])
                        for line in ff
                    ]))
            o.flush()
            os.system("twoBitToFa {twoBit} -bed={bed} {fa}".format(
                twoBit = self.twoBit,
                bed = o.name,
                fa = self.tempfile.name
            ))
        return self.tempfile
    
    def __exit__(self, *args):
        self.tempfile.close()
