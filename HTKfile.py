# HTKfile.py
#
# Provide Matlab interface to HTK files
#
# 2014-01-24 Dan Ellis dpwe@ee.columbia.edu

# from 
# http://my.fit.edu/~vkepuska/ece5526/SPHINX/SphinxTrain/python/sphinx/htkmfc.py

import numpy as np
from struct import unpack, pack

class writer(object):
    """
    Write HTK format feature files
    """
    def __init__(self, filename=None,
                 veclen=13, sampPeriod=100000,
                 paramKind = 9):
        self.veclen = veclen
        self.sampPeriod = sampPeriod
        self.sampSize = veclen * 4
        self.paramKind = paramKind
        self.dtype = 'f'
        self.filesize = 0
        self.swap = (unpack('=i', pack('>i', 42))[0] != 42)
        if (filename != None):
            self.open(filename)

    def __del__(self):
        self.close()

    def open(self, filename):
        self.filename = filename
        self.fh = file(filename, "wb")
        self.writeheader()

    def close(self):
        self.writeheader()

    def writeheader(self):
        self.fh.seek(0,0)
        self.fh.write(pack(">IIHH", self.filesize,
                           self.sampPeriod,
                           self.sampSize,
                           self.paramKind))

    def writevec(self, vec):
        if len(vec) != self.veclen:
            raise Exception("Vector length must be %d" % self.veclen)
        if self.swap:
            np.array(vec, self.dtype).byteswap().tofile(self.fh)
        else:
            np.array(vec, self.dtype).tofile(self.fh)
        self.filesize = self.filesize + self.veclen

    def writeall(self, arr):
        for row in arr:
            self.writevec(row)
