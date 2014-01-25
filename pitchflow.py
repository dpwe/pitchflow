#!/usr/bin/env python
"""
pitchflow.py

Python port of the new delta-pitch-without-pitch feature 
based on cross-correlating adjacent spectra.

2014-01-24 Dan Ellis dpwe@ee.columbia.edu
"""

import numpy as np
# For filter
import scipy.signal
# For SRI's wavreading code
import scipy.io.wavfile as wav
# For command line
import sys

def frame(x, window, hop):
    """ Convert vector x into an array of successive window-point segments, 
        stepped by hop
        Done with stride_tricks, no copying, but have to lose the final 
        part-window 
    """

    nframes = 1 + int( np.floor( (len(x)-window)/hop ))
    shape = (nframes, window)
    strides = (x.strides[0] * hop, x.strides[0])
    return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides)

def stftm(signal, nfft, window, hop):
    """ calculate the short-time fourier transform magnitude """
    frames = frame(signal, window, hop)

    # apply frame window to each frame
    window = np.hanning(window)
    wframes = frames * window

    return np.abs(np.fft.rfft(wframes, int(nfft)))

def fft2logfmx(nfft, sr=8000, nbins=0, width=1.0, fmin=50.0, bpo=12.0):
    """
    Mapping matrix to convert FFT to logf (constant-Q) spectrogram.
    2014-01-16 dpwe@ee.columbia.edu after logfsgram
    """
    # Ratio between adjacent frequencies in log-f axis
    fratio = np.power(2.0, (1.0/bpo))

    # How many bins in log-f axis
    if nbins == 0:
        # default goes all the way to nyquist
        nbins = int(np.floor( np.log((sr/2)/fmin) / np.log(fratio) ))
    else:
        nbins = int(nbins)

    # Freqs corresponding to each bin in FFT
    nfftbins = int(nfft/2+1)
    fftfrqs = np.array(range(nfftbins))*(sr/nfft)

    # Freqs corresponding to each bin in log F output
    frqs = fmin * np.exp(np.log(2.0)*np.array(range(nbins))/bpo)

    # Bandwidths of each bin in log F
    logfbws = width * frqs * (fratio - 1.0)

    # .. but bandwidth cannot be less than FFT binwidth
    logfbws = np.maximum(logfbws, sr/nfft)

    # Controls how much overlap there is between adjacent bands
    ovfctr = 0.5475   # Adjusted by hand to make sum(mx'*mx) close to 1.0

    # Weighting matrix mapping energy in FFT bins to logF bins
    # is a set of Gaussian profiles depending on the difference in 
    # frequencies, scaled by the bandwidth of that bin
    freqdiff = ( (np.tile(frqs[:, np.newaxis], (1,nfftbins)) 
                  - np.tile(fftfrqs, (nbins, 1)))
                 / (ovfctr*logfbws[:, np.newaxis]) )
    wts = np.exp( -0.5*np.square(freqdiff) )
    # Normalize rows by sqrt(E), so multiplying by mx' gets aprx orig spec back
    wts = wts / np.sqrt(2*np.sum(np.square(wts)))
    return wts, frqs


def convbyrow(H, X):
    """
    H is a vector, each row of Y is a row of X convolved with H
    """
    (nr, nc) = np.shape(X)
    lh = len(H);
    Y = np.zeros( (nr, nc + lh - 1) )
    for i in range(nr):
        Y[i,] = np.convolve(H, X[i,])
    return Y

def localmvnorm(X,W):
    """ 
    Y = localmvnorm(X,W)
    Remove local mean and divide out local SD from vector X to
    return Y.  W is the length of the local estimation window.
    For X a matrix, normalization is applied to rows.
    2014-01-17 Dan Ellis dpwe@ee.columbia.edu
    """

    (bins, frames) = np.shape(X)

    print 'size(X)=',bins,'bins x',frames,'frames W=',W

    # Pad with reflections at ends for off-the end continuity
    padpoints = np.floor(W/2)
    Xpad = np.r_[ X[padpoints:0:-1,],X,X[:-padpoints-1:-1,] ]

    print np.shape(Xpad)

    (padbins, frames) = np.shape(Xpad)

    win = np.hanning(W+2)
    # strip the zeros
    win = win[1:-1]
    win = win / np.sum(win)
    whlen = int(np.floor((W-1)/2))

    MAmean = convbyrow(win, Xpad.T).T
    print np.shape(MAmean)
    MAmean = MAmean[whlen:whlen+padbins,:]
    MAvar  = convbyrow(win, np.square(Xpad-MAmean).T).T
    print np.shape(MAvar)
    MAmean = MAmean[padpoints:padpoints+bins,:]
    MAstd  = np.sqrt(MAvar[whlen+padpoints : whlen+padpoints+bins,:])

    return (X-MAmean) /MAstd

def xcorr(a, b, halfwidth):
    """ cross-correlate a and b, returning points out to +/- halfwidth """
    r = np.correlate(a,b,'full')
    return r[len(a)-halfwidth-1:len(a)+halfwidth]

########### main function ##############

def pitchflow(d, sr):
    """
    Y = pitchflow(d,sr)
    Simply calculate a spectrogram
    and take cross-correlation of successive spectra
    to capture systematic shift in frequency
    2014-01-16 Dan Ellis dpwe@ee.columbia.edu
    """

    twin = 0.032
    thop = 0.010

    nfft = np.power(2.0, np.round(np.log(twin*sr)/np.log(2.0)))
    nhop = np.round(thop * sr)

    # Calculate base spectrogram
    D = stftm(d, nfft, nfft, nhop)

    # Convert to log-freq axis
    bpo  = 24
    fmin = 50
    fmax = 1500
    nbins = np.round(bpo * np.log(fmax/fmin)/np.log(2.0))
    width = 1.0

    wts, frqs = fft2logfmx(nfft, sr, nbins, width, fmin, bpo)

    # log-f sgram
    print np.shape(wts)
    print np.shape(D)
    DL = np.dot(wts, D.T)

    (nbins, nframes) = np.shape(DL)

    # local mean and variance normalization on each spectrum
    mvnormwin = 48
    DL = localmvnorm(DL, mvnormwin)

    # Record frame-on-frame xcorr

    delay = 2
    halfwidth = 12
    mxmd = np.empty( (2*halfwidth+1, nframes) )
    mxmw = np.empty( nframes )
    for i in range(nframes):
        # output moves with first arg (so right-shift of first arg gives
        # right-shift of output)
        dl0 = DL[:,i]
        dl1 = DL[:, np.maximum(1, i-delay)]
        mxmd[:,i] = xcorr(dl0, dl1, halfwidth)
        # normalizing constant
        mxmw[i] = np.sqrt(np.sum(np.square(dl0)*np.sum(np.square(dl1))))

    # Keep central bin as normalizing constant
    #midbin = halfwidth
    #Y0 = mxmd[midbin,]
    Y = mxmd/mxmw

    return Y

def pitchflow_collapse(Y, NBIN=0):
    """ 
    F = dpitch2_collapse(Y)
    Y is a set of column feature vectors indicating pitchflow from
    dpitch.  Collapse these 25 (?) dimensional vectors into 2 or 3
    summary dimensions - something on peakiness, something on
    center of mass, something on compactness
    2014-01-16 Dan Ellis dpwe@ee.columbia
    """

    (nr, nc) = np.shape(Y)

    maxlag = (nr - 1)/2

    # Maybe chop out just middle bins
    if NBIN > 0:
        newmaxlag = NBIN # was 8
        Y = Y[maxlag-newmaxlag : maxlag+newmaxlag+1,]
        (nr, nc) = np.shape(Y)
        maxlag = newmaxlag

    # Some preprocessing - convolve with smoothing window
    DOSMOOTH = 1
    if DOSMOOTH:
        smoohwin = 2
        smoo = 2*smoohwin+1
        smwin = np.hanning(smoo+2)
        smwin = smwin[1:-1]
        Ys = convbyrow(smwin, Y)
        Ys = Ys[:, smoohwin:-smoohwin];
    else:
        Ys = Y
  
    # Raise to a power to increase dominance of peak values
    #ep = 2.0;
    #Y = Ys.^ep;
    # Now Y is bipolar; best to make it positive before taking moments,
    # and exponentiation also emphasizes peak
    escale = 1.0
    Y = np.exp(Ys/escale)

    # First dimension - "spectral entropy" crest factor
    #p = 2
    #F0 = np.power(np.mean(np.power(Y, p)), 1/p) / np.mean(Y);

    # .. or just first moment
    F0 = np.mean(Y, axis=0)

    # Second dimension - center of mass
    lags = np.array(range(-maxlag,maxlag+1))
    F1 = np.mean(Y.T*lags, axis=1) / F0

    # Third dimension - inertial moment
    F2 = np.sqrt(np.mean(Y.T*np.square(lags), axis=1) / F0 - np.square(F1))

    ftrs = np.c_[F0,F1,F2]
    print np.shape(ftrs)

    return ftrs

############## Provide a command-line wrapper

# from 
# http://my.fit.edu/~vkepuska/ece5526/SPHINX/SphinxTrain/python/sphinx/htkmfc.py
from struct import unpack, pack

class HTKFeat_write(object):
    "Write Sphinx-II format feature files"
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

def writehtk(filename, features):
    """ write array of features to HTK file """
    rows, veclen = np.shape(features)
    writer = HTKFeat_write(filename, veclen)
    writer.writeall(features)

def main(argv):
    """ Main routine to apply pitchflow from command line """
    if len(argv) != 3:
        raise NameError( ("Usage: ", argv[0], 
                          " inputsound.wav outftrs.htk") )

    inwavfile = argv[1]
    outhtkfile = argv[2]

    # Read in wav file
    srate, wavd = wav.read(inwavfile)
    # normalize short ints to floats of -1 / 1
    data = np.asfarray(wavd) / 32768.0  

    # Apply
    ftrs1 = pitchflow(data, srate)
    ftrs = pitchflow_collapse(ftrs1)

    # Write the htk data out
    writehtk(outhtkfile, ftrs)
    print "Wrote ", outhtkfile

# Actually run main
if __name__ == "__main__":
    main(sys.argv)

