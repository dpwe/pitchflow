function [wts, frqs] = fft2logfmx(nfft, sr, nbins, width, fmin, bpo)
% [wts, frqs] = ff2logfmx(nfft, sr, nfilts, width, fmin, bpo)
%    Mapping matrix to convert FFT to logf (constant-Q) spectrogram.
% 2014-01-16 dpwe@ee.columbia.edu after logfsgram $Header: $

if nargin < 2;  sr = 8000;  end
if nargin < 3;  nbins = 0; end
if nargin < 4;  width  = 1;  end
if nargin < 5;  fmin = 50; end
if nargin < 6;  bpo = 12; end

% Construct mapping matrix

% Ratio between adjacent frequencies in log-f axis
fratio = 2^(1/bpo);

% How many bins in log-f axis
if nbins == 0
  % default goes all the way to nyquist
  nbins = floor( log((sr/2)/fmin) / log(fratio) );
end

% Freqs corresponding to each bin in FFT
fftfrqs = [0:(nfft/2)]*(sr/nfft);
nfftbins = nfft/2+1;

% Freqs corresponding to each bin in log F output
frqs = fmin * exp(log(2)*[0:(nbins-1)]/bpo);

% Bandwidths of each bin in log F
logfbws = width * frqs * (fratio - 1);

% .. but bandwidth cannot be less than FFT binwidth
logfbws = max(logfbws, sr/nfft);

% Controls how much overlap there is between adjacent bands
ovfctr = 0.5475;   % Adjusted by hand to make sum(mx'*mx) close to 1.0

% Weighting matrix mapping energy in FFT bins to logF bins
% is a set of Gaussian profiles depending on the difference in 
% frequencies, scaled by the bandwidth of that bin
freqdiff = ( repmat(frqs',1,nfftbins) - repmat(fftfrqs,nbins,1) )./repmat(ovfctr*logfbws',1,nfftbins);
wts = exp( -0.5*freqdiff.^2 );
% Normalize rows by sqrt(E), so multiplying by mx' gets approx orig spec back
wts = wts ./ repmat(sqrt(2*sum(wts.^2,2)), 1, nfftbins);
