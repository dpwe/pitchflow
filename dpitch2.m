function Y = dpitch2(d,sr,P,doplot)
% Y = dpitch2(d,sr,P,doplot)
%   Second version of delta-pitch features.
%   Simply calculate a spectrogram
%   and take cross-correlation of successive spectra
%   to capture systematic shift in frequency
% 2014-01-16 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; P = []; end
if nargin < 4; doplot = 0; end

if isfield(P, 't_win'); twin = P.t_win; else twin = 0.032; end
if isfield(P, 't_hop'); thop = P.t_hop; else thop = 0.010; end

nfft = 2^round(log(twin*sr)/log(2));
nhop = round(thop * sr);

% Calculate base spectrogram
D = specgram(d,nfft,sr,nfft, nfft-nhop);

% Convert to log-freq axis
bpo  = 24;
fmin = 50;
fmax = 1500;
nbins = round(bpo * log(fmax/fmin)/log(2));
width = 1.0;

[wts, frqs] = fft2logfmx(nfft, sr, nbins, width, fmin, bpo);

% log-f sgram
DL = wts * abs(D);

[nbins, nframes] = size(DL);

% local mean and variance normalization on each spectrum
mvnormwin = 48;
DL = localmvnorm(DL', mvnormwin)';


% Record frame-on-frame xcorr

delay = 2;
halfwidth = 12;
for i = 1:nframes
  % output moves with first arg (so right-shift of first arg gives
  % right-shift of output)
  dl0 = DL(:,i);
  dl1 = DL(:,max(1, i-delay));
  mxmd(:,i) = xcorr(dl0, dl1, halfwidth);
  % normalizing constant
  mxmw(i) = sqrt(sum(dl0.^2)*sum(dl1.^2));
end

% Keep central bin as normalizing constant
midbin = halfwidth + 1;
Y0 = mxmd(midbin, :);
Y = mxmd./repmat(mxmw, size(mxmd,1), 1);

NORMALIZE = 0;

if NORMALIZE
  % Normalize all bins by central bin
  Y = Y./repmat(Y0, 2*halfwidth+1, 1);
  % Restore the normalizing constant
  Y(halfwidth+1, :) = Y0/max(Y0);
end
