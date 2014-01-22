function F = dpitch2_collapse(Y, NBIN)
% F = dpitch2_collapse(Y)
%    Y is a set of column feature vectors indicating pitchflow from
%    dpitch.  Collapse these 25 (?) dimensional vectors into 2 or 3
%    summary dimensions - something on peakiness, something on
%    center of mass, something on compactness
% 2014-01-16 Dan Ellis dpwe@ee.columbia

if nargin < 2; NBIN = 0; end

[nr, nc] = size(Y);

maxlag = (nr - 1)/2;

% Maybe chop out just middle bins
if NBIN > 0
  newmaxlag = NBIN;  % was 8

  Y = Y(maxlag+1+[-newmaxlag:newmaxlag], :);

  nr = size(Y,1);
  maxlag = newmaxlag;
end

% Some preprocessing - convolve with smoothing window
DOSMOOTH = 1;
if DOSMOOTH
  smoohwin = 2;
  smoo = 2*smoohwin+1;
  smwin = hann(smoo+2)';
  smwin = smwin(2:end-1);
  Ys = conv2(smwin, Y);
  Ys = Ys(:,1+smoohwin:end-smoohwin);
else
  Ys = Y;
end
  
% Raise to a power to increase dominance of peak values
%ep = 2.0;
%Y = Ys.^ep;
% Now Y is bipolar; best to make it positive before taking moments,
% and exponentiation also emphasizes peak
escale = 1.0;
Y = exp(Ys/escale);

% First dimension - "spectral entropy" crest factor
%p = 2;
%F0 = (mean(Y.^p)).^(1/p) / mean(Y);

% .. or just first moment
F0 = mean(Y);

% Second dimension - center of mass
lags = [-maxlag:maxlag]';
F1 = mean(repmat(lags, 1, nc).*Y) ./ F0;

% Third dimension - inertial moment
F2 = sqrt(mean(repmat(lags.^2, 1, nc).*Y) ./ F0 - F1.^2);

F = [F0;F1;F2];
