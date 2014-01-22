function Y = localmvnorm(X,W)
% Y = localmvnorm(X,W)
%   Remove local mean and divide out local SD from vector X to
%   return Y.  W is the length of the local estimation window.
%   For X a matrix, normalization is applied to rows.
% 2014-01-17 Dan Ellis dpwe@ee.columbia.edu

% Reflect at ends for off-the end

padpoints = floor(W/2);
[Xr, Xc] = size(X);

%disp(['size(X)=',num2str(Xr),'x',num2str(Xc),' W=',num2str(W)]);

Xpad = [X(:, padpoints:-1:1),X,X(:, end:-1:end-(padpoints-1))];

win = hann(W+2);
win = win(2:end-1);
win = win / sum(win);
whlen = floor((W-1)/2);

MAmean = convbyrow(win, Xpad);
MAmean = MAmean(:, whlen+[1:size(Xpad,2)]);
MAvar  = convbyrow(win, (Xpad-MAmean).^2);
MAmean = MAmean(:, padpoints + [1:Xc]);
MAstd  = sqrt(MAvar(:, whlen + padpoints + [1:Xc]));

Y = (X-MAmean)./MAstd;

doplot = 0; %(Xr == 1);

if doplot
  subplot(311)
  tt = 1:Xc;
  plot(tt, X, tt, MAmean, tt, MAmean + MAstd)
  subplot(312)
  plot(tt, Y);
end


%%%%%%%%%%%%%%%%%%%%%%%%
function Y = convbyrow(H, X)
% H is a vector, each row of Y is a row of X convolved with H
[nr, nc] = size(X);
lh = length(H);
Y = zeros(nr, nc + lh - 1);
for i = 1:nr
  Y(i, :) = conv(H, X(i,:));
end

