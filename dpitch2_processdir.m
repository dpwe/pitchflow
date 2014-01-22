function dpitch2_processdir(indir, outdir)
% dpitch2_processdir(indir, outdir)
%   Process an entire directory of wavfiles and write out dpitch2
%   feature files.
% 2014-01-16 Dan Ellis dpwe@ee.columbia.edu

% Get a list of all the *sph files in <indir>

fdir = dir(fullfile(indir, '*.sph'));
for i = 1:length(fdir) 

  [p,n,e] = fileparts(fdir(i).name);
  inwav = fullfile(indir, [n, e]);
  outhtk = fullfile(outdir, [n, '.htk']);
  dpitch2_wav2htk(inwav, outhtk);

  disp(['Wrote ', outhtk]);
  
end
