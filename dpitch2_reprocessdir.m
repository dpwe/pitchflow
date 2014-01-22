function dpitch2_reprocessdir(indir, outdir)
% dpitch2_reprocessdir(indir, outdir)
%   Rerocess an entire directory of dpitch2 features, and collapse
%   to 3-dim moments.
% 2014-01-16 Dan Ellis dpwe@ee.columbia.edu

% Get a list of all the *htk files in <indir>

fdir = dir(fullfile(indir, '*.htk'));

for i = 1:length(fdir) 

  [p,n,e] = fileparts(fdir(i).name);
  inhtk = fullfile(indir, [n, e]);
  outhtk = fullfile(outdir, [n, '.htk']);

  % Process
  [D,FP,DT,TC,T] = readhtk(inhtk);
  dpftr = dpitch2_collapse(D');
  writehtk(outhtk, dpftr', FP, TC);

  disp(['Wrote ', outhtk]);
  
end
