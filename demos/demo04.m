% Feature extraction of FFDM images
clear, clc

% read reference image:
srcdir  = fileparts(mfilename('fullpath'));
impath  = [srcdir,'/../samples/cont_522_RMLO_0.dcm'];
info    = getinfo(impath);
imn     = ffdmRead(impath, info);
imn     = imresize(imn, info.psize/0.4);

% breast segmentation
mask = segBreast(imn, info.ismlo);

% feature extraction
features = {'imin', 'imax', 'iavg', 'ient', 'istd', 'ip05', 'ip95', 'iba1', 'iba2', 'ip30', 'ip70', 'iske', 'ikur','iran',...
'cene', 'ccor', 'ccon', 'chom', 'cent',...
'rsre', 'rlre', 'rgln', 'rrpe', 'rrln', 'rlgr', 'rhgr',...
'sgra', 'slap', 'swas', 'swav', 'swar', 'stev'};

x = xfeatures(imn, features, mask);

%feature extraction with different sampling, normalization and scales:
norm    = {'none','zscore'};
samp    = {'full', 'multi'};
scal    = {'1.0'};
res     = 0.1;
[xf, info] = ffdmFeatures(impath, res, scal, norm, samp);

fprintf('Extrated features:\n')
disp(features(:));