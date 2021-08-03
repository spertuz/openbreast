% Breast density segmentation
clear, clc

% read reference image:
srcdir  = fileparts(mfilename('fullpath'));
impath  = [srcdir,'/../samples/cont_522_RMLO_0.dcm'];
info    = getinfo(impath);
im      = ffdmRead(impath, info);
im      = imresize(im, info.psize/0.4);

% density segmentation using LIBRA
[pd, mask1] = pdensity(impath);
subplot(121)
showseg(im, fliplr(mask1), .25);
title('Density (LIBRA)')

% morphology-based density estimation 
mask = segBreast(im, info.ismlo);
seg = mseg(im, mask, 0.4);
subplot(122), showseg(im, seg, .25)