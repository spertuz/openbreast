% Breast density segmentation
clear, clc

% read reference image:
srcdir  = fileparts(mfilename('fullpath'));
impath  = [srcdir,'/../samples/cont_522_RMLO_0.dcm'];
info    = getinfo(impath);
im      = ffdmRead(impath, info);
im      = imresize(im, info.psize/0.4);

% morphology-based density estimation 
mask = segBreast(im, info.ismlo);
seg = mseg(im, mask, 0.4);
showseg(im, seg, 0.25);