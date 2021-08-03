% Demo on Breast Segmentation

% read FFDM image :
srcdir = fileparts(mfilename('fullpath'));
impath = [srcdir,'/../samples/cont_522_RMLO_0.dcm'];

% retrieve image info
info = getinfo(impath);

% read image
im = ffdmRead(impath, info);

% standardize size to 0.4mm/pixel
im = imresize(im, info.psize/0.4);

tic
%segment boundary and chest wall:
[mask, contour, cwall] = segBreast(im, info.ismlo);
toc

%display result (image is flipped!)
figure
showseg(fliplr(im), fliplr(mask)), hold on
title('FFDM image')

%detect reference points (nipple is p0)
rpoints = refpoints(contour, cwall, info.ismlo);
plot(rpoints.p0.x, rpoints.p0.y, 's', 'markersize', 10)
legend({'Nipple'})