% Demo on Breast Segmentation

% read FFDM image from BCDR dataset:
srcdir = fileparts(mfilename('fullpath'));
impath = [srcdir,'/../samples/img_143_195_1_LO.tif'];
im = imread(impath);        %read image
im = im2double(im);         %convert to double
im = imresize(im, .25);     %re-size for speed up
isMLO = true;               %true for MLO view
isFFDM = true;              %true only for FFDM images

%segment boundary and chest wall:
[mask, contour, cwall] = segBreast(im, isMLO, isFFDM);

%display result
figure
showseg(im, mask)
title('FFDM image')

%detect reference points (nipple is p0)
rpoints = refpoints(contour, cwall, isMLO);
figure, imshow(im), hold on
plot(rpoints.p0.x, rpoints.p0.y, 's', 'markersize', 10)