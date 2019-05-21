% Demo on ST-mapping
clear, clc

% read reference image:
srcdir = fileparts(mfilename('fullpath'));
impath = [srcdir,'/../samples/img_143_195_1_LO.tif'];
im = imread(impath);        %read image
im = im2double(im);         %convert to double
im = imresize(im, .25);     %re-size for speed up
isMLO = true;               %true for MLO view
isFFDM = true;              %true only for FFDM images

% Segment breast boundary and chest wall:
[~, contour, cwall] = segBreast(im, isMLO, isFFDM);

% Generate and display ST mapping
subplot(1,3,1)
rpoints = refpoints(contour, cwall, isMLO); %find reference points
map = stmap(contour, rpoints);              %find parameters
showmap(im, map)
title('Mapping')

subplot(1,3,2)
st = mapping(im, map);                      %generate ST-map
imshow(st)
title('ST map')

% Perform inverse mapping
subplot(1,3,3)
[imr, isvalid] = imapping(st, map);
showseg(imr, isvalid)
title('Inverse ST map')

% Re-map template to a reference anatomy
map_ref = map;
im_ref = im;
impath = [srcdir,'/../samples/img_251_334_1_LO.tif'];
im = imresize(im2double(imread(impath)), .25);
[~, contour, cwall] = segBreast(im, isMLO, isFFDM);
rpoints = refpoints(contour, cwall, isMLO);
map = stmap(contour, rpoints);
[imr, isvalid] = remapping(im, map, map_ref);
figure
subplot(131), imshow(mat2gray(im)), title('Template')
subplot(132), imshow(mat2gray(im_ref)), title('Reference')
subplot(133), showseg(mat2gray(imr), isvalid), title('Template (re-mapped)')
