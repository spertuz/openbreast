% ROI-detection

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
[mask, contour, cwall] = segBreast(im, isMLO, isFFDM);

% Detect maximum squared ROI:
mask_SQ = sqmax(mask);

%Detect multiple-ROIs of size 32x32 pixels:
[x, y] = rcenters(64, 32, mask);
[~, ~, mask_multi] = xySampling(im, x, y, 32);

% Detect retroareolar region
rpoints = refpoints(contour, cwall, isMLO);
map = stmap(contour, rpoints);
mask_RA = st2mask(map, [.1, .5], [.1 .9]);

% show results
subplot(141), showseg(im, mask, .3), title('Full breast')
subplot(142), showseg(im, mask_SQ, .3), title('Squared ROI')
subplot(143), showseg(im, mask_multi, .3), title('Multi-ROI')
subplot(144), showseg(im, mask_RA, .3), title('RA-region')