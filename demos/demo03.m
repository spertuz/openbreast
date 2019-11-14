% ROI-detection

clear, clc

% read reference image:
srcdir  = fileparts(mfilename('fullpath'));
impath  = [srcdir,'/../samples/cont_522_RMLO_0.dcm'];
info    = getinfo(impath);
im      = ffdmRead(impath, info);
im      = imresize(im, info.psize/0.4);

% Segment breast boundary and chest wall:
[mask, contour, cwall] = segBreast(im, info.ismlo);

% Detect maximum squared ROI:
mask_SQ = sqmax(mask);

%Detect multiple-ROIs of size 32x32 pixels:
[x, y] = rcenters(64, 32, mask);
[~, ~, mask_multi] = xySampling(im, x, y, 32);

% Detect retroareolar region
rpoints = refpoints(contour, cwall, info.ismlo);
map = stmap(contour, rpoints);
mask_RA = st2mask(map, [.1, .5], [.1 .9]);

% show results
subplot(141), showseg(im, mask, .3), title('Full breast')
subplot(142), showseg(im, mask_SQ, .3), title('Squared ROI')
subplot(143), showseg(im, mask_multi, .3), title('Multi-ROI')
subplot(144), showseg(im, mask_RA, .3), title('RA-region')