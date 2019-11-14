% Demo on ST-mapping
clear, clc

% read reference image:
srcdir  = fileparts(mfilename('fullpath'));
impath  = [srcdir,'/../samples/cont_522_RMLO_0.dcm'];
info    = getinfo(impath);
im      = ffdmRead(impath, info);
im      = imresize(im, info.psize/0.4);

% Segment breast boundary and chest wall:
[~, contour, cwall] = segBreast(im, info.ismlo);

% Generate and display ST mapping
subplot(1,3,1)
rpoints = refpoints(contour, cwall, info.ismlo); %find reference points
map = stmap(contour, rpoints);              %find mapping parameters
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
% (previously saved as map_ref)
load([srcdir,'/../samples/map_ref.mat'])
im_ref = imread([srcdir,'/../samples/ref.jpg']);
[imr, isvalid] = remapping(im, map, map_ref);
figure
subplot(131), imshow(mat2gray(im)), title('Template')
subplot(132), imshow(mat2gray(im_ref)), title('Reference')
subplot(133), showseg(mat2gray(imr), isvalid), title('Template (re-mapped)')
