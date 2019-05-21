% Feature extraction of FFDM images
clear, clc

clear, clc

% read reference image:
srcdir = fileparts(mfilename('fullpath'));
impath = [srcdir,'/../samples/img_143_195_1_LO.tif'];
im = imread(impath);        %read image
im = im2double(im);         %convert to double
imn = imresize(im, .25);     %re-size for speed up
isMLO = true;               %true for MLO view
isFFDM = true;              %true only for FFDM images

% breast segmentation
mask = segBreast(imn, isMLO, isFFDM);
mask = imresize(mask, size(im));

% feature extraction
features = {'imin', 'imax', 'iavg', 'ient', 'istd', 'ip05', 'ip95', 'iba1', 'iba2', 'ip30', 'ip70', 'iske', 'ikur','iran',...
'cene', 'ccor', 'ccon', 'chom', 'cent',...
'rsre', 'rlre', 'rgln', 'rrpe', 'rrln', 'rlgr', 'rhgr',...
'sgra', 'slap', 'swas', 'swav', 'swar', 'stev'};

x = xfeatures(im, features, mask);
disp(table(features(:), x(:), 'variablenames',{'Feature','Value'}))