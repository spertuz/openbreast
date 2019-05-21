function mask = roi2mask(xpoints, ypoints, imsize)
% Covert ROI to binary mask
% Sintax:
%     mask = roi2mask(xpoints, ypoints, imsize)
% Inputs:
%     xpoints,    Nx1 array with x-coordinates of ROI edges
%     ypoints,    Nx1 array with y-coordinates of ROI edges
%     imsize,     1x2 array with size of output mask
% Outputs:
%     mask,       MxN binary mask with selected ROI
%     
% S. Pertuz
% Jan09/2018
dummy = 0;
mask = poly2mask(xpoints, ypoints, imsize(1), imsize(2));
mask = imdilate(mask, ones(9));

