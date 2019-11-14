function [seg, pd, mask] = mseg(im, mask, pixel_size)
% Morphological PD estimation
% Sintax:
%     [seg, pd, mask, im] = mseg(im, pixel_size)
% Inputs:
%     im,           input image
%     mask,         breast mask
%     pixel_size,   pixel size in mm
% 
% Outputs:
%     seg,        segmented dense tissue
%     pd,         percent density
% 
% S. Pertuz & F. Torres
% Nov/2019

% Segmentation parameters:
sigma   = 0.8;  %std of smoothing gaussian filter
n       = 25;   %number of gray levels

% Post-processing parameters
skin_gap    = 8;    % skin gap in mm
area_th     = 16;   % area threshold in mm^2


% Smoothing filter
h = gausswin(5*sigma+1);

% Intensity levels
imin = min(im(mask));
imax = max(im(mask));
ivalues = linspace(imin, imax, n);

% Compute morphological area curve
area = zeros(n, 1);
for k = 1:n
    seg = (im>=ivalues(k))&mask;
    area(k) = sum(seg(:));
end
 
% Compute first morphological area gradient (MAG)
mag = diff(area);

% Smooth MAG to remove noise
mag = conv(mag, h, 'same');
 
% Minimize MAG
[~, i] = min(mag);

% Segment image
seg = (im>=ivalues(i+1))&mask;

% Post-process image
seg = pprocess(seg, mask, skin_gap, area_th, pixel_size);

% Percent density
pd = sum(seg(:))/sum(mask(:));

