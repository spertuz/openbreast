function [x, y] = rcenters(spacing, wsize, mask)
% Get region centers in XY coordinates
% Sintax:
%     [x, y] = rcenters(spacing, wsize, mask);    
% inputs:
%     spacing,    an integer Sx with the spacings of the
%                 regions in the x- and y- directions.
%     wsize,      window size of each region
%     mask,       an optional binary mask with the valid
%                 breast area. Regions outside the mask
%                 will not be considered.
%                 
% Outputs:
%     x,          a Px1 vector with the x-coordinates of
%                 the ROIs
%     y,          a Px1 vector with the y-coordinates of
%                 the ROIs
%                 
% S. Pertuz
% Jul13/2017

% Find coordinates of samples:
[m, n] = size(mask);
xs = round(1:spacing:n);
ys = round(1:spacing:m);
[xs, ys] = meshgrid(xs, ys);

% Remove regions with pixels outside the image:
delta = floor(0.5*wsize);
xmin = xs(:) - delta;               
xmax = xs(:) + delta;
ymin = ys(:) - delta;
ymax = ys(:) + delta;
remov = (xmin<1)|(xmax>n)|(ymin<1)|(ymax>m);
xs(remov) = [];
ys(remov) = [];


% remove regions with pixels outside the mask. All corners of the ROI
% must be inside mask.
i1 = sub2ind([m, n], ymin(~remov), xmin(~remov));   % upper left corner
i2 = sub2ind([m, n], ymin(~remov), xmax(~remov));   % upper right corner
i3 = sub2ind([m, n], ymax(~remov), xmin(~remov));   % lower left corner
i4 = sub2ind([m, n], ymax(~remov), xmax(~remov));   % lower right corner
remov = ~(mask(i1)&mask(i2)&mask(i3)&mask(i4));
x = xs(~remov);
y = ys(~remov);

if ~(isempty(x)||isempty(y))
    return
else
    %If the image is two small use only one window
    %located at the center of the breast
    stats = regionprops(mask, 'centroid');
    x = min([round(stats.Centroid(1)), n-delta]);
    y = min([round(stats.Centroid(2)), m-delta]);
    x = max([x, 1+delta]);
    y = max([y, 1+delta]);
end




% figure, imshow(mask), hold on, plot(x, y, 'r+')