function  [mask, cwall] = segChest(im, contour)
% Chest wall detection in FFDM image
% Sintax:
%     [mask, cwall] = ffdmChest(im, mask)
% 
% Inputs:
%     im,     MxN grayscaly mammography image
%     c,      Px2 vector with x-y coordinates
%             of breast contour (as returned
%             by the function ffdmForeground)
% Outputs:
%     mask,   MxN binary mask with the segmented
%             chest wall
%     cwall,  structure with the slope and intersect of the
%             pectoral line
% 
% S. Pertuz
% Jul03/2017


%Pre-process image:
imc = imdilate(im, strel('disk',8,0));
imc = imfilter(imc, fspecial('gaussian',[5 5],1));

%%%% Detect pectoral line using Hough transform %%%%
%Find edges:
edge_map = edge(imc, 'canny', [], 2.0);

% Crop edge map according to contour:
% ROI1:
% xmax = round(max(contour.x));
% ymax = round(max(contour.y));
% imc = im(1:ymax, 1:xmax);
% ROI2:

ymax = round(0.6*max(contour.y));
xmax = round(min(contour.x(contour.y<ymax)));
edge_map = edge_map(1:ymax, 1:xmax);

%remove lower diagonal:
[m,n] = size(edge_map);
[x, y] = meshgrid(1:n, 1:m);
yref = m - (m - 1)*(x-1)/(n-1);
edge_map(y>yref) = 0;

% %detect chest wall using hough transform:
% theta = linspace(25, 85, 256);
% res   = 0.5*norm(size(edge_map))/256;
% [H, T, R] = hough(edge_map, ...
%     'rhoResolution', res, 'theta', theta);
% P = houghpeaks(H, 1);
% line = houghlines(edge_map, T, R, P);
% b = line.rho/sind(line.theta);
% m = -cosd(line.theta)/sind(line.theta);

% Find edges coordinates:
[y, x] = find(edge_map);
%Quantize parameters space:
N = 128; %Number of quantization points
% rho_max = sqrt(sum(size(im).^2));
rho_max = min(size(edge_map));
rho_min = 1;
theta_min = 20*pi/180;
theta_max = 45*pi/180;

%Compute accumulation array A
theta = linspace(theta_min, theta_max, N);
rho = linspace(rho_min, rho_max, N);
rho_k = x*cos(theta) + y*sin(theta);
A = histc(rho_k, rho);

%Find maximum in accumulation array:
[~,imax] = max(A(:));
[i,j] = ind2sub(size(A), imax);

%Get rect
T = theta(j);
R = (rho_max - rho_min)*(i - 1)/(N-1) + rho_min;
b = R/sin(T);
m = -cos(T)/sin(T);

%Get mask
[x, y] = meshgrid(1:size(im,2), 1:size(im,1));
mask = true(size(im));
mask(y<b+m*x) = false;
cwall.m = m;
cwall.b = b;
end