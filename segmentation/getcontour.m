function contour = getcontour(mask, cflag)
% Retrieve breast contour
% Sintax:
%     contour = getcontour(mask, cflag)
% Inputs:
%     mask,       MxN binary segmentation mask
%     cflag,      Binary flag. True for contour
%                 clipping based on curvature. Default 
%                 is false.
% Outputs:
%     contour,    structure with contour data
% 
% S. Pertuz
% Nov06/2017

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_th    = 0.05; % Curvature threshold (higher is more strict)
                % this works only if ismlo = true.
elim  = 8;      % No of pixels to edge limit
npts    = 100;  % number of contour points to return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find contour points:
B = bwboundaries(mask);
remov = (B{1}(:,1) <= elim)|(B{1}(:,2) <= elim)|(B{1}(:,1)>=(size(mask,1)-elim));

ys = B{1}(~remov,1);
xs = B{1}(~remov,2);

% sub-sample and smooth contour:
n = length(xs);
xs = smooth(interp1(xs, linspace(1, n, npts), 'pchip'));
ys = smooth(interp1(ys, linspace(1, n, npts), 'pchip'));

% Crop contour by curvature analysis
[xc, yc, ycut] = cropContour(xs, ys, K_th);

if ~cflag
    contour.x = xs(:);
    contour.y = ys(:);
else
        
    contour.x = xc(:);
    contour.y = yc(:);
end

contour.ycut = ycut;
contour.size = size(mask);

end
function [xc, yc, ycut] = cropContour(xs, ys, K_th, im)
% Crop breast contour using curvature thresholding
% Sintax:
%     [xc, yc] = cropContour(xs, ys, C_th)
% 
% Inputs:
%     xs,     Nx1 vector with the X-coordinates of the breast contour
%     ys,     Nx1 vector with the Y-coordinates of the breast contour
%     Kth,    Curvature threshold
% Ouputs:
%     xc,     Mx1 vector with X-coordinates of cropped contour
%     yc,     Mx1 vector with Y-coordinates of cropped contour
%     
% S. Pertuz
% Jul03/2017

dispflag = (nargin>3);

%compute curvature k:
dx1 = diff(xs);    %dx/dt
dx2 = diff(dx1);    %d2x/dt2
dy1 = diff(ys);    %dy/dt
dy2 = diff(dy1);    %d2y/dt2
dx1(1) = [];
dy1(1) = [];

k = (dx1.*dy2 - dy1.*dx2)./((dx1.^2 + dy1.^2).^(1.5));

% Cut contour points with curvature
% above the threshold K_th

[kmin, i] = min(k);

if (abs(kmin)> K_th)&& (xs(i)<0.4*max(xs))&&(ys(i)>0.5*max(ys))
    ycut = floor(ys(i));
    xs(i+1:end) = [];
    ys(i+1:end) = [];
else
    ycut = floor(max(ys));
end    

xc = xs;
yc = ys;

if dispflag
    figure
    subplot(121)
    imshow(mat2gray(im)), hold on
    plot(xc, yc, 'g.')
    subplot(122)
    plot(k), title('Curvature'), hold on
    line([1 length(k)],[-K_th, -K_th], 'color', 'r')
end
end