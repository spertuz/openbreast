function [im, isvalid] = imapping(st, mapp, s, t)
% Inverse ST-mapping
% SINTAX:
%     im = imapping(st, mapp)
%     im = imapping(st, mapp)
% Inputs:
%     st,     NTxNS matrix with st representation of mammogram
%     mapp,   structure with parameters of the ST-coordinate
%             system (as returned by stmap)
%     s,      Optional NSx1 vector with s-coordinates of ROI
%     t,      Optional NTx1 vector with t-coordinates of ROI
%    
% Outputs:
%     im,     M x N output image where Of is the
%             oversampling factor of mapp
% 
% Said Pertuz
% Mar13/2013


NY = mapp.size(1);   %No. of points in Y-grid
NX = mapp.size(2);   %No. of points in X-grid

% Compute sample coordinates:
if nargin<3
    s = linspace(0, 1, size(st,2));
    t = linspace(0, 1, size(st,1));
end

[S, T] = meshgrid(s, t);
theta = (mapp.theta(2) - mapp.theta(1))*T + mapp.theta(1);
PA = T.*polyval(mapp.P1, S) + polyval(mapp.P2, S);
PB = T.*polyval(mapp.P3, S) + polyval(mapp.P4, S);
x = PA.*cos(theta) - PB.*sin(theta) + mapp.x0;
y = PA.*sin(theta) + PB.*cos(theta) + mapp.y0;

% Create sampling grid
xi = 1:NX;
yi = 1:NY;
[XI, YI] = meshgrid(xi, yi);

%Compute image using linear interpolation:
warning('off', 'all');
im = griddata(x(:), y(:), st(:), XI, YI, 'linear');
isvalid = st2mask(mapp, [0,1], [0,1]);

if mapp.flip, im = fliplr(im);
end
im(~isvalid|isnan(im)) = 0;