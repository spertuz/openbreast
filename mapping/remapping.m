function [Ir, isvalid] = remapping(im, MAPP, RMAPP)
% Remap mammography image using reference ST map (implicit
% image registration)
% 
% Sintax:
%   Ir = remapping(im, mapp, map_ref);

% Inputs:
%     im,     MxN input image
%     mapp,   structure with ST mapp of input image
%     rmapp,  structure with reference ST map for remapping
% 
% Outputs:
%     Ir,     output image remapped to rmapp
%     
% Said Pertuz
% Mar15/2013


%Find (x1,y1) coordinates in ST-space of input image:


NR = size(im, 2);   %Number of points (R-grid)
NT = size(im, 1);   %Number of points (T-grid)
s = linspace(0, 1, NR);
t = linspace(0, 1, NT);
[S,T] = meshgrid(s, t);

theta = (MAPP.theta(2) - MAPP.theta(1))*T + MAPP.theta(1);
PA = T.*polyval(MAPP.P1, S) + polyval(MAPP.P2, S);
PB = T.*polyval(MAPP.P3, S) + polyval(MAPP.P4, S);
x1 = round(PA.*cos(theta) - PB.*sin(theta) + MAPP.x0);
y1 = round(PA.*sin(theta) + PB.*cos(theta) + MAPP.y0);

%Find target (x2,y2) coordinates in reference ST-space:
theta = (RMAPP.theta(2) - RMAPP.theta(1))*T + RMAPP.theta(1);
PA = T.*polyval(RMAPP.P1, S) + polyval(RMAPP.P2, S);
PB = T.*polyval(RMAPP.P3, S) + polyval(RMAPP.P4, S);
x2 = round(PA.*cos(theta) - PB.*sin(theta) + RMAPP.x0);
y2 = round(PA.*sin(theta) + PB.*cos(theta) + RMAPP.y0);

%Remove points out of the image limits
Xmax = size(im, 2);
Ymax = size(im, 1);
remov1 = (x1>Xmax)|(x1<1)|(y1>Ymax)|(y1<1);
remov2 = (x2>Xmax)|(x2<1)|(y2>Ymax)|(y2<1);
remov = remov1|remov2;
x1(remov) = [];
x2(remov) = [];
y1(remov) = [];
y2(remov) = [];
Ir = zeros(size(im));
index2 = sub2ind([Ymax, Xmax], y2, x2);
index1 = sub2ind([Ymax, Xmax], y1, x1);

%Assign corresponding points: (x1,y1) --> (x2,y2)
if MAPP.flip
    im = fliplr(im);
end

Ir(index2) = im(index1);
isvalid = st2mask(RMAPP, [0,1], [0,1]);

if RMAPP.flip
    Ir = fliplr(Ir);
    isvalid = fliplr(isvalid);
end