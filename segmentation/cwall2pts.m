function pts = cwall2pts(cwall)
% Convert chest wall to points
% Sintax:
%     pts = cwall2pts(cwall)
% Inputs:
%     cwall,      a structure with intercept and slope of
%                 chest wall (as returned by segBreast).
% Outputs:
%     pts,        100x2 array with x-y coordinates of points
%                 on the chest wall.
%                 
% S. Pertuz
% Jan09/2018

%chest wall points:
x = linspace(1, floor(-cwall.b/cwall.m), 100);
y = cwall.b + cwall.m*x;
x(1) = 0;
y(end) = 0;

pts = [x(:), y(:)];