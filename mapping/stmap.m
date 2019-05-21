function MAPP = stmap(contour, rpoints)
% Compute parameters of ST mapping
% 
% SINTAX:
%       MAPP = stmap(contour, rp);
%
% DESCRIPTION:
% MAPP = ST_mapping(contour, refPoints) returns a struct that includes the
% parameters of the ST coordinate system such as the coefficients of 4
% polynomials (P1, P2, P3 and P4), the Cartesian coordinates of the origin
% (x0, y0), and the angles of the two extreme points (theta). This coordinate
% system is generated from the contour and the reference points.
%
% INPUTS:
% contour:      is the struct that defines the contour of the breast. It
%               must include the fields 'xpoints' and 'ypoints', which
%               define the Cartesian points of the contour. Also it must
%               include the field: 'Size'
%
% rp:           is a structure that contains the Cartesian coordinates of
%               three reference points: the origin (P0) and two extreme
%               points (P1 and P2). It must include the fields: 'P0', 'P1',
%               and 'P2'.
% im,           optional MxN grayscale matrix with the input image. If this
%               argument is passed, the st-mapping is displayed for 
%               visualization.
%
%
% REFERENCES:
% [1] Pertuz Said, Carme Julia, and Domenec Puig. "A novel mammography image
%     representation framework with application to image registration." 
%     Pattern Recognition (ICPR), 2014 22nd International Conference on. 
%     IEEE, 2014.
%
% Copyright 2017, German F. Torres and Said Pertuz.


% Definitions
D = 5;                                  % Degree of polynomials
xp = contour.x;                         % x-coordinates of the contour
yp = contour.y;                         % y-coordinates of the contour
x0 = rpoints.p0.x;                    % x-coordinates of p0
y0 = rpoints.p0.y;                    % x-coordinates of p1

% Indices of refence points
i0 = find((xp==x0)&(yp==y0));                      
i1 = find((xp==rpoints.p1.x)&(yp==rpoints.p1.y));
i2 = find((xp==rpoints.p2.x)&(yp==rpoints.p2.y));

% Auxiliary curves
x_phy0 = xp(i0:-1:i1);
y_phy0 = yp(i0:-1:i1);
x_phy1 = xp(i0:i2);
y_phy1 = yp(i0:i2);

% Surface distance
s0 = cumsum(sqrt(diff(x_phy0).^2+diff(y_phy0).^2));
s0 = [0;s0(:)]/max(s0);
s1 = cumsum(sqrt(diff(x_phy1).^2+diff(y_phy1).^2));
s1 = [0;s1(:)]/max(s1);

% Fit polynomials:
Px0 = [pfit(s0, x_phy0-x0, D);0];
Px1 = [pfit(s1, x_phy1-x0, D);0];
Py0 = [pfit(s0, y_phy0-y0, D);0];
Py1 = [pfit(s1, y_phy1-y0, D);0];

% Find extreme angles:
theta0 = atan2(polyval(Py0, 1), polyval(Px0, 1))+2*pi;
theta1 = atan2(polyval(Py1, 1), polyval(Px1, 1));

% Coefficients
C0 = cos(-theta0);
C1 = cos(-theta1);
S0 = sin(-theta0);
S1 = sin(-theta1);

% Mapping parameters
MAPP.P1 = C1*Px1 - S1*Py1 - C0*Px0 + S0*Py0;
MAPP.P2 = C0*Px0 - S0*Py0;
MAPP.P3 = S1*Px1 + C1*Py1 - S0*Px0 - C0*Py0;        
MAPP.P4 = S0*Px0 + C0*Py0;
MAPP.x0 = x0;
MAPP.y0 = y0;
MAPP.x1 = rpoints.p1.x;
MAPP.y1 = rpoints.p1.y;
MAPP.x2 = rpoints.p2.x;
MAPP.y2 = rpoints.p2.y;
MAPP.theta = [theta0 theta1];
MAPP.size = contour.size;
MAPP.flip = contour.flip;
end

function P = pfit(x, y, n)
% Custom polynomial Fit.

% Construct Vandermonde matrix V:
V(:,n) = x(:);
for j = n-1:-1:1
    V(:,j) = x(:).*V(:,j+1);
end

% Solve least squares problem:
[Q,R] = qr(V,0);
P = R\(Q'*y);
end