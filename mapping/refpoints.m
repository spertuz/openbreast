function rpoints = refpoints(contour, cwall, ismlo)
% Localization of the three reference points (the nipple and two extreme
% points)
% 
% SINTAX:
%       rp = refPoints(contour, pecLine);
%
% DESCRIPTION:
% rp = refPoints(Img, contour) returns a struct that includes
% the Cartesian coordinates for three reference points as fields: p0 
% (the nipple), p1 and p2 (superior and inferior extremes). The nipple 
% corresponds to the farthest point from the pectoral line, which is
% identified using the straight line Hough transform. The
% extreme points are located at normalized distances, 
% measured from p0.
%
% INPUTS:
% 
% contour:      is the struct that defines the contour of the breast. It
%               must include the fields 'x' and 'y' with the
%               Cartesian coordinates of the contour points
% pecLine,      is a structure with the slope 'm' and intercept 'b'
%               of the pectoralis line, as returned by the functions
%               'ffdmChest' or 'sfmChest'.
%
% 
% REFERENCES:
% [1] Torres, German F., and Pertuz, S. "Automatic Detection of the Retroareolar
%     Region in X-Ray Mammography Images." In VII Latin American Congress on 
%     Biomedical Engineering CLAIB 2016, Bucaramanga, Santander, Colombia,
%     October 26th-28th, 2016, pp. 157-160. Springer, Singapore, 2017.
%
% Copyright 2017, German F. Torres and Said Pertuz.

%---------------- Assign Defaults ----------------
% Parameters of the extreme point (Results of Statistical investigation)
if nargin<3
    ismlo = false;
end
if ~ismlo
    rsup = 0.4800;
    rinf = 0.4800;
else
    rsup = 0.4573;
    rinf = 0.3500;
end

% Find the nipple (P0)
nipple = find_nipple(contour, cwall);

% Find the extreme points (P1 and P2)
[supPoint, infPoint] = extremePoints(nipple, contour, rsup, rinf);

% Output struct
rpoints = struct('p0', nipple, 'p1', supPoint, 'p2', infPoint);

end

function [supPoint, infPoint] = extremePoints(nipple, contour, rsup, rinf)
% Function for the localization of the extreme points

% Cartesian coordinates of the points on the contour
xc = contour.x;
yc = contour.y;

% Index of the nipple point in the contour curve
in = find( (xc==nipple.x)&(yc==nipple.y));  

% Computation of the auxiliar curves phy0 and phy1
x_phy0 = xc(in:-1:1);
y_phy0 = yc(in:-1:1);
x_phy1 = xc(in:end);
y_phy1 = yc(in:end);

% Calculating the length of the curves phy0 and phy1
r0 = cumsum(sqrt(diff(x_phy0).^2+diff(y_phy0).^2));
r1 = cumsum(sqrt(diff(x_phy1).^2+diff(y_phy1).^2));

% Total lenght of the contour
rmax = max(cumsum(sqrt(diff(xc).^2+diff(yc).^2)));

% Normalized curve
r = [-r0(end:-1:1);0;r1]/rmax;

% Finding the superior point
k1 = dsearchn(r, -rsup);
xsup = xc(k1);
ysup = yc(k1);

% Finding the inferior point
k0 = dsearchn(r, rinf);
xinf = xc(k0);
yinf = yc(k0);

% Output structs
supPoint.x = xsup;
supPoint.y = ysup;
infPoint.x = xinf;
infPoint.y = yinf;   

end

function nipple = find_nipple(contour, pecLine)
% Detection of nipple point on the contour

% Cartesian coordinates of the points on the contour
xc = contour.x;
yc = contour.y;

% Detection of the pectoral line
% [m, b, mask] = pecLineDetection(I0, contour);
m = pecLine.m;
b = pecLine.b;

% Distance from each point of the contour to the pectoral line
p_distance = abs(m*xc-yc+b)/sqrt(m^2+1);
p_distance(yc>0.8*max(yc)) = -inf;
% Index of the maximum distance
[~, id] = max(p_distance);

% Coordinates of the nipple point
nipple.x = xc(id);
nipple.y = yc(id);

% Pectoral line parameters
% pecLine.mask = mask;
% pecLine.m = m;
% pecLine.b = b;   
end