function [x,y] = st2xy(s, t, MAPP, sc)
% Convert ST-cordinates to XY coordinates
% Sintax:
%     [x, y] = st2xy(s, t, mapp)
% Inputs:
%     s,      Px1 vector with s-coordinates
%     t,      Px1 vector with x-coordinates
%     mapp,   structure with st-mapping of an image
%             (as returned by stmap).
% Outputs:
%     x,      Px1 vector with corresponding x-coordinates
%     y,      Px1 vector with corresponding y-coordinates
%     
% S. Pertuz
% Jul13/2017

if nargin<4, sc = 1;
end

theta = (MAPP.theta(2) - MAPP.theta(1))*t + MAPP.theta(1);
PA = t.*polyval(MAPP.P1, s) + polyval(MAPP.P2, s);
PB = t.*polyval(MAPP.P3, s) + polyval(MAPP.P4, s);
x = round(PA.*cos(theta) - PB.*sin(theta) + MAPP.x0);
y = round(PA.*sin(theta) + PB.*cos(theta) + MAPP.y0);

%Check whether there are point outside the image and replace them
% by the neareast neighbor

x = sc*x;
y = sc*y;

no_samples = length(x);
if any(x>sc*MAPP.Size(2))
    x = min([x(:), sc*MAPP.Size(2)*ones(no_samples, 1)], [], 2);
end
if any( x<1 )
    x = max([x(:), ones(no_samples, 1)], [], 2);
end
if any(y>sc*MAPP.Size(1))
    y = min([y(:), sc*MAPP.Size(1)*ones(no_samples, 1)], [], 2);
end
if any(y<1)
    y = min([y(:), ones(no_samples, 1)], [], 2);
end