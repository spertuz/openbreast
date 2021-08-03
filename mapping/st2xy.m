function [x,y] = st2xy(s, t, MAPP, refsize)
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

if nargin<4
    sc = [1 1];
else
    sc = refsize./MAPP.size;
end

theta = (MAPP.theta(2) - MAPP.theta(1))*t + MAPP.theta(1);
PA = t.*polyval(MAPP.P1, s) + polyval(MAPP.P2, s);
PB = t.*polyval(MAPP.P3, s) + polyval(MAPP.P4, s);
x = round(PA.*cos(theta) - PB.*sin(theta) + MAPP.x0);
y = round(PA.*sin(theta) + PB.*cos(theta) + MAPP.y0);


% re-scale coordinates and flip if necessary
x = round(sc(2)*(x-1) + 1);
y = round(sc(1)*(y-1) + 1);


%Check whether there are point outside the image and replace them
% by the neareast neighbor
no_samples = length(x);
if MAPP.flip
    x = sc(2)*MAPP.size(2)-x+1;
end

if any(x>sc(2)*MAPP.size(2))
    x = min([x(:), sc(2)*MAPP.size(2)*ones(no_samples, 1)], [], 2);
end
if any( x<1 )
    x = max([x(:), ones(no_samples, 1)], [], 2);
end
if any(y>sc(1)*MAPP.size(1))
    y = min([y(:), sc(1)*MAPP.size(1)*ones(no_samples, 1)], [], 2);
end
if any(y<1)
    y = min([y(:), ones(no_samples, 1)], [], 2);
end