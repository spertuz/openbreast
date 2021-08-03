function [si, ti] = xy2st(x, y, MAPP, sc)
% Convert XY coordinates to ST coordinates
% Sintax:
%     [s, t] = xy2st(x, y, mapp)
% Inputs:    
%     x,      Px1 vector with x-coordinates
%     y,      Px1 vector with y-coordinates
%     mapp,   structure with st-mapping of an image
%             (as returned by stmap).
% Outputs:
%     s,      Px1 vector with corresponding s-coordinates
%     t,      Px1 vector with corresponding t-coordinates
%     
% S. Pertuz
% Jul13/2017

if nargin<4, sc = 1;
end

NS = MAPP.size(2);
NT = MAPP.size(1);
s = linspace(0, 1, NS);
t = linspace(0, 1, NT);

% Create sampling grid:
[S, T] = meshgrid(s, t);
theta = (MAPP.theta(2) - MAPP.theta(1))*T + MAPP.theta(1);
PA = T.*polyval(MAPP.P1, S) + polyval(MAPP.P2, S);
PB = T.*polyval(MAPP.P3, S) + polyval(MAPP.P4, S);
X = PA.*cos(theta) - PB.*sin(theta) + MAPP.x0;
Y = PA.*sin(theta) + PB.*cos(theta) + MAPP.y0;

if MAPP.flip
    x = sc*MAPP.size(2)-x+1;
end
warning off
si = griddata(X, Y, S, x/sc, y/sc);
ti = griddata(X, Y, T, x/sc, y/sc);
warning on