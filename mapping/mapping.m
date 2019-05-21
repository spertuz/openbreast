function [st, S, T] = mapping(im, MAPP, s, t)
% ST mapping of input image
% SINTAX:
%     st = mapping(im, mapp)
%     st = mapping(im, mapp)
%     st = mapping(im, mapp, s, t)
% 
% Inputs:
%     im,     MxN grayscale input image
%     mapp,   structure with parameters of ST-cordinate
%             system as returned by stmap.
%     s,      NSx1 vector with s-coordinates of region of
%             interest. Default is linspace(0, 1, Of*N).
%     t,      NTx1 vector with t-coordinates of region of
%             interest. Default is linspace(0, 1, Of*M).
%             Of is an oversampling factor of 4.
% 
% Outputs:
%     st,     P x P output ST representation
%             of im. 
% Said Pertuz
% Mar13/2013


if isinteger(im), im = im2double(im); 
end

if MAPP.flip, im = fliplr(im);
end

if nargin<3
    NS = size(im, 2);
    NT = size(im, 1);
    s = linspace(0, 1, NS);
    t = linspace(0, 1, NT);
end

% Create sampling grid:
[S, T] = meshgrid(s, t);
theta = (MAPP.theta(2) - MAPP.theta(1))*T + MAPP.theta(1);
PA = T.*polyval(MAPP.P1, S) + polyval(MAPP.P2, S);
PB = T.*polyval(MAPP.P3, S) + polyval(MAPP.P4, S);
Xi = PA.*cos(theta) - PB.*sin(theta) + MAPP.x0;
Yi = PA.*sin(theta) + PB.*cos(theta) + MAPP.y0;

%Sample with splines:
st = interp2(im, Xi, Yi);
end