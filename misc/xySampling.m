function [f, r, mask] = xySampling(im, x, y, wsize, dispflag)
% Sample image regions using XY coordinates
% Sintax:
%     [f, r] = xySampling(im, x, y, wsize)
%     [f, r] = xySampling(im, x, y, wsize, dflag)
% 
% Inputs:
%     im,         MxN input grayscale image
%     x,          Px1 vector with the x-coordinates of each ROI
%     y,          Px1 vector with the y-coordinates of each ROI% 
%     dflag,      Optional Boolean flag to display the selected
%                 regions. Default is false.
%                 
% Outputs:
%     r,        A Px1 cell array with each image region
%     f,        A Px1 vector with the intensity values
%               at each x-y position
%                 
% S. Pertuz
% Jul03/2017

dispflag = (nargin>4);

mask = false(size(im));
delta = floor(0.5*wsize);
nrois = length(x);
if nargout<2&&~dispflag
    i = sub2ind(size(im), y(:), x(:));
    f = im(i);
else
    r = cell(nrois, 1);
    f = zeros(nrois, 1);
    for n = 1:nrois        
        r{n} = im(y(n)-delta: y(n)+delta, x(n)-delta:x(n)+delta);        
        f(n) = im(y(n), x(n));
        mask(y(n)-delta: y(n)+delta, x(n)-delta:x(n)+delta) = true;
    end
end


%display sampling
if dispflag
    figure
    imshow(mat2gray(im)), hold on
    nregions = length(r);    
    for n = 1:nregions
        R = [x(n)-delta, y(n)-delta, wsize, wsize];
        rectangle('position', R, 'edgecolor', 'g')
        plot(x(n), y(n), 'r+')
    end
end