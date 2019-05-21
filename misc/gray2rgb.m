function im_rgb = gray2rgb(im, cmap, clims)
% Convert grayscale image to RGB
% Sintax:
%     im_rgb = gray2rgb(im);
%     im_rgb = gray2rgb(im, cmap);
%     im_rgb = gray2rgb(im, cmap, clims);
% Inputs:
%     im,     MxN intensity matrix
%     cmap,   Nx3 RGB colormap (default is parula)
%     clims,  Vector [Imin, Imax] with the
%             maximum and minimum intensities
%             to be mapped (default is 
%             [min(im(:)), max(im(:))]
% S. Pertuz
% Apr02/2017

im = double(im);

if nargin<3
    clims = [min(im(:)), max(im(:))];
end

if nargin<2
    cmap = colormap(parula(64));
end


im = (im-clims(1))/(clims(2)-clims(1));
im(im<0)=0;
im(im>1)=1;
ncolors = size(cmap, 1) - 1;
idx = uint8(ncolors*im) + 1;
idx(isnan(idx)) = 1;
im_rgb = ind2rgb(idx, cmap);