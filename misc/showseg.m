function varargout = showseg(im, mask, alpha)
% Overlay segmentation on image.
% Sintax:
%     imout = showseg(im, mask);
%     imout = showseg(im, mask, alpha)
%
% Inputs:
%     im,     MxN input image
%     mask,   MxN segmentation labels
%     alpha,  alpha value for mask overlay. If this argment
%             is not passed or if alpha=0, only the mask
%             borders are shown
% Ouputs:
%     imout,  MxN image with overlayed 
%             segmentation
% 
% S. Pertuz
% Mar29/2014

%Parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
border = 4;     %border width (in pixels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resize image if necessary
if (size(im,1)~=size(mask,1))||(size(im,2)~=size(mask,2))
    im = imresize(im, size(mask));
end

% convert image to double
if isinteger(im)
    im = double(im);
end

% extract color planes
% gmin = min(im(mask~=0));
% gmax = max(im(mask~=0));
% im = mat2gray(im,  [gmin, gmax]);
im = uint8(255*mat2gray(im));
if ndims(im)==3
    imR = im(:,:,1);
    imG = im(:,:,2);
    imB = im(:,:,3);
else
    imR = im;
    imG = im;
    imB = im;
end

%show segmentation mask:
if (nargin<3)||isempty(alpha)||(alpha==0)
    %show mask edges:    
    maskout = bwperim(mask);
    maskout = imdilate(maskout, ones(border));
    
    
    imR(maskout) = 255;
    imout = cat(3, imR, imG, imB);
else
    %mask overlay
%     imR = im;
%     imG = im;
%     imB = im;
    imG(mask) = imG(mask) + alpha*255;
    imout = cat(3, imR, imG, imB);
end

% generate output
if nargout==1
    varargout{1} = imout;
else
    imshow(imout)
end
end