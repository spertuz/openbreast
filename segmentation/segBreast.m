function [mask, contour, cwall] = segBreast(im, ismlo, isffdm)
% Breast segmentation from mammography image
% Sintax:
%     mask = segBreast(im, imlo)
%     mask = segBreast(im, ismlo, isffdm)
% Inputs:
%     im,     MxN grayscale image as returned by
%             ffdmRead
%     ismlo,  Boolean flag. True is view is MLO (turns chest
%             wall detection on/off).
%     isffdm, Boolean flag. True for FFDM images (default: true)
% 
% Outputs:
%     mask,       MxN binary segmentation mask with
%                 breast region.
%     contour,    structure with contour data.
%     cwall,      structure with chest wall data.
%     
% S. Pertuz
% Nov09/2017

% By default, assume that input is FFDM
if nargin<3 
    isffdm = true;
end

%flip image if necessary
isflipped = isright(im);
if isflipped
    im = fliplr(im);
end

%Breast boundary detection:
if isffdm
    [mask, contour] = ffdmForeground(im, ismlo);
else
    [mask, contour] = sfmForeground(im, ismlo);
end

contour.flip = isflipped;

%Breast chest wall detection:
if ismlo
    [cmask, cwall] = segChest(im, contour);
else
    cmask = true(size(mask));
    cwall = struct('m', size(mask, 1), 'b', 0);
end

mask = mask&cmask;
% mask(1,:) = false;
% mask(:,1) = false;
% mask(end,:) = false;

%flip back mask if necessary
if isflipped
    mask = fliplr(mask);
end