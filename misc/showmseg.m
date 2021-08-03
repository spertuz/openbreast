function [ims, mask1, pd] = showmseg(segdata)

% display manual segmentation
impath = segdata.path{1};
im = ffdmRead(impath, getinfo(impath));
mask0 = seg2mask(segdata.contour{1},...
    segdata.cwall{1});
maske = seg2mask(segdata.contour{1},...
    segdata.cwall{1},...
    [], segdata.gap1,...
    segdata.gap2);


mask1 = im>segdata.th1;
mask1 = imresize(mask1, size(maske));
mask1 = mask1&maske;
mask1 = bwareaopen(mask1, segdata.th2);
pd = sum(mask1(:))/sum(mask0(:));

im = adapthisteq(mat2gray(im));
ims = showseg(im, mask0);
ims = showseg(ims, mask1, 0.25);
if nargout==0
    imshow(ims);
end