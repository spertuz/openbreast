function [M, C] = mpSampling(imdata, params)
% Random micro-pattern sampling of mammography image
% Sintax:
%     [M, C] = psampling(imdata, params)
% 
% Inputs:
%     imdata,   1xK structure with fields:
%               .path,       string with path to images
%               .cwall,      structure with chest wall data
%               .contour,    structure with contour data
%               .class,      boolean flag with type of image
% 
%     params,   a structure with the following fields:  
%               .nsamples,  integer with the number of micro-pattern samples 
%                           per image.
%               .mpsize,    integer with micro-pattern size in pixels (odd 
%                           number).
%               .rflag,     boolean flag to apply feature extraction in the
%                           RA region. Default is false.
%               .nflag,     boolean flag to perform feature normalization.
%                           Default is false.
% Outputs:
%     M,        array of size (mpsize^2)x K*nsamples with extracted patches
%     C,        array of size 1x K*nsamples with corresponding class
% 
% S. Pertuz
% Nov24/2017

resolution = 0.05;
    
% Apply recursively %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_images = size(imdata, 1);
if no_images>1
    M = zeros(params.mpsize^2, params.nsamples*no_images);
    C = false(1, params.nsamples*no_images);
    fprintf('Sampling images     ')
    for n = 1:no_images
        m0 = (n-1)*params.nsamples + 1;
        m1 = n*params.nsamples;
        [M(:,m0:m1), C(m0:m1)] = mpSampling(imdata(n,:), params);
        fprintf('\b\b\b\b\b[%02d%%]', floor(100*n/no_images))
    end
    fprintf('\n')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if istable(imdata)
    imdata = table2struct(imdata);
end

% load image
info = getinfo(imdata.path);
im = ffdmRead(imdata.path, info);
im = imresize(im, info.psize/resolution);

%load mask
[mask, ~] = seg2mask(imdata.contour, imdata.cwall);
mask = imresize(mask, size(im));

% %remove cancer region (if any)
% if imdata.class
%     maskr = roi2mask(imdata.xpoints, imdata.ypoints, size(im));    
%     mask = mask&~maskr;
% end

mask = imerode(mask, ones(params.mpsize+1));

% find mask indici
idx = find(mask(:));

% select random locations
rng(1);
i = idx(randperm(length(idx), params.nsamples));

% initialize output:
M = zeros(params.mpsize^2, params.nsamples);

% sample micro-patterns:
[y, x] = ind2sub(size(mask), i);
delta = floor(0.5*params.mpsize);
for n = 1:length(x)
    tmp = im(y(n)-delta:y(n)+delta, x(n)-delta:x(n)+delta);   
    M(:,n) = tmp(:);
end

% apply Webber's law if requested
if params.nflag 
    C = repmat(sqrt(sum(M.^2)), [size(M, 1), 1]);
    M = M.*log(1+C/0.03)./(C+eps);
end

% output labels
C = repmat(imdata.class, [1, size(M,2)]);