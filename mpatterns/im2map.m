function [map0,map1,immap] = im2map(imdata, model, w)
% Convert input mammogram to MP map
% Sintax:
%     map = im2map(imdata, model, DH)
% Inputs:
%     imdata, data structure with input mammogram
%     model,  MP model     
%     w,      1xN vector with prototype weights
%     
% Output:
%     map,    HxW output map
%     
% S. Pertuz
% Feb 16/2017

%maximum number of samples per batch. This
% is necessary due to memmory limitations:
nsamples = 5e5; %number of samples per batch
alpha = 0.3;    %color saturation 
sgap = 8;       %skin gap (mm)

% read image
info = getinfo(imdata.path);
im = ffdmRead(imdata.path, info);

%load mask
if model.params.rflag
    [~, mask] = seg2mask(imdata.contour, imdata.cwall, imdata.rpoints);
else
    [mask, ~] = seg2mask(imdata.contour, imdata.cwall);
end
mask = imresize(mask, size(im));
sgap = round(sgap/info.psize);
mask = imerode(mask, ones(sgap));

% find mask indici
idx_N0 = find(mask(:));
prototype_N0 = zeros(size(idx_N0));

% initialize ouput
map = zeros(size(mask));

nbatches = ceil(length(idx_N0)/nsamples);
for k = 1:nbatches
    fprintf('Batch %d/%d\n', k, nbatches)
    n0 = (k-1)*nsamples + 1;
    n1 = min([k*nsamples, length(idx_N0)]);
    
    %select indici for this batch
    idx = idx_N0(n0:n1);
    
    % initialize output:
    S = zeros(model.params.mpsize^2, length(idx));
    
    % sample micro-patterns:
    [y, x] = ind2sub(size(mask), idx);
    delta = floor(0.5*model.params.mpsize);
    fprintf('Sampling       ')
    for n = 1:length(idx)
        tmp = im(y(n)-delta:y(n)+delta, x(n)-delta:x(n)+delta);
        map(y(n)-delta:y(n)+delta, x(n)-delta:x(n)+delta) = n0+n-1;
        S(:,n) = tmp(:);
        fprintf('\b\b\b\b\b[%02d%%]', floor(100*n/length(x)))
    end
    
    % apply Webber's law if requested
    if model.params.nflag
        fprintf('\nNormalize ')
        c = sqrt(sum(S.^2));
        S = S.*log(1+c/0.03)./(c+eps);
    end
    
    %compute distances of every-pixel response to each texton:
    nprototypes = size(model.P, 2);
    d = zeros(nprototypes, length(idx));
    fprintf('[100%%]\nDistances      ')
    for n = 1:nprototypes    
        p = model.P(:,n);
        d(n,:) = sqrt(sum((S-p).^2));
        fprintf('\b\b\b\b\b[%02d%%]', floor(100*n/nprototypes))
    end
    
    %find closest texton to each pixel:
    [~, prototype_N0(n0:n1)] = min(d);
    fprintf('\n')
end
map(map~=0) = prototype_N0(map(map~=0));
map(map~=0) = w(map(map~=0));
map0 = map.*double(map>0);
map1 = abs(map).*double(map<=0);
immap = cat(3, im, im, im);
immap(:,:,1) = immap(:,:,1) + alpha*mat2gray(map1);
immap(:,:,2) = immap(:,:,2) + alpha*mat2gray(map0);