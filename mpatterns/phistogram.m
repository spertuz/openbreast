function h = phistogram(imdata, P, params)
% Texton histogram of input image
% Sintax:
%     h = phistogram(imdata, P, params)
% Inputs:
%     im,       structure with input image data:
%                 .path,       string with path to images
%                 .cwall,      structure with chest wall data
%                 .contour,    structure with contour data
%                 .class,      boolean flag with type of image
% 
%     P,        NxNP array with MP protorypes. Each colum is 
%               one prototype.
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
%     h,      1xNP histogram with relative frequencies of
%             prototypes from P found in input image.
%             
% S. Pertuz
% Nov25/2017


%Extract MP samples from this image:
S = mpSampling(imdata, params);

%compute distances of every-pixel response to each texton:
nsamples = size(S, 2);
nprototypes = size(P, 2);
d = zeros(nprototypes, nsamples);
for n = 1:nprototypes
    %p = repmat(P(:,n), [1, nsamples]);
    p = P(:,n);
    d(n,:) = sqrt(sum((S-p).^2));
end

%find closest texton to each pixel:
[~, p_index] = min(d);

%generate histogram:
h = zeros(1, nprototypes);
for n = 1:nprototypes
    h(n) = sum(p_index==n);
end
h = h/sum(h);