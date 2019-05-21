function [f, full_list] = xfeatures(im, fname, mask, psize)
% Extract image features
% Sintax:
%     f = xfeatures(im, fname, mask)
%     f = xfeatures(im, fname, mask, psize)
%     f = xfeatures(im, fname, [], psize)
% Inputs:
%     im,     MxN grayscale mammographic image region or a Px1 cell array
%             with image regions (as returned by xySampling).
%     fname,  A string (or Qx1 cell array of strings)
%             with the feature(s) to compute.
%     mask,   an optional MxN binary mask with the region
%             of interest where the features will be 
%             computed. If this argument is not passed,
%             the whole image is processed (default). If im is
%             a cell array, this argument is ignored.
%     psize,  a double with the pixel size. This argument is required only
%             for feature extraction using 'fdim'.
% 
% Outputs:
%     f,      PxQ feature vector
%     
% S. Pertuz
% Jul10/2017

% Note: this is just a wrapper for efficiecy. This function simply groups
% features according to their working principles into four groups: co-occurence
% matrix (GLCM), run-length matrix (GLRL), histogram analysis (GLHA) and
% shaprness measure (GLSM). Each group of features is then computed using
% specific funcions.


% Compute iteratively for multiple regions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(im)
    P = length(im);
    Q = length(fname);
    f = zeros(P, Q);
    for p = 1:P
        f(p,:) = xfeatures(im{p}, fname, [], psize);
    end
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<3)||isempty(mask), mask = true(size(im));
end

if any(size(mask)~=size(im))
    mask = imresize(mask, size(im), 'nearest');
end

%intensity-based features
list1 = {'imin', 'imax', 'iavg',...
    'ient', 'istd', 'ip05',...
    'ip95', 'iba1', 'iba2',...
    'ip30', 'ip70', 'iske',...
    'ikur', 'iran'};

%Co-occurrence matrix
list2 = {'cene', 'ccor', 'ccon',...
    'chom', 'cent'};

%Run length matrix
list3 = {'rSRE', 'rLRE', 'rGLN',...
    'rRPE', 'rRLN', 'rLGR', 'rHGR'};

%Sharpness measures
list4 = {'sgra', 'slap', 'swas', 'swav', 'swar', 'stev'};

% %other structural features:
list5 = {'fdim'};

[flist1, i1] = feature_list(list1, fname);
[flist2, i2] = feature_list(list2, fname);
[flist3, i3] = feature_list(list3, fname);
[flist4, i4] = feature_list(list4, fname);
[flist5, i5] = feature_list(list5, fname);

f1 = [];
f2 = [];
f3 = [];
f4 = [];
f5 = [];

if ~isempty(flist1) % Graylevel histogram analysis
    f1 = features_GLHA(im, flist1, mask);
end

if ~isempty(flist2) %Graylevel co-occurrence matrix
    f2 = features_GLCM(im, flist2, mask);
end

if ~isempty(flist3) %Graylevel run-length
    f3 = features_GLRL(im, flist3, mask);
end

if ~isempty(flist4) %Gray-level sharpness measures
    f4 = features_GLSM(im, flist4, mask);
end

if ~isempty(flist5)
    f5 = features_FDIM(im, mask, psize);
end

%concatenate all features:
ftot = cat(1, f1, f2, f3, f4, f5);
itot = cat(1, i1(:), i2(:), i3(:), i4(:), i5(:));

%Sort features in the same order as the input list
f(itot) = ftot(:);

%Return list of features that were successfully computed:
full_list(itot) = [flist1(:); flist2(:); flist3(:); flist4(:); flist5(:)];

% Check whether there were features that could not be computed:
sd = setdiff(1:length(fname), itot);
if ~isempty(sd)
    full_list(sd) = [];
    f(sd) = [];
    warning('the following features are unknown:')
    disp(fname(sd))
end
end

function [flist, i] = feature_list(ref_list, flist)
% Return the list of strings in 'flist' that are present
% in 'ref_list'.

remov = false(size(flist));
for n = 1:length(flist)
    remov(n) = ~any(strcmpi(ref_list, flist{n}));
end
i = find(~remov);
flist(remov) = [];

end