function f = features_GLRL(im, flist, mask, par)
% Image features from Gray Level Run Lengths  [1]
%
% SINTAX:
%       F = features_GLRL(im, flist)
%       F = features_GLRL(im, flist, mask)
%       F = features_GLRL(im, flist, mask, params)
% 
% flist,    A cell array of strings specifying the different
%           measures you want to get. The valid strings are the
%           following:
%                   'rSRE'   Short Runs Emphasis
%                   'rLRE'   Long Runs Emphasis
%                   'rGLN'   Gray Level Nonuniformity                   
%                   'rRPE'   Run Percentage
%                   'rRLN'   Run Length Nonuniformity
%                   'rLGR'   Low gray-level run emphasis
%                   'rHGR'   High gray-level run emphasis
%
%                   Features:
% {'rSRE', 'rLRE', 'rGLN', 'rRPE', 'rRLN', 'rLGR', 'rHGR'};
%
% REFERENCES:
% [1] Galloway, Mary M. "Texture analysis using gray level run lengths."
%     Computer graphics and image processing 4.2 (1975): 172-179.
% [2] https://se.mathworks.com/matlabcentral/fileexchange/52640-gray-level-run-length-image-statistics
%
% Copyright 2017, German F. Torres and Said Pertuz.

% Turning off the warning
warning('off','images:removing:function');

%Default parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    par.nlevels = 256;
end

if nargin<3
    mask = true(size(im));
end

f = zeros(length(flist), 1);
gmin = min(im(mask));
gmax = max(im(mask));
% gmin = 0;
% gmax = 255;
offset = 1:4;
GLRL  = grayrlmatrix(im, offset, par.nlevels, [gmin, gmax], mask);
[SRE, LRE, GLN, RLN, RP, LGRE, HGRE] = grayrlprops(GLRL, sum(mask(:)));

for n = 1:length(flist)
    switch upper(flist{n})
        case 'RSRE'
            f(n) = mean(SRE);
        case 'RLRE'
            f(n) = mean(LRE);
        case 'RGLN'
            f(n) = mean(GLN);
        case 'RRPE'
            f(n) = mean(RP);
        case 'RRLN'
            f(n) = mean(RLN);
        case 'RLGR'
            f(n) = mean(LGRE);
        case 'RHGR'
            f(n) = mean(HGRE);
    end
end