function [f, info] = ffdmFeatures(impath, res)
% Extract full set of features using
% different scales, sampling strategies,
% and normalization methods
% Sintax:
%     [f, info] = ffdmFeatures(impath)
%     [f, info] = ffdmFeatures(impath, res)
% Inputs:
%     impath,     path to DICOM file
%     res,        resolution for analysis. Default is 0.07 mm/pixel
% Outputs:
%     f,          Nx1 vector with N features
%     info,       structure with the following fields:
%           .fnames,    Nx1 cell array with feature names
%           .eflag,     true if there was an error
%           .norm,      Nx1 array with normalization 1: none, 2:zscore, 3:SMF
%           .samp,      Nx1 array with sampling: 1: full, 2: RA, 
%                       3:Multi, 4:SQ
%           .scal,      Nx1 array with scale. 1: 1x, 2:0.5x
% 
% S. Pertuz
% Feb14/2018
    
 %global parameters:
if nargin<2   
    res = 0.07;  %image resolution (mm/pixel)
end

% list of features to compute (all):
flist = {'imin', 'imax', 'iavg', 'ient', 'istd', 'ip05', 'ip95', 'iba1', 'iba2', 'ip30', 'ip70', 'iske', 'ikur','iran',...
'cene', 'ccor', 'ccon', 'chom', 'cent',...
'rsre', 'rlre', 'rgln', 'rrpe', 'rrln', 'rlgr', 'rhgr',...
'sgra', 'slap', 'swas', 'swav', 'swar', 'stev','fdim'};

scales = {'x1', 'x3'};
normal = {'n1', 'n2', 'n3'};
sampl = {'s1','s2','s3','s4'};

% parameters for RA region:
slims       = [.10 .50];
tlims       = [.10 .90];

try
% retrieve image info
info = getinfo(impath);

%read image:
[imn, im] = ffdmRead(impath, info);

% if ~info.israw
%     normal(3) = []; 
% else
%     imr = imresize(im, info.psize/res);    
%     im = imn;    
% end
normal(3) = [];
if info.israw
    imr = imresize(im, info.psize/res);
    im = imn;
end

%breast mask
[mask, contour, cwall] = segBreast(imresize(mat2gray(imn),1/8), info.ismlo);

%re-scale if necessary
im      = imresize(im, info.psize/res);

mask0   = imresize(mask, size(im), 'nearest');
mask0   = imerode(mask0, ones(31));

%st-mapping:
rpoints     = refpoints(contour, cwall);
mapp        = stmap(contour, rpoints);
mask_RA0    = st2mask(mapp, slims, tlims);                    
mask_SQ0    = sqmax(mask0);

%initialize output
no_feats    = length(flist)*length(scales)*length(sampl)*length(normal);
f               = zeros(1, no_feats);
clear info
info.fnames     = cell(1, no_feats);
info.eflag      = false;
info.scal       = zeros(1, no_feats);
info.norm       = zeros(1, no_feats);
info.samp       = zeros(1, no_feats);

k = 1;
for nn = 1:length(normal)
    switch nn 
        case 1 %n1: no normalization
            imn = im;
        case 2 %n2: zscore
            imn = (im-mean(im(mask0)))/std(im(mask0));
        case 3 %n3: SMF
            imn = ffdmSMF(imr, info, mask0);            
    end
    imin = min(imn(:));
    imax = max(imn(:));
    for nx = 1:length(scales)
        switch nx
            case 1 %x1: full scale
                imx = imn;                
            case 2 %x2: 1/2 scale
                imx = imresize(imn, 1/2);
                imx(isnan(imx))= imin;
                imx(imx<imin) = imin;
                imx(imx>imax) = imax;
            case 3 %x3: 1/4 scale
                imx = imresize(imn, 1/4);                
                imx(isnan(imx)) = imin;
                imx(imx<imin) = imin;
                imx(imx>imax) = imax;
            case 4 %x4: 1/8 scale
                imx = imresize(imn, 1/8);                
                imx(isnan(imx)) = imin;
                imx(imx<imin) = imin;
                imx(imx>imax) = imax;
        end
        mask = imresize(mask0, size(imx), 'nearest');        
        mask_RA = imresize(mask_RA0, size(imx));
        mask_SQ = imresize(mask_SQ0, size(imx));
        
        for ns = 1:length(sampl)
            %feature indici:
            k0 = (k-1)*length(flist) + 1;
            k1 = k*length(flist);
            k = k + 1;
            
            %feature names:
            for nf = 1:length(flist)
                info.fnames{k0+nf-1} = sprintf('%s%s%s_%s',normal{nn},scales{nx},sampl{ns},flist{nf});
                
            end
            info.norm(k0:k1) = nn;
            info.scal(k0:k1) = nx;
            info.samp(k0:k1) = ns;
            
            switch ns
                case 1 %s1: full breast
                    f(k0:k1) = xfeatures(imx, flist, mask, res);
                case 2 %s2: RA region                                        
                    f(k0:k1) = xfeatures(imx, flist, mask_RA, res);                    
                case 3 %s3: multi-window 63x63                    
                    [x,y]       = rcenters([128 128], 63, mask);                    
                    [~, r]      = xySampling(imx, x, y, 63);
                    f(k0:k1)    = mean(xfeatures(r, flist, [], res), 1, 'omitnan');                    
                    
                case 4 %s4: squared ROI
                    f(k0:k1) = xfeatures(imx, flist, mask_SQ, res);                    
            end
        end
    end
end

catch
    info.eflag = true;
    f = zeros(1, 528);
end