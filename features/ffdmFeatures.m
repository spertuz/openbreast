function [f, info] = ffdmFeatures(impath, res, scales, normal, sampl)
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
%           .norm,      Nx1 cell array with normalization 'none' or 'zscore'
%           .samp,      Nx1 cell array with sampling: 'full', 'RA', 
%                       'multi' or 'SQ'
%           .scal,      Nx1 array with scale. '1.0', '0.5' or '0.25' 
% 
% S. Pertuz
% Feb14/2018
    
 %global parameters:
if (nargin<2)||isempty(res)   
    res = 0.07;  %image resolution (mm/pixel)
end

if (nargin<3)||isempty(scales)
    scales = {'1.0', '0.5','0.25'};
end

if (nargin<4)||isempty(normal)
    normal = {'none', 'zscore'};
end

if (nargin<5)||isempty(sampl)
    sampl = {'full','RA','multi','SQ'};
end

% list of features to compute (all):
flist = {'imin', 'imax', 'iavg', 'ient', 'istd', 'ip05', 'ip95', 'iba1', 'iba2', 'ip30', 'ip70', 'iske', 'ikur','iran',...
'cene', 'ccor', 'ccon', 'chom', 'cent',...
'rsre', 'rlre', 'rgln', 'rrpe', 'rrln', 'rlgr', 'rhgr',...
'sgra', 'slap', 'swas', 'swav', 'swar', 'stev','fdim'};


% parameters for RA region:
slims       = [.10 .50];
tlims       = [.10 .90];

% try
% retrieve image info
image_info = getinfo(impath);

%read image:
[imn, im] = ffdmRead(impath, image_info);

if image_info.israw
    imr = imresize(im, image_info.psize/res);
    im = imn;
elseif length(normal)>2
    normal(3) = [];
end

%number of features:
no_feats    = length(flist)*length(scales)*length(sampl)*length(normal);

%breast mask
[mask, contour, cwall] = segBreast(imresize(mat2gray(imn),0.25), image_info.ismlo);

%re-scale if necessary
im      = imresize(im, image_info.psize/res);

mask0   = imresize(mask, size(im), 'nearest');
mask0   = imerode(mask0, ones(31));

%st-mapping:
rpoints     = refpoints(contour, cwall);
mapp        = stmap(contour, rpoints);
mask_RA0    = st2mask(mapp, slims, tlims);                    
mask_SQ0    = sqmax(mask0);

%initialize output
f               = zeros(1, no_feats);
info.eflag      = [];
info.fnames     = cell(1, no_feats);
info.scal       = cell(1, no_feats);
info.norm       = cell(1, no_feats);
info.samp       = cell(1, no_feats);

%flip side if necessary
if isright(imn)
    im = fliplr(im);    
    mask0 = fliplr(mask0);
    mask_RA0 = fliplr(mask_RA0);
    mask_SQ0 = fliplr(mask_SQ0);    
end

k = 1;
for nn = 1:length(normal)
    switch lower(normal{nn})
        case 'none' % no normalization
            imn = im;
        case 'zscore' % zscore
            imn = (im-mean(im(mask0)))/std(im(mask0));
%         case 'smf' %n3: SMF
%             imn = ffdmSMF(imr, info, mask0);            
    end
    imin = min(imn(:));
    imax = max(imn(:));
    for nx = 1:length(scales)
        switch scales{nx}
            case '1.0' %x1: full scale
                imx = imn;                
            case '0.5' %x2: 1/2 scale
                imx = imresize(imn, 1/2);
                imx(isnan(imx))= imin;
                imx(imx<imin) = imin;
                imx(imx>imax) = imax;
            case '0.25' %x3: 1/4 scale
                imx = imresize(imn, 1/4);                
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
            info.fnames(k0:k1) = flist(:);            
            info.norm(k0:k1) = normal(nn);
            info.scal(k0:k1) = scales(nx);
            info.samp(k0:k1) = sampl(ns);
            
            switch lower(sampl{ns})
                case 'full' %s1: full breast
                    f(k0:k1) = xfeatures(imx, flist, mask, res);
                case 'ra' %s2: RA region                                        
                    f(k0:k1) = xfeatures(imx, flist, mask_RA, res);                    
                case 'multi' %s3: multi-window 63x63                    
                    [x,y]       = rcenters([128 128], 63, mask);                    
                    [~, r]      = xySampling(imx, x, y, 63);
                    f(k0:k1)    = mean(xfeatures(r, flist, [], res), 1, 'omitnan');                                        
                case 'sq' %s4: squared ROI
                    f(k0:k1) = xfeatures(imx, flist, mask_SQ, res);                    
            end
        end
    end
end
info.eflag = 'success';
% catch
%     info.eflag = lasterr;
%     f = zeros(1, no_feats);
% end