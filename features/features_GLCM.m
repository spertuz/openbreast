function f = features_GLCM(im, flist, mask, par)
% Texture features from GLCM
% Sintax:
%     f = features_GLCM(im, fname)
%     f = features_GLCM(im, fname, mask)
% 
% Features:
% {'cene', 'ccor', 'ccon',...
%     'chom', 'cent'};
% 
% S. Pertuz, G. F. Torres
% Jul11/2017

% Parameters for GLCM %%%%%%%%
if nargin<4
    par.length = 10; % by Hai; See the paper: Parenchymal texture analysis in digital mammography: A fully automated pipeline for breast cancer risk assessment
    par.nlevels = 128;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    mask = true(size(im));
end

f = zeros(length(flist), 1);
im(~mask) = NaN;
offset = par.length*[0 1;-1 1;-1 0;-1 -1];

warning('off','Images:graycomatrix:scaledImageContainsNan');
g = graycomatrix(im, 'Offset', offset, 'NumLevels', par.nlevels,...
    'Symmetric', false, 'GrayLimits', [min(im(mask)), max(im(mask))]); %by Hai, change 'Symmetric'=true to false
warning('on','Images:graycomatrix:scaledImageContainsNan');

for n = 1:length(flist)
    switch flist{n}
        case 'cene' %Energy (Angular second moment)
            p = graycoprops(g, 'energy');
            f(n) = mean(p.Energy);
        case 'ccor' %Correlation
            p = graycoprops(g, 'correlation');
            f(n) = mean(p.Correlation);
            %if isnan(f(n)), f(n)=1; end
        case 'ccon' %Contrast (Inertial)
            p = graycoprops(g, 'contrast');
            f(n) = mean(p.Contrast);
        case 'chom' %Homogeneity (Inverse difference moment)
            p = graycoprops(g, 'homogeneity');
            f(n) = mean(p.Homogeneity);
        case 'cent' %Entropy
            G = sum(sum(g));
            gsize = [par.nlevels, par.nlevels];
            G = repmat(G, gsize);       % Normalization Factor
            p = g./G;                   % Normalized GLCM
            p(p==0) = 1;
            E = -sum(sum( p.*log2(p) ));
            f(n) = mean(E);            
    end
end
