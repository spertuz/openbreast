function [mask, contour] = ffdmForeground(im, cflag)
% Detect breast foreground in FFDM image
% Sintax:
%     mask = ffdmForeground(im, cflag)
% 
% Inputs:
%     im,     MxN grayscale mammography image
%     cflag,  binary flag. True activates contour cutting
%             using curvature. Default is false.
% Outputs
%     mask,       MxN binary mask with breast
%                 foreground.
%     contour,    structure with contour data
% 
% S. Pertuz
% Oct31/2017

%Parameters: %%%%
nbins   = 1000; % Number of bins for image histogram
%%%%%%%%%%%%%%%%%

if (nargin<2)||isempty(cflag), cflag = false;
end

% find intensity threshold:
warning off
xth = getLower(im(:), nbins);
warning on


% find mask:
mask0 = (im>=max([xth, 0]));

% remove artifacts and holes in the mask:
mask = cleanMask(mask0);

%compute contour:
contour = getcontour(mask, cflag);
contour.th = xth;
mask(round(max(contour.y)):end,:) = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xth = getLower(x, nbins)

%Minimum and maximum intesities:
x_min = min(x);
x_max = max(x);

%Relative frequency:
xi = linspace(x_min, x_max, nbins)';
n = histc(x, xi);

%Smooth histogram:
n = conv(n, gausswin(25), 'same');
n = n/max(n);

% find threshold by fitting a gaussian to the
% histogram peak(s) located below the mean.

xsup = min(mean(x), max(prctile(x, 30), .2));
[~,ipeaks] = findpeaks(n.*(xi<=xsup), 'MinPeakHeight', .35);

if (numel(ipeaks)==1)   %only one peak
    select = (n>.35)&(xi<xsup);
    f = fit(xi(select), n(select), 'gauss1');
    xth = f.b1 + sqrt(f.c1^2*(log(f.a1)-log(0.05)));        
    
elseif (numel(ipeaks)>1) %two peaks
    %find minimum between peaks
    n_min = min(n(min(ipeaks):max(ipeaks)));
    i_min = find(n==n_min, 1);
    
    %adjust second peak
    select = (xi>=xi(i_min))&(n>.35)&(xi<xsup);
    f = fit(xi(select), n(select), 'gauss1');
    xth = f.b1 + sqrt(f.c1^2*(log(f.a1)-log(0.05)));
    
elseif isempty(ipeaks)  %no peaks
    n_max = max(n);
    i_th = find(n<0.05*n_max, 1);
    xth = xi(i_th);
end

end

function mask = cleanMask(mask0)

% remove first and last row
mask0(1,:) = false;
mask0(end,:) = false;
mask0 = imerode(mask0, ones(5));

% keep biggest region
cc = bwconncomp(mask0);
stats = regionprops(cc, 'area');
idx = find([stats.Area]==max([stats.Area]), 1);
mask0 = ismember(labelmatrix(cc), idx);

% remove spurious holes:
mask0 = imdilate(mask0, ones(5));
mask0 = imclose(mask0, ones(5));
mask = imfill(mask0, 'holes');
end