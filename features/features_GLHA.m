function f = features_GLHA(im, flist, mask)
% Gray-level histogram analysis 
% Sintax:
%     f = features_GLHA(im, flist)
%     f = features_GLHA(im, flist, mas)
%     
% Features:
% {'imin', 'imax', 'iavg',...
%     'ient', 'istd', 'ip05',...
%     'ip95', 'iba1', 'iba2',...
%     'ip30', 'ip70', 'iske',...
%     'ikur','iran'};
% 
% S. Pertuz
% Jul11/2017

if nargin<3
    mask = true(size(im));
end

f = zeros(length(flist), 1);
x = im(mask(:));

for n = 1:length(flist)
    switch lower(flist{n})
        case 'imin'
            f(n) = min(x);
        case 'imax'
            f(n) = max(x);
        case 'iavg'
            f(n) = mean(x);
        case 'ient'
            c = hist(x, 256);
            p = c/sum(c);
            f(n) = -sum(p(p~=0).*log2(p(p~=0)));
        case 'istd'
            f(n) = std(x);
        case 'ip05'
            f(n) = prctile(x, 5);
        case 'ip95'
            f(n) = prctile(x, 95);
        case 'ip30'
            f(n) = prctile(x, 30);
        case 'ip70'
            f(n) = prctile(x, 70);
        case 'iba1'
            p05 = double(prctile(x, 5));
            p95 = double(prctile(x, 95));
            u   = double(mean(x));
            f(n) = (p95 -u + eps)/(u-p05 + eps);
        case 'iba2'
            p30 = double(prctile(x, 30));
            p70 = double(prctile(x, 70));
            u   = double(mean(x));
            f(n) = (p70 -u + eps)/(u-p30 + eps);
        case 'iske'
            f(n) = skewness(double(x));
            %if isnan(f(n)), f(n)=0; end
        case 'ikur'
            f(n) = kurtosis(double(x));
            %if isnan(f(n)), f(n)=0; end
        case 'iran'
            f(n) = range(double(x));
        otherwise
            error('unknown feature %s', upper(flist{n}))
    end
end