function f = features_GLSM(im, flist, mask)
% Gray-level sharpness measure
% Sintax:
%     f = features_GLSM(im, flist)
%     f = features_GLSM(im, flist, mask)
%     
% Features:
% {'sgra', 'slap', 'swas',...
%     'swav', 'swar', 'stev'};
% 
% S. Pertuz
% Jul11/2017

if nargin<3
    mask = true(size(im));
end

f = zeros(length(flist), 1);

for n = 1:length(flist)
    switch lower(flist{n})
        case 'sgra'
            fx = imgradient(im);            
            f(n) = mean(fx(mask));
        case 'slap'
            h = [-1 2 -1];
            ly = imfilter(im, h, 'replicate', 'conv');
            lx = imfilter(im, h', 'replicate', 'conv');
            fx = abs(lx) + abs(ly);            
            f(n) = mean(fx(mask));
        case 'swas'
            [C, S] = wavedec2(im, 1, 'db6');
            H = wrcoef2('h', C, S, 'db6', 1);   
            V = wrcoef2('v', C, S, 'db6', 1);   
            D = wrcoef2('d', C, S, 'db6', 1);   
            fx = abs(H) + abs(V) + abs(D);            
            f(n) = mean(fx(mask));
        case 'swav'
            [C,S] = wavedec2(im, 1, 'db6');
            H = abs(wrcoef2('h', C, S, 'db6', 1));
            V = abs(wrcoef2('v', C, S, 'db6', 1));
            D = abs(wrcoef2('d', C, S, 'db6', 1));            
            f(n) = std(H(mask))^2 + std(V(mask))^2 +...
                std(D(mask))^2;            
        case 'swar'
            [C,S] = wavedec2(im, 3, 'db6');
            H = abs(wrcoef2('h', C, S, 'db6', 1));
            V = abs(wrcoef2('v', C, S, 'db6', 1));
            D = abs(wrcoef2('d', C, S, 'db6', 1));
            A1 = abs(wrcoef2('a', C, S, 'db6', 1));
            A2 = abs(wrcoef2('a', C, S, 'db6', 2));
            A3 = abs(wrcoef2('a', C, S, 'db6', 3));
            A = A1 + A2 + A3;
            WH = H.^2 + V.^2 + D.^2;
            WH = mean(WH(mask));
            WL = mean(A(mask));
            f(n) = (WH + eps)/(WL + eps);        
        case 'stev'
            g = imgradient(im);
            f(n) = std(g(mask));
	case 'suni' %uniformity
		fx = imgradient(im);
		h = hist(fx(mask), 256);
		h = h/sum(h);
		f(n) = sum(h.^2);
	case 'ssmo' %smoothness
		fx = imgradient(im);
		f(n) = 1/(1+ std(fx(mask))^2);
	case 'sske' %skewness
		fx = imgradient(im);
		f(n) = skewness(fx(mask));
	case 'sent' %entropy
		fx = imgradient(im);
		c = hist(fx(mask), 256);
        	p = c/sum(c);
            	f(n) = -sum(p(p~=0).*log2(p(p~=0)));
	case 'sfdi' %fractal dimension
		psize = 0.1; %pixel size in mm
		epsilon = psize*[1 2 4 8 16 32 64 128]';
		scale = psize*(1./epsilon);
		% Find surface area
		A = zeros(length(epsilon), 1);
		for n = 1:length(epsilon)
			im_e = imresize(im, scale(n));
			mask_e = imresize(mask, scale(n));
			diff_x = imfilter(im_e, [1, -1]);
			diff_y = imfilter(im_e, [1; -1]);
			Ax = abs(diff_x(mask_e(:))) + abs(diff_y(mask_e(:)));
			A(n) = epsilon(n)^2 + epsilon(n)*sum(Ax(:));
		end
		%Find slope of log(A) vs log(e):
		p = polyfit(log(epsilon), log(A), 1);
		f(n) = -p(1);
        otherwise            
            error('unknown method %s', flist{n})
    end
    
end
