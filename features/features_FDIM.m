function f = features_FDIM(im, mask, psize, ftype)
% Fractal dimension:
% Sintax:
% 	f = features_FDIM(im)
% 	f = features_FDIM(im, mask, psize, ftype)
% Inputs:
% 	im,		MxN input image
% 	mask,	MxN optional binary mask (default is true(M,N))
% 	psize,	pixel pitch. Default is 0.1 mm/pixel
% 	ftype,	type of fractal dimension: 'boxc' for box counting
% 			or 'beta' for espectral power coefficient. Default
% 			is boxc.
% Outputs:
% 	f,		fractal dimension.
% 
% Reference:
% [1] Li, H., Giger, M. et al., Fractal Analysis of Mammographic Parenchymal
% Patterns in Breast Cancer Risk Assessment, Academic Radiology, 2007

if (nargin<2)||isempty(mask)
    mask = true(size(im));
end
if (nargin<3)||isempty(psize)
    psize = 0.1; %0.1mm (pixel size)
end
if (nargin<4)||isempty(ftype)
    ftype = 'boxc';
end

switch ftype
    case 'boxc'
        mask = imerode(mask, ones(2));
        Lmax = max(size(im));
        nmax = floor(log(Lmax)/log(2))-2;
        n = 0:1:min([7,nmax]);
        epsilon = psize*(2.^n(:));
        scale = psize*(1./epsilon);
        
        % Find surface area
        A = zeros(length(epsilon), 1);
        
        for n = 1:length(epsilon)
            im_e = imresize(im, scale(n));
            mask_e = imresize(mask, size(im_e));
            diff_x = [diff(im_e, 1, 2), zeros(size(im_e, 1), 1)];
            diff_y = [diff(im_e, 1, 1); zeros(1, size(im_e, 2))];
            Ax = abs(diff_x(mask_e)) + abs(diff_y(mask_e));
            A(n) = epsilon(n)^2 + epsilon(n)*sum(Ax(:));
        end
        
        %Find slope of log(A) vs log(e):
        p = polyfit(log(epsilon), log(A), 1);
        f = -p(1);
    case 'beta'
        nbands = 16;
        sc = mean(size(mask)/size(im));
        dx = round(sc*5/psize); %security ring of 1cm
        mask(:,1) = false;
        mask(1,:) = false;
        mask(end,:) = false;
        mask(:,end) = false;
        mask = imerode(mask, fspecial('disk', dx)~=0);        
        [~, rect] = sqmax(mask);
        rect = round(rect/sc);
        L = rect(3);
        imc = imcrop(im, rect);
        win = hanning(L+1);
        win = win*win';
        win = win/max(win(:));
        imc = imc.*win;
        imc = imc-mean2(imc);
        F = fftshift(fft2(imc));
        P = abs(F).^2;
        fs = 1/psize;
        f = linspace(-0.5*fs, 0.5*fs, L);
        [u, v] = meshgrid(f, f);
        f = sqrt(u.^2 + v.^2);
        flim = linspace(0, 0.5*fs, nbands + 1);
        Pn = zeros(nbands, 1);
        fn = zeros(nbands, 1);
        for n = 1:nbands
            select = (f<flim(n+1))&(f>flim(n));
            Pn(n) = mean(P(select));
            fn(n) = 0.5*(flim(n)+flim(n+1));
        end        
        p = polyfit(log(fn), log(Pn), 1);
        f = -p(1);
end