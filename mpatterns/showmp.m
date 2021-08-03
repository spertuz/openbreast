function im = showmp(P, C, psize)
% Visualize MP prototypes
% Sintax:
%     im = showmp(P, C)
% Inputs:
%     P,      MxN matrix with N prototypes
%     C,      1xN class labels
% Output:
%     im,     output image
%     
% S. Pertuz
% Feb06/2018

L = 700;
edge = 10;
side_y = floor((L/2-(psize(2)+1)*edge)/psize(2));
side_x = floor((L-(psize(1)+1)*edge)/psize(1));
side = min([side_x, side_y]);
LX = side*psize(1) + edge*(psize(1)+1);
LY = side*psize(2) + edge*(psize(2)+1);

im_up = 0.9*ones(LY, LX);
P0 = P(:,~C);
P1 = P(:,C);
k = 1;
npix = sqrt(size(P0, 1));
yf = edge;
for m = 1:psize(2)
    yo = yf + edge;
    yf = yo + side - 1;
    xf = edge;
    for n = 1:psize(1)
        xo = xf + edge;
        xf = xo + side - 1;        
        im_up(yo:yf, xo:xf) = mat2gray(imresize(reshape(P0(:,k), [npix npix]), [side, side]));
        k = k + 1;
    end        
end

k = 1;
im_do = 0.9*ones(LY, LX);
yf = -edge + 1;
for m = 1:psize(2)
    yo = yf + edge;
    yf = yo + side - 1;
    xf = edge;
    for n = 1:psize(1)
        xo = xf + edge;
        xf = xo + side - 1;
        im_do(yo:yf, xo:xf) = mat2gray(imresize(reshape(P1(:,k), [npix npix]), [side, side]));
        k = k + 1;
    end    
end
im = [im_up; im_do];
imshow(im)