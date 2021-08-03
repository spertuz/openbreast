function imm = map2im(im, map0, map1, alpha)
% MP-map to image
% Sintax:
%     imm = map2im(im, map0, map1)
%     imm = map2im(im, map0, map1, alpha)
% 
% S. Pertuz
% Nov07/2018

if (nargin<4)
    alpha = 0.5;
end
im = imresize(im, size(map0));
imm = cat(3, im, im, im);
imm(:,:,1) = imm(:,:,1)+alpha*(map1);
imm(:,:,2) = imm(:,:,2)+alpha*(map0);