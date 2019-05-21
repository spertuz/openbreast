function d = dsc(mask1, mask2)
% Dice's similarity coefficient

num = 2*sum(mask1(:)&mask2(:));
den = sum(mask1(:)) + sum(mask2(:)); 
d = num/den;