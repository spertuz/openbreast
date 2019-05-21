function flag = isright(im)
% Determine whether mammogram is left or right.
% TRUE if the mammogram is in the right side of the image.

im(im>.95*max(im(:))) = 0;
s = sum(im);
n = round(0.5*size(im, 2));

flag = sum(s(1:n))<sum(s(n+1:end));
    