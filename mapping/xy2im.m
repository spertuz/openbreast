function imout = xy2im(x, y, f, im, mask)

sc  = size(mask)./size(im);
im  = imresize(im, size(mask));
y   = floor(sc(1)*(y-1) + 1);
x   = floor(sc(2)*(x-1) + 1);

[xi, yi] = meshgrid(1:size(mask,2), 1:size(mask,1));
i = sub2ind(size(mask), y, x);
mask(i) = false;
im(i) = f;
imout = griddata(xi(~mask(:)), yi(~mask(:)), im(~mask(:)), xi, yi);
