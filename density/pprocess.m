function seg = pprocess(seg, mask, skin_gap, area_th, pixel_size)
% Post process density segmentation

% Convert parameters to pixels
skin_gap = round(2*skin_gap/pixel_size);
area_threshold = round(area_th/(pixel_size^2));

% Erode breast mask:
maske = imerode(mask, ones(skin_gap));

% Apply area filter:
seg = bwareaopen(seg&maske, area_threshold);

