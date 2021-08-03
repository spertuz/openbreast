function [imn, im] = ffdmRead(impath, info)
% Read FFDM image from DICOM file
% Sintax:
%     [imn, im] = ffdmRead(impath, info)
% Inputs:
%     impath,     path to DICOM file
%     info,       structure as returned by getinfo
% Outputs:
%     imn,        output intensity image with values
%                 normalized in [0,1]
%     im,         output with original intensity 
%                 values
%                 
% S. Pertuz
% Oct27/2017

% info = getinfo(impath);
im = double(dicomread(impath));

%find vendor:
vendor = upper(info.source);
is_negative     = contains(vendor, 'FUJIFILM')||...
    contains(vendor, 'SECTRA')||...
    contains(vendor, 'PHILIPS');
is_afga         = contains(vendor, 'AGFA');


%normalize:
if info.israw
    gmax = 2^14-1;
    gmin = 1;
    im(im<gmin) = gmin;    
    im(im>gmax) = gmax;
    imn = log(im);
    imn = (max(gmax) - imn).^2;
    imn = mat2gray(imn);
elseif is_negative
    gmax = max(im(:));
    gmin = 0;
    im = gmax-im;
    im(im<gmin) = gmin;
    im(im>gmax) = gmax;
    imn = mat2gray(im);
elseif is_afga
    gmax = max(im(:));
    gmin = 0;
    im = gmax-im;
    im(im<gmin) = gmin;
    im(im>gmax) = gmax;
    imn = mat2gray(im);
else
    gmax = 2^12-1;
    gmin = 0;
    imn = im;
    imn(im>gmax) = gmax;
    imn(im<gmin) = gmin;
    imn = mat2gray(imn);
end

% display
if nargout<1
    imshow(imn)
end