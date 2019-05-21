function [mask1, mask2] = seg2mask(contour, cwall, rpoints, cgap, sgap)
% Convert segmentation to mask
% Sintax:
%     mask1 = seg2mask(contour)
%     mask1 = seg2mask(contour, cwall)
%     [mask1, mask2] = seg2mask(contour, cwall, rpoints)
% Inputs:
%     contour,    a structure with the breast contour (as returned
%                 by segBreast)
%     cwall,      100x2 array with x-y coordinates of points on
%                 the chest wall (as returned by cwall2pts)
%     rpoints,    structure with reference points (as returned
%                 by refpoints)
% Outputs:
%     mask1,      MxN binary mask wtih full breast
%     mask2,      MxN binary mask with RA region
%     
% S. Pertuz
% Jan09/2018

if nargin<3, rpoints = [];
end
if nargin<2, cwall = [];
end

% parameters for RA region:
slims       = [.2 .8];
tlims       = [.1 .9];

%Breast contour:
M = contour.size(1);
N = contour.size(2);
x = contour.x;
y = contour.y;

x = [0; x(1); x(:); x(end); 0];
y = [min(y); min(y); y(:);  y(end); y(end)];
mask1 = poly2mask(x, y, M, N);

%apply skin gap if necessary
if (nargin>4)&&~isempty(sgap)&&(sgap~=0)
    mask1 = imerode(mask1, ones(sgap));
end
        
%Chest wall:
if (nargin>1)&&~isempty(cwall)
    x = cwall(:,1);
    y = cwall(:,2);
    x = [0; x(:)];
    y = [0; y(:)];
    mask2 = poly2mask(x, y, M, N);
    mask2(1:end,1:2) = true;
    
    %apply chest gap if necessary
    if (nargin>3)&&~isempty(cgap)&&(cgap~=0)
        mask2 = imdilate(mask2, ones(cgap));
    end
    
    mask1 = mask1&~mask2;
end

%RA mask
if (nargin>2)&&(nargout==2)&&~isempty(rpoints)
    %Adjust reference points to contour:
    x = contour.x;
    y = contour.y;
    [~, i] = min(sqrt( (rpoints.p0.x-x).^2 + (rpoints.p0.y - y).^2));
    rpoints.p0.x = x(i);
    rpoints.p0.y = y(i);
    [~, i] = min(sqrt( (rpoints.p1.x-x).^2 + (rpoints.p1.y - y).^2));
    rpoints.p1.x = x(i);
    rpoints.p1.y = y(i);
    [~, i] = min(sqrt( (rpoints.p2.x-x).^2 + (rpoints.p2.y - y).^2));
    rpoints.p2.x = x(i);
    rpoints.p2.y = y(i);
    
    mapp = stmap(contour, rpoints);
    mask2   = st2mask(mapp, slims, tlims);
    
    %remove borders:
    mask2(:,1) = false;
    mask2(:,end) = false;
    mask2(1,:) = false;
    mask2(end,:) = false;
end

%Remove borders
mask1(:,1) = false;
mask1(:,end) = false;
mask1(1,:) = false;
mask1(end,:) = false;

if contour.flip
    mask1 = fliplr(mask1);    
end