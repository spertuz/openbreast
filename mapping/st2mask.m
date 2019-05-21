function mask = st2mask(MAPP, slims, tlims)
% Convert st-map to binary mask
% Sintax:
%     mask = st2mask(MAPP, slims, tlims)
% 
% Inputs:
%     MAPP,   ST map as returned by the function stmap
%     slims,  a vector [smin, smax] with the s-limits
%             of the region of interest
%     tlims,  a vector [tmin, tmax] with the t-limits
%             of the region of interest
% 
% Outputs:
%     mask,   an MxN binary mask with the ROI delimited
%             by slims and tlims
%             
% S. Pertuz, F. Torres
% Jul04/2017

NS = MAPP.size(2);   %Number of points (S-grid)
NT = MAPP.size(1);   %Number of points (T-grid)
s = linspace(slims(1), slims(2), NS);
t = linspace(tlims(1), tlims(2), NT);
[S,T] = meshgrid(s, t);

x1 = round((T.*polyval(MAPP.P1, S) + polyval(MAPP.P2, S)).*...
    cos((MAPP.theta(2) - MAPP.theta(1))*T + MAPP.theta(1)) -...
    (T.*polyval(MAPP.P3, S) + polyval(MAPP.P4, S)).*...
    sin((MAPP.theta(2) - MAPP.theta(1))*T + MAPP.theta(1)) +...
    MAPP.x0);
y1 = round((T.*polyval(MAPP.P1, S) + polyval(MAPP.P2, S)).*...
    sin((MAPP.theta(2) - MAPP.theta(1))*T + MAPP.theta(1)) +...
    (T.*polyval(MAPP.P3, S) + polyval(MAPP.P4, S)).*...
    cos((MAPP.theta(2) - MAPP.theta(1))*T + MAPP.theta(1)) +...
    MAPP.y0);
mask = false(MAPP.size);
remov = (x1>MAPP.size(2))|(x1<1)|(y1>MAPP.size(1))|(y1<1);
x1(remov) = [];
y1(remov) = [];

index1 = unique(sub2ind(MAPP.size, y1, x1));

mask(index1) = true;
mask = imdilate(mask, [1 1;1 1]);

if MAPP.flip
    mask = fliplr(mask);
end
end