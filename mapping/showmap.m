function showmap(im, MAPP)
% Show st-mapping
% Sintax:
%     showmap(im, MAPP)
% Inputs:
%     im,     MxN input image
%     MAPP,   structure with parameters of ST mapping
%     
% S. Pertuz
% Jan/2018


M = 30;             %No. of points in Y-grid
N = 30;             %No. of points in X-grid

imshow(mat2gray(im)), hold on

s = linspace(0, 1);
P1 = polyval(MAPP.P1, s);
P2 = polyval(MAPP.P2, s);
P3 = polyval(MAPP.P3, s);
P4 = polyval(MAPP.P4, s);

%Show t-grid:
t = linspace(0, 1, M);
theta = (MAPP.theta(2) - MAPP.theta(1))*t + MAPP.theta(1);
for m = 1:M
    PA = t(m)*P1 + P2;
    PB = t(m)*P3 + P4;
    x = PA*cos(theta(m)) - PB*sin(theta(m));   
    y = PA*sin(theta(m)) + PB*cos(theta(m));
    plot(x + MAPP.x0, y + MAPP.y0, 'g')
end

%Show s-grid:
s = linspace(0, 1, N);
t = linspace(0, 1);
theta = (MAPP.theta(2) - MAPP.theta(1))*t + MAPP.theta(1);
for n = 1:N
    PA = t*polyval(MAPP.P1, s(n)) + polyval(MAPP.P2, s(n));
    PB = t*polyval(MAPP.P3, s(n)) + polyval(MAPP.P4, s(n));    
    x = PA.*cos(theta) - PB.*sin(theta);   
    y = PA.*sin(theta) + PB.*cos(theta);
    plot(x + MAPP.x0, y + MAPP.y0, 'g')
end

% %Valid region (in red):
% %First boundary:
% s = linspace(0, 1);
% x1 = polyval(Px0, s);
% y1 = polyval(Py0, s);
% 
% %Second boundary:
% x2 = polyval(Px1, s);
% y2 = polyval(Py1, s);
% 
% %Third boundary:
% t = linspace(0, 1);
% theta = (theta1 - theta0)*t + theta0;
% PA = t*polyval(MAPP.P1, 1) + polyval(MAPP.P2, 1);
% PB = t*polyval(MAPP.P3, 1) + polyval(MAPP.P4, 1);
% x3 = PA.*cos(theta) - PB.*sin(theta);
% y3 = PA.*sin(theta) + PB.*cos(theta);


% Show reference points
plot([MAPP.x0, MAPP.x1, MAPP.x2], [MAPP.y0, MAPP.y1, MAPP.y2],...
'yo', 'linewidth',2)
