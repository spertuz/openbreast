function [sq, rect] = sqmax(mask)
% Find largest circunscribed squared ROI
% Sintax:
%     sq = sqmax(mask)
% Inputs:
%     mask,   MxN input binary mask with ROI
% Outputs:
%     sq,     MxN binary mask with the largest
%             square circunscribed in the input ROI
%             
% This function is based on:
% [1] Jarek Tuszynski, Instribed Rectangle Package, Mathworks,
%  online: https://se.mathworks.com/matlabcentral/fileexchange/
%  28155-inscribed-rectangle. Last visited: Jan. 2017.
%  
% S. Pertuz
% Jan09/2018


%parameters:
sq = false(size(mask));
s = FindLargestSquares(mask);
[y, x] = find(s==max(s(:)), 1);
sq(y:y+s(y,x), x:x+s(y,x)) = true;
rect = [x, y, s(y,x)-1, s(y,x)-1];