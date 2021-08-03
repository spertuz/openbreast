function Ps = mpSort(P, C)
% Sort micro-pattern prototypes
% Sintax:
%     Ps = mpSort(P, C)
% Sorts micro-pattern prototypes according to Fisher's
% linear discriminant score
% 
% Inputs:
%     P,      MxN array with N M-dimensional prototypes
%     C,      1xN boolean vector with class of each prototype
%     
% Outputs:
%     Ps,     MxN array with sorted prototypes. They are
%             sorted according to increasing Fisher's discriminat
%             score
%             
% S. Pertuz
% Jan29/2018

mdl = fitcdiscr(P', C(:), 'OptimizeHyperparameters', 'auto',...
    'scoreTransform', 'logit');
[~, sc] =  predict(mdl, P');
[~, i] = sort(sc(:,2), 'ascend');
Ps = P(:, i);