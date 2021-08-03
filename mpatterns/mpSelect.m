function [P, CP] = mpSelect(M, C, np, pmax)
% Select micro-pattern prototypes
% Sintax:
%     [P, CP] = mpSelect(M, C)
%     [P, CP] = mpSelect(M, C, NP, pmin)
% Inputs:
%     M,      MxK array with micro-pattern samples. K is the 
%             total number of samples.
%     C,      1xK vector of class label of each micro-pattern
%     np,     is the sought number of micro-pattern prototypes
%     pmax,   maximum percentage of samples after anomaly-based pruning
% 
% Outputs:
%     P,      Nx2NP array with N-dimensional prototypes, where NP is the 
%             number of prototypes per class.
%             
% S. Pertuz
% Sep30/2017

% set defaults
if (nargin<4)||isempty(pmax)
    pmax = 0.2;
end

if (nargin<3)||isempty(np)
    np = 40;
end

P1 = M(:,C==1);     %class 1
P0 = M(:,C==0);     %class 0

if pmax<1   %prune M using anomaly detection
    nmax = round(pmax*size(M, 2)/2);
    mu = mean(P1, 2)';
    sd = cov(P1');
    p = mvnpdf(P1', mu, sd);
    [~, i] = sort(p, 'descend');
    P1 = P1(:,i(1:nmax));
    
    mu = mean(P0, 2)';
    sd = cov(P0');
    p = mvnpdf(P0', mu, sd);
    [~, i] = sort(p, 'descend');
    P0 = P0(:,i(1:nmax));
end
   
% select final prototypes using k-means:
warning off
[~, c_0] = kmeans(P0', np, 'MaxIter', 300);
[~, c_1] = kmeans(P1', np, 'MaxIter', 300);
warning on

% assign output
P = [c_0', c_1'];
CP = [false(1, size(c_0, 1)), true(1, size(c_1, 1))];
