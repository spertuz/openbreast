function [Z, meanX, stdX] = NormalizeFeatures(X)
% Standarize an input feature matrix
%
% SINTAX:
%       Z = NormalizeFeatures(X)
%       [Z, meanX, stdX] = NormalizeFeatures(X)
%
% DESCRIPTION:
% Z = NormalizeFeatures(X) returns an output feature matrix whose mean and
% standard deviation are 0 and 1, respectively, for every feature in X. The
% additional outputs meanX and stdX are arrays that contain the mean and
% the standard deviation of every feature in X.
%
% INPUTS:
% X:        is a numeric matrix that contains the values of features. The
%           features are presented as columns in X.
%
% Copyright 2017, German F. Torres and Said Pertuz.

% Check X
validateattributes(X,{'numeric'},{'2d'},mfilename,'X',1);

% Number of observations
N = size(X,1);
% Number of variables
M = size(X,2);

% Output array of normalized values
Z = zeros(N,M);

% Subtract mean of each Column from data
mu = repmat(mean(X),N,1);
Z = X-mu;

% Normalize each observation by the standard deviation of that variable
sigma = repmat(std(X,0,1),N,1);
Z = Z./sigma;

% Additional outputs
if nargout>1
    meanX = mean(X);
    stdX = std(X,0,1);
end
end