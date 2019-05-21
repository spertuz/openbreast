function [Z, Pz, Beta, stats] = LR_Pz_choice( Y,Y_choice,Disease, varargin );
% *********************  LR_Pz_choice  ****************************
%   (c) John W Pickering, 2009
%     Christchurch Kidney Research Group
%     University of Otago Christchurch
%     New Zealand
%
%   Last update:  17 July 2012
%
% LR_Pz_choice by John Pickering is licensed under a 
%	Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in 
%     the documentation and/or other materials provided with the distribution
%
% Attribution to John Pickering. 
% *************************************************************************

%  PURPOSE
%   LR_Pz calculates the logistic regression for parameters in Y defined by
%   Y_choice and outputs the probablity vector P(z) and the beta
%   coefficients in the structure Pz_structure
%   Diseases, boolean vector defining outcome
%   varargin: boolean vectors defining the patients to include.

%  INPUTS
%  Y: Matrix of input variables - each column represents a new variable
%  Y_choice: array representing which columns to put in the model
%  Disease: boolean array where 1 represents those with the
%  disease/event/outcome of interest
%  varargin: boolean arrays for choosing the cohort of interest

%  OUTPUTS
%   Z: Z values for each subject
%   Pz: Probablity of Disease for each subject Pz=1/(1+ e^(-z))
%   Beta:  beta values of the model

%  Choose only the columns in Y defined by the vector Y_choice and only the patients defined by varargin
bcohort=choice(Disease,varargin{:});

X=Y(:,Y_choice);

excl=~choice(varargin{:});
X(any(excl,2),:)=[];
bcohort(any(excl,2),:)=[];

% logistic regression
[b,dev,stats]=glmfit(X,bcohort,'binomial');  %put here the model that best suits your data

n=size(b,1);
myY=zeros(size(X,1),1);
myY(:,1)=b(1);
for i=2:n
    myY=myY+b(i)*X(:,i-1);
end

%Caclulate the AUC for the logistic regression model
Z=myY;
Pz=1./(1+exp(-myY));
Beta=b;

end