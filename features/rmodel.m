function model = rmodel(x_train, y_train)
% Obtain risk model from training set
% Sintax:
%     model = rmodel(x_train, y_train)
%     model = rmodel(x_train, y_train)
% Train risk model using stepwise feature selection and
% logistic regression.
% 
% Inputs:
%     x_train,    NxM feature matrix. Each row is a 1xM feature vector
%                 corresponding to one image.
%     y_train,    class labels.
%     
% Outputs:
%     model,      a structure with output risk model
%     
% S. Pertuz
% Feb 05/2018


% standardize features:
model.mu = mean(x_train);
model.std = std(x_train);
x_train = (x_train-model.mu)./model.std;

% feature selection:
[~,~,~, model.select] = stepwisefit(x_train, y_train, 'display', 'off');

% logistic regression (fit):
warning off
model.weights = glmfit(x_train, y_train, 'binomial', 'link', 'logit');
warning on

% logistic regression (eval):
model.scores = glmval(model.weights, x_train, 'logit');
model.class = y_train;    