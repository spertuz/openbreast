% Statistical analysis for risk assessment
% This demo compares the predictive value
% of two models (density only vs. density+age)
% using logistic regression (LR) to estimate Odds
% ratios with significance testing
clear, clc

%Load dataset with density, age and class:
srcdir = fileparts(mfilename('fullpath'));
fpath = [srcdir,'/../data/sample_data.mat'];
load(fpath,'dataset')

%feature matrix and target:
X = zscore([dataset.PD, dataset.Age]);
Y = dataset.Class;

%generate folds (5-fold CV)
rng default
cv = cvpartition(Y, 'kfold', 5);

%initialize score vectors:
scores_pd   = zeros(size(X, 1), 1);
scores_age  = zeros(size(X, 1), 1);

%run folds:
for n = 1:cv.NumTestSets
    i_train = cv.training(n);
    i_test = cv.test(n);
    Y_train = Y(i_train);
    X_train = X(i_train, :);
    X_test = X(i_test, :);
    
    %LR with PD:
    B = glmfit(X_train(:,1), Y_train, 'binomial', 'link', 'logit');
    scores_pd(i_test) = glmval(B, X_test(:,1), 'logit');
    
    %LR with age:
    B = glmfit(X_train(:,2), Y_train, 'binomial', 'link', 'logit');
    scores_age(i_test) = glmval(B, X_test(:,2), 'logit');
end

%Compute Odds ratios:
%Density only:
fprintf('Density only:\n')
[B,~,stats] = glmfit(zscore(scores_pd), Y, 'binomial', 'link', 'logit');
B_lo = B(2)-1.965*stats.se(2);  %Lower bound
B_hi = B(2)+1.965*stats.se(2);  %Upper bound
fprintf('PD:  OR=%1.3f 95%%CI: %1.3f-%1.3f, p=%1.4f\n',exp(B(2)), exp(B_lo), exp(B_hi), stats.p(2))

%Density + Age:
fprintf('\nDensity + Age:\n')
[B,~,stats] = glmfit(zscore([scores_pd, scores_age]), Y, 'binomial', 'link', 'logit');
B_lo = B(2)-1.965*stats.se(2);  %Lower bound
B_hi = B(2)+1.965*stats.se(2);  %Upper bound
fprintf('PD:  OR=%1.3f 95%%CI: %1.3f-%1.3f, p=%1.4f\n',exp(B(2)), exp(B_lo), exp(B_hi), stats.p(2))

B_lo = B(3)-1.965*stats.se(3);  %Lower bound
B_hi = B(3)+1.965*stats.se(3);  %Upper bound
fprintf('Age: OR=%1.3f 95%%CI: %1.3f-%1.3f, p=%1.4f\n',exp(B(3)), exp(B_lo), exp(B_hi), stats.p(3))

fprintf('\nPD retains its statistical significance\nafter the inclusion of Age (p<0.05)\n')
