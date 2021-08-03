function [auc, or, auc_CI, or_CI] = getperf(labels, scores, dispflag)
% Compute binary classification performance
% Sintax:
%     [AUC, OR] = getperf(labels, scores, dispflag)
%     [AUC, OR, aucCI, orCI] = getperf(labels, scores)
% Inputs:
%     labels, boolean vector with correct labels (true/false)
%     scores, scalar vector with prediction scores
% Outputs:
%     AUC,    area under the ROC curve
%     OR,     odds per standard deviation (OPERA)
%     aucCI,  AUC's 95% confidence interval [low, hi, pvalue]
%     orCI,   OR's 95% confidence interval [low, hi, pvalue]
%     
% S. Pertuz
% Oct24/2019

x = zscore(scores(:));
y = labels(:);

%Compute OR's with CI's
[B, ~, stats] = glmfit(x, y, 'binomial', 'link', 'logit');
B_lo = B(2)-1.965*stats.se(2);  %Lower bound
B_hi = B(2)+1.965*stats.se(2);  %Upper bound
or = exp(B(2));
or_CI(1) = exp(B_lo); %CI lower limit
or_CI(2) = exp(B_hi); %CI higher limit   
or_CI(3) = stats.p(2); %pvalue

%Compute AUC with CI's
res = AUC_compare_correlated([x, x], y);
auc = res{1}(1);
auc_CI = res{1}(2:3);

if dispflag||(nargout==0)
    out = table([or; auc], [or_CI(1,1:2); auc_CI(1,1:2)], [or_CI(3); nan], ...
        'VariableNames',{'value', 'CI', 'pvalue'},'RowNames', {'OR', 'AUC'});
    disp(out);
end