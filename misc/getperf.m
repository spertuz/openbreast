function [auc, acc, fsc] = getperf(labels, scores, Y, dispflag)

if nargin<4
    dispflag = true;
end

acc = 100*sum(labels(:)==Y(:))/numel(labels);
[~,~,~,auc] = perfcurve(Y(:), scores, true);


p = sum(Y(:));
tp = sum(Y(:)&labels(:));
fp = sum(~Y(:)&labels(:));
prec = tp/(tp+fp);
rec = tp/p;
fsc = 2*prec*rec/(prec+rec);

if dispflag
    fprintf('Acc : %.2f\n', acc)
    fprintf('AUC : %1.3f\n', auc)
    fprintf('F-1 : %1.3f\n', fsc)
end

