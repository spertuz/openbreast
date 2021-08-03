function best_acc = bestStump(x, y, dispflag)
if nargin<3
    dispflag = false;
end
xth = linspace(min(x), max(x));
acc = zeros(size(xth));
for n = 1:length(xth)
    acc(n) = sum(double(x>xth(n))==double(y))/numel(x);
end

best_acc = max(acc);

if dispflag
    figure, subplot(121), plot(xth, 100*acc)
    xlabel('Frequency threshold'), ylabel('Accuracy')
    xth = linspace(min(x), max(x), 15);
    subplot(122), bar(xth, [hist(x(y), xth); hist(x(~y), xth)]')
    xlabel('Frequency threshold'), ylabel('Number of MP samples')
    legend({'Case','Controls'})
end