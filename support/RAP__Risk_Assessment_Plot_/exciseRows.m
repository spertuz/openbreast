function X = exciseRows(X)
X(any(isnan(X),2),:) = [];