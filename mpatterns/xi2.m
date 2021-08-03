function G = xi2(U, V)
% Xi-squared kernel

% G = 1 - 0.5*sum( (U-V).^2./(U+V) );

m = size(U, 1);
n = size(V, 1);
G = zeros(m, n);

for i = 1:m
    u = U(i,:);
    G(i,:) = 1-0.5*sum( (u-V).^2./(u+V), 2)';
end
