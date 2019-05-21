function showcorr(X, fnames)
% Visualize correlation among features
% Sintax:
%     showcorr(X, fnames)
% Inputs:
%     X,      MxN matrix of N-dimensional feature vectors
%     fnames, 1xN cell array with feature names
% 
% S. Pertuz
% Feb08/2018

L = 1024;
[R,p] = corrcoef(X);
[~, i] = sort(R(1,:), 'descend');
R = corrcoef(X(:,i));
xlims = linspace(1, L, size(R, 1)+1);
xi = 0.5*(xlims(1:end-1)+xlims(2:end));
imagesc(imresize(R, [L, L], 'nearest'));
set(gca, 'ytick', xi, 'yticklabel', fnames(i))
set(gca, 'xtick', xi, 'xticklabel', fnames(i), 'xticklabelrotation', 45)
axis on
daspect([1 1 1])
colorbar
 set(gcf,'units','centimeters',...
            'position',[2 2 18 16],...
            'paperunits','centimeters',...
            'papersize',[18 16])
 set(gca, 'fontname', 'monospaced')
s = zeros(size(p));
s(p<0.05) = 128;
s(p<0.01) = 192;
figure, imagesc(imresize(s, [L, L], 'nearest'), [0 255])
colormap gray
set(gca, 'ytick', xi, 'yticklabel', fnames(i))
set(gca, 'xtick', xi, 'xticklabel', fnames(i), 'xticklabelrotation', 45)
axis on
daspect([1 1 1])
cb = colorbar;
 set(gcf,'units','centimeters',...
            'position',[2 2 18 16],...
            'paperunits','centimeters',...
            'papersize',[18 16])
 set(gca, 'fontname', 'monospaced')
 pos = get(gca, 'position');
 delete(cb)
 set(gca, 'position', pos)
 
 n_weak = (sum(p(:)<0.05&abs(R(:))<0.3)-length(fnames))/2;
 n_moderate = (sum(p(:)<0.05&(abs(R(:))>=0.3 & abs(R(:))<0.7))-length(fnames))/2;
 n_large = (sum(p(:)<0.05&(abs(R(:))>=0.7))-length(fnames))/2;
 fprintf('Weak correlation: %d\n', round(n_weak))
 fprintf('Moderate correlation: %d\n', round(n_moderate))
 fprintf('Large correlation: %d\n', round(n_large))

