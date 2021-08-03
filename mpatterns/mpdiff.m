function DP = mpdiff(model, dataset)
% Discriminative power of MP patterns
% Sintax:
%     cdd = mpdiff(model)
%     cdd = mpdiff(model, dataset)
% Inputs:
%     dataset,    structure with test images
%     model,      structure with parameters of the model
% 
% Outputs:
%     D,          NxN difference matrix
%     
% S. Pertuz
% Feb 08/2018

fullflag = (nargin>1)&&~isempty(dataset);

if fullflag
    no_images = length(dataset);
    D = zeros(size(model.P, 2), size(model.P, 2), no_images);
    
    %For multiple images
    if no_images>1
        fprintf('Computing      ')
        for n = 1:no_images
            D(:,:,n) = mpdiff(model, dataset(n));
            fprintf('\b\b\b\b\b[%02d%%]', floor(100*n/no_images))
        end
        fprintf('\n')
        Y = [dataset.class];
        DP = zeros(size(D, 1), size(D, 2));
        for m = 1:size(D, 1)
            for n = m + 1:size(D, 2)
                X = squeeze(D(m,n,:));
                DP(m,n) = bestStump(X(:), Y(:));
                DP(n,m) = DP(m,n);
            end
        end
        return
    end
else
    DP = zeros(1, size(model.H, 2));
    Y = [model.class];
    for n = 1:size(model.H, 2)
        X = model.H(:,n);
        DP(n) = bestStump(X(:), Y(:));
    end
    return
end

% for only one image (full DP):
S = mpSampling(dataset, model.params);
no_samples = size(S, 2);
DP = zeros(size(model.P,2));
for m = 1:size(model.P,2)
    for n = m+1:size(model.P,2)
        h1 = sqrt(sum((S-model.P(:,m)).^2));
        h2 = sqrt(sum((S-model.P(:,n)).^2));
        [~, i] = min([h1; h2]);
        n1 = sum(i==1)/no_samples;
        n2 = sum(i==2)/no_samples;
        DP(m,n) = n1 - n2;
    end
end
end