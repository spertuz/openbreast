function [c, d, l] = pclassify(dataset, model, K)
% Classify mamography image using KNN
% Sintax:
%     [c, f] = pclasssify(dataset, model, K)
% Inputs:
%     dataset,  Mx1 cell array with image paths
%     model,    classification model obtained from hmodel.
%     K,        number of neareast neighbors per class used
%               in the classification on the input image.
% Outputs:
%     c,      Mx1 vextor with predicted class for each input image.
%     d,      MxK feature vector computed for input image. The feature
%             vector corresponds to the distances to the K closest 
%             centroids in the model.%             
%     l,      MxK labels corresponding to each distance in d
%             
% S. Pertuz
% Sep30/2017

% Apply recursively for multiple images %%%%%%%%%%%%%%%%%%%%%%%%
if length(dataset)>1
    m = length(dataset);
    c = zeros(m, 1);
    d = zeros(m, 2);
    l = zeros(m, K);
    fprintf('Classifying      ')
    for n = 1:m        
        [c(n), d(n,:), l(n,:)] = pclassify(dataset(n), model, K);
        fprintf('\b\b\b\b\b[%02d%%]', floor(100*n/m))
    end
    fprintf('\n')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%obtain texton histogram from input image:
h = phistogram(dataset, model.P, model.params);
% Compute Xi-squared distance to each 
% histogram in the model:
% H = repmat(h, [size(model.H, 1), 1]);
Xi = 0.5*sum( ((model.H - h).^2+eps)./(model.H + h +eps), 2);

%find K nearest neighbors overall
[~, i] = sort(Xi, 'ascend');

% distances of k-closest centroids:
d(1) = max([min(Xi(~model.class(i(1:K)))), 0]);         
d(2) = max([min(Xi(model.class(i(1:K)))), 0]);

% class of k-closest centroids
l = model.class(i(1:K)); 

% Most frequent class in k-closest centroid
c = mode(l);