function model = pmodel(dataset, P, params)
% MP-based statistical population model from image dataset
% Sintax:
%     model = pmodel(dataset, P, params)
% Inputs:
%     dataset,    Mx1 structure with image data:
%                 .path,   string with image path
%                 .class,  integer with class label
%     params,   a structure with the following fields:  
%               .nsamples,  integer with the number of micro-pattern samples 
%                           per image.
%               .mpsize,    integer with micro-pattern size in pixels (odd 
%                           number).
%               .rflag,     boolean flag to apply feature extraction in the
%                           RA region. Default is false.
%               .nflag,     boolean flag to perform feature normalization.
%                           Default is false.
%     P,          NxNP array with N-dimensional MP prototypes
% 
% Outputs:
%     model,      a structure with the following fields:
%                 .H,       MxNP array with the response histogram
%                           of each image in dataset
%                 .C,       the same as the dataset.class
%                 .params,  the same as the input
%                 .P,       the same as the input
%                 
% S. Pertuz
% Sep30/2017

no_images = size(dataset, 1);
no_prototypes = size(P, 2);

%initialize output:
model.H         = zeros(no_images, no_prototypes);
model.class     = dataset.class(:);
model.params    = params;
model.P         = P;
model.C         = length(unique(dataset.class(:)));

fprintf('Building model      ')
for n = 1:no_images
        
    %compute histogram:
    model.H(n,:) = phistogram(dataset(n,:), P, params);
    fprintf('\b\b\b\b\b[%02d%%]', floor(100*n/no_images))
end
fprintf('\n')