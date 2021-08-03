function r = mpRank(model)
% Rank MP patterns
% Sintax: 
%     r = mpRank(model)
% Inputs:
%     model,  MP model as returned by pmodel
% Outputs:
%     r,      ranking of MP prototypes in model
%     
% S. Pertuz
% Feb08/2018

r = ILFS(model.H, double(model.class), 6);