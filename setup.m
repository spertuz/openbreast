function setup()
% Add OpenBreast toolbox to the path
% 
% Author: Said Pertuz
% Copyright (c) 2018 Said Pertuz
% All rights reserved
% 
% Jan08/2018

[root, ~, ~] = fileparts(mfilename('fullpath'));

addpath(fullfile(root, 'demos'))
addpath(fullfile(root, 'segmentation'))
addpath(fullfile(root, 'misc'))
addpath(fullfile(root, 'mapping'))
addpath(fullfile(root, 'features'))
addpath(fullfile(root, 'mpatterns'))
addpath(fullfile(root, 'cumulus'))
addpath(fullfile(root, 'libra'))
addpath(fullfile(root, 'data'))
addpath(fullfile(root, 'density'))
addpath(genpath(fullfile(root, 'support')))

fprintf('OpenBreast v1.5 is ready\n')