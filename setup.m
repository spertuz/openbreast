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
addpath(genpath(fullfile(root, 'support')))

fprintf('OpenBreast v0.1 is ready\n')