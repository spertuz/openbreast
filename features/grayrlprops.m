function [SRE, LRE, GLN, RLN, RP, LGRE, HGRE] = grayrlprops(g, np)
% Inputs:
%     g,      cell array as returned by grayrlmatrix
%     np,     total number of pixels in the image
%     
% Modified from: https://se.mathworks.com/matlabcentral/fileexchange/17482-gray-level-run-length-matrix-toolbox
%  -------------------------------------------
%  Current supported statistics include:
%  -------------------------------------------
%   Short Run Emphasis (SRE)
%   Long Run Emphasis (LRE)
%   Gray-Level Nonuniformity (GLN)
%   Run Length Nonuniformity (RLN)
%   Run Percentage (RP)
%   Low Gray-Level Run Emphasis (LGRE)
%   High Gray-Level Run Emphasis (HGRE)
%   Short Run Low Gray-Level Emphasis (SRLGE)
%   Short Run High Gray-Level Emphasis (SRHGE)
%   Long Run Low Gray-Level Emphasis (LRLGE)
%   Long Run High Gray-Level Emphasis (LRHGE)
%  --------------------------------------------
%  Reference:
%  --------------------------------------------
%   Xiaoou Tang,Texture Information in Run-Length Matrices
%   IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL.7, NO.11,NOVEMBER 1998
% ---------------------------------------------
%  See also GRAYRLMATRIX.
% ---------------------------------------------
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% ---------------------------------------------
% History:
% ---------------------------------------------
% Creation: beta         Date: 01/10/2007
% Revision: 1.0          Date: 12/11/2007
% 1.Accept cell input now
% 2.Using MATLAB file style
% 3.Fully vectorized programming
% 4.Fully support the IEEE reference
% 5. ...


% Initialization
SRE = zeros(length(g), 1);
LRE = zeros(length(g), 1);
GLN = zeros(length(g), 1);
RLN = zeros(length(g), 1);
RP = zeros(length(g), 1);
LGRE = zeros(length(g), 1);
HGRE = zeros(length(g), 1);


for n = 1 : length(g)
    
    %Current GLRL matrix:
    p = g{n};
    
    %Row and column indici
    [j, i] = meshgrid(1:size(p, 2), 1:size(p, 1));
    
    % Total number of runs
    N_runs = sum(sum(p));

        %------------------------Statistics-------------------------------
    % 1. Short Run Emphasis (SRE)
    SRE(n) = (1/N_runs)*sum(sum(p./(j.^2), 2), 1);
    
    % 2. Long Run Emphasis (LRE)
    LRE(n) = (1/N_runs)*sum(sum(p.*(j.^2), 2), 1);
    
    % 3. Gray-Level Nonuniformity (GLN)
    GLN(n) = (1/N_runs)*sum(sum(p, 2).^2, 1);
    
    % 4. Run Length Nonuniformity (RLN)
    RLN(n) = (1/N_runs)*sum(sum(p, 1).^2, 2);
    
    % 5. Run Percentage (RP)
    RP(n) = N_runs/np;
    
    % 6. Low Gray-Level Run Emphasis (LGRE)
    LGRE(n) = (1/N_runs)*sum(sum(p./(i.^2), 2), 1);    
        
    % 7. High Gray-Level Run Emphasis (HGRE)
    HGRE(n) = (1/N_runs)*sum(sum(p.*(i.^2), 2), 1);
    
%     % 8. Short Run Low Gray-Level Emphasis (SRLGE)
%     SGLGE =calculate_SGLGE(tGLRLM,r_matrix',c_matrix',N_runs);
%     % 9. Short Run High Gray-Level Emphasis (SRHGE)
%     SRHGE =calculate_SRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
%     % 10. Long Run Low Gray-Level Emphasis (LRLGE)
%     LRLGE =calculate_LRLGE(tGLRLM,r_matrix',c_matrix',N_runs);
%     % 11.Long Run High Gray-Level Emphasis (LRHGE
%     LRHGE =calculate_LRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
%     %----------------insert statistics----------------------------
%     stats(p,:)=[SRE LRE GLN RLN  RP LGRE HGRE SGLGE SRHGE LRLGE  LRHGE ];
end % end all run length matrixs

%   ----------------------Utility functions--------------------
%-----------------------------------------------------------------------------
function SGLGE =calculate_SGLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run Low Gray-Level Emphasis (SRLGE):

term = tGLRLM./((r_matrix.*c_matrix).^2);
SGLGE= sum(sum(term))./N_runs;

%------------------------------------
function  SRHGE =calculate_SRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run High Gray-Level Emphasis (SRHGE):
%
term  = tGLRLM.*(r_matrix.^2)./(c_matrix.^2);
SRHGE = sum(sum(term))/N_runs;
%------------------------------------
function   LRLGE =calculate_LRLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run Low Gray-Level Emphasis (LRLGE):
%
term  = tGLRLM.*(c_matrix.^2)./(r_matrix.^2);
LRLGE = sum(sum(term))/N_runs;
%---------------------------------------
function  LRHGE =calculate_LRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run High Gray-Level Emphasis (LRHGE):
%
term  = tGLRLM.*(c_matrix.^2).*(r_matrix.^2);
LRHGE = sum(sum(term))/N_runs;
%----------------------------------------