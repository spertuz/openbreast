function [SE, CIMin95, CIPlus95, Area]=CIAUC(C,D,A);
% *********************  CIAUC  ****************************
%   (c) John W Pickering, Novemeber 2009
%     Christchurch Kidney Research Group
%     University of Otago Christchurch
%     New Zealand
%
%   Last update:  17 July 2012
%
%	Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in 
%     the documentation and/or other materials provided with the distribution
%
% Attribution to John Pickering.  
% *************************************************************************

% PURPOSE
% Calculates the Standard Error (SE)  and the  95% Confidence interval 
% (CIMin to CIPlus) of an Area under the curve with AUC of A (=Area).
%
% INPUTS
% C is the number of controls/normals
% D is the number of diseased/abnormals
% A is the AUC

% OUTPUTS
% SE: Standard Error
% CIMin95:  95% CI for Standard Error
% CIPlus95:  95% CI for Standard Error
% Area: A

Q1=A/(2-A);
Q2=2*A*A/(1+A);

SE=sqrt((A*(1-A)+(D-1)*(Q1-A*A)+(C-1)*(Q2-A*A))/(C*D));

CIMin95=A-1.96*SE;
CIPlus95=A+1.96*SE;
Area=A;
end
