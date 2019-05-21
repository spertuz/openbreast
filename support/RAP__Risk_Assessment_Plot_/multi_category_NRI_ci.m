function [multi_category_NRI_structure]=multi_category_NRI_ci( Reference, New_model,  Risk_groups, bEvent,  varargin)
% *********************  multi_category_NRI_ci  ****************************
%   (c) John W Pickering, August 2011
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
% Attribution to John Pickering.  Publications to reference
%  Pickering JW, Endre ZH. 
%   New Metrics for Assessing Diagnostic Potential of Candidate Biomarkers. 
%   Clin J Am Soc Nephro 2012, On line ahead of print, doi:10.2215/CJN.09590911;
%   http://cjasn.asnjournals.org/content/early/2012/05/24/CJN.09590911.full.pdf+html
%    (Published Open Access)
% *************************************************************************

% PURPOSE% Calculates the 95% confidence interval for the multi_category NRI

bchoice=choice(varargin{:});
[NRI_output]=multi_category_NRI( Reference, New_model,  Risk_groups, bEvent,  bchoice);
tic 
ci=bootci(2000,{@multi_category_NRI, Reference, New_model,  Risk_groups, bEvent,  bchoice}); % 2000 replications will result in an ~1% accuracy.  See  Hedges, SB Mol. Biol. Evol. 1992: 9, 366-369.    
 
NRI_structure.NRI_event=NRI_output(1);
NRI_structure.NRI_nonevent=NRI_output(2);
NRI_structure.NRI=NRI_output(3);

toc
temp1=[num2str(NRI_output(1),3),' (',num2str(ci(1,1),3)  ,' to ', num2str(ci(2,1),3) ,')'];
multi_category_NRI_structure.NRI_event_ci={temp1};
temp2=[num2str(NRI_output(2),3),' (',num2str(ci(1,2),3)  ,' to ', num2str(ci(2,2),3) ,')'];
multi_category_NRI_structure.NRI_nonevent_ci={temp2};
temp3=[num2str(NRI_output(3),3),' (',num2str(ci(1,3),3)  ,' to ', num2str(ci(2,3),3) ,')'];
multi_category_NRI_structure.NRI_ci={temp3};

end
