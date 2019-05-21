function [cfNRI_structure]=Category_Free_NRI_ci(Reference, New_model,  bEvent,  varargin)
% *********************  Category_Free_NRI_ci  ****************************
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

% PURPOSE
% Calculates the 95% confidence interval for the category-free NRI and IDI
% using the Category_Free_NRI function

bchoice=choice(varargin{:});
[cfNRI_output]=Category_Free_NRI(Reference, New_model,  bEvent,  bchoice);
tic
ci=bootci(2000,{@Category_Free_NRI,Reference, New_model,  bEvent,  bchoice}); % 2000 replications will result in an ~1% accuracy.  See  Hedges, SB Mol. Biol. Evol. 1992: 9, 366-369.    

cfNRI_structure.cfNRI_event=cfNRI_output(1);
cfNRI_structure.cfNRI_nonevent=cfNRI_output(2);
cfNRI_structure.cfNRI=cfNRI_output(3);
cfNRI_structure.IP_ref=cfNRI_output(4);
cfNRI_structure.IP_new=cfNRI_output(5);
cfNRI_structure.IS_ref=cfNRI_output(6);
cfNRI_structure.IS_new=cfNRI_output(7);
cfNRI_structure.IDI_event=cfNRI_output(8);
cfNRI_structure.IDI_nonevent=cfNRI_output(9);
cfNRI_structure.IDI=cfNRI_output(10);
cfNRI_structure.RelativeIDI=cfNRI_output(11);

toc
temp1=[num2str(cfNRI_output(1),3),' (',num2str(ci(1,1),3)  ,' to ', num2str(ci(2,1),3) ,')'];
cfNRI_structure.cfNRI_event_ci={temp1};
temp2=[num2str(cfNRI_output(2),3),' (',num2str(ci(1,2),3)  ,' to ', num2str(ci(2,2),3) ,')'];
cfNRI_structure.cfNRI_nonevent_ci={temp2};
temp3=[num2str(cfNRI_output(3),3),' (',num2str(ci(1,3),3)  ,' to ', num2str(ci(2,3),3) ,')'];
cfNRI_structure.cfNRI_ci={temp3};
temp4=[num2str(cfNRI_output(4),3),' (',num2str(ci(1,4),3)  ,' to ', num2str(ci(2,4),3) ,')'];
cfNRI_structure.IP_ref_ci={temp4};
temp5=[num2str(cfNRI_output(5),3),' (',num2str(ci(1,5),3)  ,' to ', num2str(ci(2,5),3) ,')'];
cfNRI_structure.IP_new_ci={temp5};
temp6=[num2str(cfNRI_output(6),3),' (',num2str(ci(1,6),3)  ,' to ', num2str(ci(2,6),3) ,')'];
cfNRI_structure.IS_ref_ci={temp6};
temp7=[num2str(cfNRI_output(7),3),' (',num2str(ci(1,7),3)  ,' to ', num2str(ci(2,7),3) ,')'];
cfNRI_structure.IS_new_ci={temp7};
temp8=[num2str(cfNRI_output(8),3),' (',num2str(ci(1,8),3)  ,' to ', num2str(ci(2,8),3) ,')'];
cfNRI_structure.IDI_event_ci={temp8};
temp9=[num2str(cfNRI_output(9),3),' (',num2str(ci(1,9),3)  ,' to ', num2str(ci(2,9),3) ,')'];
cfNRI_structure.IDI_nonevent_ci={temp9};
temp10=[num2str(cfNRI_output(10),3),' (',num2str(ci(1,10),3)  ,' to ', num2str(ci(2,10),3) ,')'];
cfNRI_structure.IDI_ci={temp10};
temp11=[num2str(cfNRI_output(11),3),' (',num2str(ci(1,11),3)  ,' to ', num2str(ci(2,11),3) ,')'];
cfNRI_structure.RelativeIDI_ci={temp11};
end
