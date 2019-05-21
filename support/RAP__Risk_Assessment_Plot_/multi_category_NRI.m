function [NRI_output]=multi_category_NRI( Reference, New_model,  Risk_groups, bEvent,  varargin)
% *********************  multi_category_NRI  ****************************
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
%
%
%  REFERENCES
%  1. Pencina, MJ et al. 
%   Evaluating the added predictive ability of a new marker: From area under
%   the ROC curve to reclassification and beyond. 
%   Stat Med 2008; 27:157-172. doi:10.1002/sim.2929
%  2. Pencina et al. 
%   Extensions of net reclassification improvement calculations to measure 
%   usefulness of new biomarkers. 
%   Stat Ned 2011; 30:11-21. doi:10.1002/sim.4085
% *************************************************************************

% PURPOSE
% Calculates the Net Reclassification Improvement(NRI), see Ref 1, for any
% number of risk groups

% INPUTS
% Reference:  The reference model (matrix of variables) eg APACHE II score, Age, Sex, Log(urinary
% biomarker) etc ) with n rows where n is the number of subjects.
% New_model:  Additional variables (eg one more urinary biomarker). 
% Risk_groups:  risk thresholds eg Risk_groups=[0.2,0.7] defines three groups "low" with risk
% probability<0.2, "Medium" 0.2 to <0.7, and "High" >0.7
% bEvent:  Outcome vector (1=Event, 0=non event with 1 column).
% varargin:  Boolean vectors enabling choice of subcategories.  Eg bFemales, bOlderthan60

% OUTPUTS
% NRI_events: the NRI for those with the event (as a percentage)
% NRI_nonevents: the NRI for those without the event (as a percentage)
% NRI=NRI_events+NRI_nonevents

% boolean vector for selecting a cohort of subjects
if isempty(varargin)
    bchoice=ones(size(bEvent,1),1);
else
    bchoice=choice(varargin{:});
end

% Calculate for each model the Pzs etc
Y=[Reference, New_model];
Y_choice=zeros(2, size(Y,2));
for a=1:2
    for b=1:size(Y,2)-1;
        Y_choice(a,b)=1;
    end
end
Y_choice(1, size(Y,2))=0;
Y_choice(2, size(Y,2))=1;

for i=1:size(Y_choice,1) %ie each model
    
    Y_choice_model=find(Y_choice(i,:)==1);
    
    [Zi, Pzi, B, stats] = LR_Pz_choice( Y,Y_choice_model,bEvent, bchoice );  %Logistic regression model incorporating only those covariates in Y_choice
    Z(:,i)=Zi;
    Pz(:,i)=Pzi; %Probability for each patient for model i to have the outcome.
    
end

% how many risk groups?
n_risk=size(Risk_groups,2)+1;
Risk_class=1:n_risk;

% find which risk group in
Risk=zeros(size(Pz,1),size(Pz,2));

Risk(find(Pz<Risk_groups(1)))=1;
Risk(find(Pz>=Risk_groups(n_risk-1)))=n_risk;
for j=2:n_risk-1;
    Risk(find(Pz>=Risk_groups(j-1) & Pz<Risk_groups(j)))=j;
end

% For each model construct 2 matrices compared with each other model.

bOut=bEvent(find(bchoice==1));

c=combnk(1:size(Y_choice,1),2);
c=sortrows(c,[1,2]);
for z=1:size(c,1)
    
    R1_index=find(Risk(:,c(z,1))>0);
    R2_index=find(Risk(:,c(z,2))>0);
    R_index=intersect(R1_index,R2_index);
    bOut_valid=bOut(R_index);
    bOut_Disease_index=find(bOut_valid==1);
    n_disease=size(bOut_Disease_index,1);
    bOut_NotDisease_index=find(bOut_valid==0);
    n_Notdisease=size(bOut_NotDisease_index,1);
    R_temp=Risk(R_index,:);
    
    %Calculate the NRIs (See Pencina et al Stat Med 2008; 27:157-172)
    Nup_events=size(find(R_temp(bOut_Disease_index,c(z,2))>R_temp(bOut_Disease_index,c(z,1))),1);
    Pup_events=Nup_events/n_disease;
    
    Ndown_events=size(find(R_temp(bOut_Disease_index,c(z,2))<R_temp(bOut_Disease_index,c(z,1))),1);
    Pdown_events=Ndown_events/n_disease;
    
    Nup_nonevents=size(find(R_temp(bOut_NotDisease_index,c(z,2))>R_temp(bOut_NotDisease_index,c(z,1))),1);
    Pup_nonevents=Nup_nonevents/n_Notdisease;
    
    Ndown_nonevents=size(find(R_temp(bOut_NotDisease_index,c(z,2))<R_temp(bOut_NotDisease_index,c(z,1))),1);
    Pdown_nonevents=Ndown_nonevents/n_Notdisease;
    
    NRI=(Pup_events-Pdown_events)-(Pup_nonevents-Pdown_nonevents);
    
    % EVENTS
    NRI_events=(Pup_events-Pdown_events);
    
    % NON-EVENTS
    NRI_nonevents=(Pdown_nonevents-Pup_nonevents);
    
    NRI_output=[100*NRI_events, 100*NRI_nonevents, 100*NRI];
    
end
