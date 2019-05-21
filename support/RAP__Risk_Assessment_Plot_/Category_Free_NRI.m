function [cfNRI_output]=Category_Free_NRI(Reference, New_model,  bEvent,  varargin)
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
% Calculates the Category Free/Continuous Net Reclassification Improvement
% (NRI), see Ref 2 and the Integrated Discrimination Improvement (Ref 1)
% 
%
% INPUTS
% Reference:  The reference model (matrix of variables) eg APACHE II score, 
%  Age, Sex, etc ) with n rows where n is the number of participants.
% New_model:  Additional variables (eg one more urinary biomarker).  Array
%  or matrix
%
% bEvent:  Outcome vector (1=Event, 0=non event with 1 column).
% varargin:  Boolean vectors enabling choice of subcategories.  Eg Females
% and CKD patients.
%
% OUTPUTS
% cfNRI_output=[cfNRI_event, cfNRI_nonevent, cfNRI, IP_ref, IP_new, IS_ref,
%                      IS_new, IDI_event, IDI_nonevent, IDI, RelativeIDI ];
%
% cfNRI_event: The category free NRI for those with the event
% cfNRI_nonevent: The category free NRI for those without the event
% cfNRI= cfNRI_event + cfNRI_nonevent;
% IP_ref: Reference model Integrated 1-Specificity
% IP_new: New model Integrated 1-Specificity
% IS_ref: Reference model Integrated Sensitivity
% IS_new: New model Integrated Sensitivity
% IDI_event: Integrated Discrimination Improvement for those with the event
% IDI_nonevent: IDI for those without the event
% IDI=IDI_event + IDI_nonevent
% RelativeIDI:  Improvement in discrimination slope 
%               (IDI/Discrimination slope of reference)


% Calculate for each model the Pzs etc
if isempty(varargin)
    bchoice=ones(size(bEvent,1),1);
else
    bchoice=choice(varargin{:});
end

Y=[Reference, New_model];
Y_choice=zeros(2, size(Y,2));
for a=1:2
    for b=1:size(Y,2)-1;
        Y_choice(a,b)=1;
    end
end
Y_choice(1, size(Y,2))=0;
Y_choice(2, size(Y,2))=1;


for i=1:size(Y_choice,1) %ie each model (normally = 2)
    
    Y_choice_model=find(Y_choice(i,:)==1);
    
    [Zi, Pzi, B, stats] = LR_Pz_choice( Y,Y_choice_model,bEvent, bchoice );  %Logistic regression model incorporating only those covariates in Y_choice
    % Z(:,i)=Zi;
    Pz(:,i)=Pzi; %Probability for each patient for model i to have the outcome.
    % Beta(i)={B};
    % expBeta(i)={exp(B)};
    % CIplus(i)={exp(B+1.96*stats.se)};
    % CIminus(i)={exp(B-1.96*stats.se)};
    % CI95=strcat('(',num2str(CIminus,2)   ,' to ', num2str(CIplus,2) ,')');
end

bEvent=bEvent(find(bchoice==1));

n_event=sum(bEvent);
n_nonevent=sum(~bEvent);

% Probability movement (New model - old model)
dPz=Pz(:,2)-Pz(:,1);
n=n_event+n_nonevent;
n_up=size(find(dPz>0),1);                      %Number with increased risk
n_down=size(find(dPz<0),1);                 %Number with decreased risk
b_up=zeros(n,1); b_down=zeros(n,1);
b_up(find(dPz>0))=1;  b_down(find(dPz<0))=1;   %Individuals up or down

P_event=n_event/n;
P_nonevent=n_nonevent/n;
P_up=n_up/n;
P_down=n_down/n;
P_event_up=  (sum(bEvent.*b_up)/n)/P_up;
P_event_down=(sum(bEvent.*b_down)/n)/P_down;
P_nonevent_up=(sum((~bEvent).*b_up)/n)/P_up;
P_nonevent_down=(sum((~bEvent).*b_down)/n)/P_down;

% Category Free NRI metrics
cfNRI_event=100*(P_event_up*P_up-P_event_down*P_down)/P_event; %expressed as a percentage
cfNRI_nonevent=100*(P_nonevent_down*P_down-P_nonevent_up*P_up)/P_nonevent;
cfNRI=cfNRI_event+cfNRI_nonevent;

Events_ref_pz=Pz(find(bEvent),1);
mP_events_ref=nanmean(Events_ref_pz);
Events_new_pz=Pz(find(bEvent),2);
mP_events_new=nanmean(Events_new_pz);
Nonevents_new_pz=Pz(find(~bEvent),2);
mP_Nonevents_new=nanmean(Nonevents_new_pz);
Nonevents_ref_pz=Pz(find(~bEvent),1);
mP_Nonevents_ref=nanmean(Nonevents_ref_pz);

[f_ne_new, x_ne_new]=ecdf(Nonevents_new_pz);
[f_ne_ref, x_ne_ref]=ecdf(Nonevents_ref_pz);
[f_e_new, x_e_new]=ecdf(Events_new_pz);
[f_e_ref, x_e_ref]=ecdf(Events_ref_pz);

%Areas under the empirical cumulative distribution-risk curves
IP_new= trapz([0; x_ne_new; 1],[1; 1-f_ne_new; 0]);%Integrated 1-specificity
IP_ref=trapz([0; x_ne_ref; 1],[1; 1-f_ne_ref; 0]);
IS_new= trapz([0; x_e_new; 1],[1; 1-f_e_new; 0]); %integrated sensitivity
IS_ref= trapz([0; x_e_ref; 1],[1; 1-f_e_ref; 0]);

% IDI metrics
IDI_event=mP_events_new-mP_events_ref;
IDI_nonevent=mP_Nonevents_ref-mP_Nonevents_new;
%IDI_event=IS_new-IS_ref;    %Alternative form (See REFERENCES)
%IDI_nonevent=IP_ref-IP_new; %Alternative form
IDI=IDI_event+IDI_nonevent;
RelativeIDI=100*(((mP_events_new-mP_Nonevents_new)/(mP_events_ref-mP_Nonevents_ref))-1);

% Output
cfNRI_output=[cfNRI_event, cfNRI_nonevent, cfNRI, IP_ref, IP_new, IS_ref, IS_new, IDI_event, IDI_nonevent, IDI, RelativeIDI ];

end


