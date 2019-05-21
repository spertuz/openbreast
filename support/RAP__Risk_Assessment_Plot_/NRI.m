function [NRI_structure]=NRI(Model_Names, Reference, New_model,  Risk_groups, bEvent,  varargin)
% *********************  NRI  ****************************
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
%  3. DeLong E, DeLong D, Clarke-Pearson D. 
%   Comparing the areas under 2 or more correlated receiver operating 
%   characteristic curves - a nonparametric approach. 
%   Biometrics 1988;44(3):837?45. 
% *************************************************************************

% PURPOSE
% Calculates the Net Reclassification Index (NRI) and IDI for each possible
% pair of models (See Ref 1).  Calculates the Area Under the Reciever
% Operator Characteristic Curve (AUC) for both the reference and new model
% and compares them using the De Long method (Ref 3).

% INPUTS
% Model_Names: Cell array of model names
% Reference:  The reference model (matrix of variables) eg APACHE II score, Age, Sex, Log(urinary
% biomarker) etc ) with n rows where n is the number of participants.
% New_model:  Additional variables (eg one more urinary biomarker).  Array
% or matrix
% Risk_groups:  Thresholds for risk groups eg [0.2, 0.4] defines three
% groups <0.2 (low), 0.2 to <0.4 (medium), >=0.4 (high)
% bEvent:  Outcome vector (1=Event, 0=non event with 1 column).
% varargin:  Boolean vectors enabling choice of subcategories.  Eg Females
% and CKD patients.

% OUTPUTS
% Structure containing
%   Positive: Reclassification table for those with the event (final column
%   and final rows are summations)
%   Negative:  Reclassification table for those without the event
%   bOut: Vector of 1=subject had event, 0=subject didn't have event for
%   the subjects in the chosen cohort.
%   Combination: Cell array of model names
%   Beta: Cell array of the beta coefficients of each logistic regression
%   model for the reference and new models.
%   expBeta: exponential of beta coefficients (odds ratio)
%   NRI: Net Reclassification Improvement
%   NRI_events: Net Reclassification Improvement for those with the event
%   NetChangePerClass_Events: Net (%) change in reclassification for those
%   from each risk group with event 
%   NRI_nonevents: Net Reclassification Improvement for those without the
%   event
%   NetChangePerClass_nonEvents: Net (%) change in reclassification for those
%   from each risk group without event
%   Pup_events: Proportion of those with the event moving to a higher risk
%   group
%   Pdown_events: Proportion of those with the event moving to a lower risk
%   group
%   Pup_nonevents:Proportion of those without the event moving to a higher risk
%   group
%   Pdown_nonevents: Proportion of those without the event moving to a lower risk
%   group
%   Nup_events: Number of those with the event moving to a higher risk
%   group
%   Ndown_events: Number of those with the event moving to a lower risk
%   group
%   Nup_nonevents: Number of those without the event moving to a higher risk
%   group
%   Ndown_nonevents: Number of those without the event moving to a lower risk
%   group
%   AUC: Cell array of the AUCs for each model, the difference in the AUC
%   and a p value for the difference calculated using the DeLong method
%   (Ref 3).
%   Pz: Probability of having the event calculated by the reference model
%   (column 1) and the new model (column 2)
%   Risk: Risk group calculated by the reference model(column 1) and the
%   new model (column 2).  1 is the lowest risk group.

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


%Logisitic regression for each model
for i=1:size(Y_choice,1) 
    
    Y_choice_model=find(Y_choice(i,:)==1);
    
    [Zi, Pzi, B, stats] = LR_Pz_choice( Y,Y_choice_model,bEvent, bchoice );  
    Z(:,i)=Zi;
    Pz(:,i)=Pzi; %Probability for each patient for model i to have the outcome.
    Beta(i)={B};
    expBeta(i)={exp(B)};
    CIplus(i)={exp(B+1.96*stats.se)};
    CIminus(i)={exp(B-1.96*stats.se)};
    
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

R=Risk;

% for each model construct 2 matrices compared with each other model.

bOut=bEvent(find(bchoice));

c=combnk(1:size(Y_choice,1),2);
c=sortrows(c,[1,2]);
for z=1:size(c,1)
    
    R1_index=find(R(:,c(z,1))>0);
    R2_index=find(R(:,c(z,2))>0);
    R_index=intersect(R1_index,R2_index);
    bOut_valid=bOut(R_index);
    bOut_Disease_index=find(bOut_valid==1);
    n_disease=size(bOut_Disease_index,1);
    bOut_NotDisease_index=find(bOut_valid==0);
    n_Notdisease=size(bOut_NotDisease_index,1);
    R_temp=R(R_index,:);
    
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
    NRI_non_events=(Pdown_nonevents-Pup_nonevents);
    
    
    % Sum in each cell in an n_risk X n_risk matrix
    Positive_Out=zeros(n_risk+1,n_risk+1);
    Negative_Out=zeros(n_risk+1,n_risk+1);
    
    for j=1:size(R,1)
        if R(j,c(z,1))> 0 & R(j,c(z,2)) > 0
            if bOut(j)==1
                Positive_Out(R(j,c(z,1)),R(j,c(z,2)))=Positive_Out(R(j,c(z,1)),R(j,c(z,2)))+1;
            else
                Negative_Out(R(j,c(z,1)),R(j,c(z,2)))=Negative_Out(R(j,c(z,1)),R(j,c(z,2)))+1;
            end
        end
    end
    
    %Complete the classification tables with the sums
    for k=1:n_risk
        for l=1:n_risk
            Positive_Out(l,n_risk+1)=Positive_Out(l,n_risk+1)+Positive_Out(l,k);
            Positive_Out(n_risk+1,k)=Positive_Out(n_risk+1,k)+Positive_Out(l,k);
            Negative_Out(l,n_risk+1)=Negative_Out(l,n_risk+1)+Negative_Out(l,k);
            Negative_Out(n_risk+1,k)=Negative_Out(n_risk+1,k)+Negative_Out(l,k);
        end
    end
    Positive_Out(n_risk+1,n_risk+1)=sum(Positive_Out(:,n_risk+1));
    Negative_Out(n_risk+1,n_risk+1)=sum(Negative_Out(:,n_risk+1));
    
    for m=1:n_risk
        NetChangePerClass_Events(m)=100*(sum(Positive_Out(1:n_risk,m))-sum(Positive_Out(m,1:n_risk)))/ Positive_Out(n_risk+1,n_risk+1);
        NetChangePerClass_Non_Events(m)=100*(sum(Negative_Out(1:n_risk,m))-sum(Negative_Out(m,1:n_risk)))/ Negative_Out(n_risk+1,n_risk+1);
    end
    
    % AUC comparison
    [AUC_compare]=AUC_compare_correlated([Pz(:,c(z,1)), Pz(:,c(z,2))], bEvent(find(bchoice)), ones(size(bEvent(find(bchoice)),1),1));
    
    NRI_structure.Postive(:,:,z)=Positive_Out;
    NRI_structure.Negative(:,:,z)=Negative_Out;
    NRI_structure.bOut=bOut;
    
    NRI_structure.Combination(z,:)=strcat(Model_Names(c(z,1)),'-', Model_Names(c(z,2)));
    NRI_structure.Beta=Beta;
    NRI_structure.expBeta=expBeta;
    
    NRI_structure.NRI(:,z)=NRI;
    
    NRI_structure.NRI_events(:,z)=NRI_events;
    
    NRI_structure.NetChangePerClass_Events=NetChangePerClass_Events;
    
    NRI_structure.NRI_non_events(:,z)=NRI_non_events;
    
    NRI_structure.NetChangePerClass_Non_Events=NetChangePerClass_Non_Events;
    
    NRI_structure.Pup_events(:,z)=Pup_events;
    NRI_structure.Pdown_events(:,z)=Pdown_events;
    NRI_structure.Pup_nonevents(:,z)=Pup_nonevents;
    NRI_structure.Pdown_nonevents(:,z)=Pdown_nonevents;
    NRI_structure.Nup_events(:,z)=Nup_events;
    NRI_structure.Ndown_events(:,z)=Ndown_events;
    NRI_structure.Nup_nonevents(:,z)=Nup_nonevents;
    NRI_structure.Ndown_nonevents(:,z)=Ndown_nonevents;
    NRI_structure.n_disease(:,z)=n_disease;
    NRI_structure.n_Notdisease(:,z)=n_Notdisease;
    
    NRI_structure.AUC(:,z)= AUC_compare;
end

%NRI_structure.names=Model_Names;
NRI_structure.Pz=Pz;
NRI_structure.Risk=R;

end
