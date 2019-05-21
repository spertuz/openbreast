
function RAP_stats=Risk_Assessment_Plot(Model_Names, Reference, BM, Risk, bEvent,  varargin)
% *********************  Risk_Assessment_Plot  ****************************
%   (c) John W Pickering, August 2011
%     All rights reserved.
%     Christchurch Kidney Research Group
%     University of Otago Christchurch
%     New Zealand
%
%   Last update:  15 August 2012
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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% 
% Attribution to John Pickering.  Publications to reference
%  Pickering JW, Endre ZH.
%   New Metrics for Assessing Diagnostic Potential of Candidate Biomarkers.
%   Clin J Am Soc Nephro 2012, On line ahead of print, doi:10.2215/CJN.09590911;
%   http://cjasn.asnjournals.org/content/early/2012/05/24/CJN.09590911.full.pdf+html
%    (Published Open Access)
%
%   PLEASE NOTE:  The Risk Assessment Plots are based on a combination of several
%                 classification plots (see ref 4) and were first described and
%                 published in the above publication. That publication is
%                 also an introduction to the appropriate application of 
%                 the NRI and IDI
%                 The NRI and IDI were first described in the two Pencina
%                 publications below
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
%  4.Pepe MS, Feng Z, Huang Y, et al. 
%   Integrating the predictiveness of a marker with its performance as a classifier.
%   Am J Epidemiol 2008;167(3):362?8. 
% *************************************************************************

% PURPOSE
% Calculates the Net Reclassification Index (NRI) for a two category model
% with cutpoint at Risk (Ref 1), the continuous or category free NRI (Ref
% 2) and the IDI (Ref 1). Also calculates the AUC for each model and the
% difference between (Ref 3).  Confidence intervals for the NRI and IDI are
% generated using a boot strap method.
% Finally this plots a Risk Assessment Plot.  


% INPUTS
% Model_Names: Cell array of model names
% Reference:  The reference model (matrix of variables) eg APACHE II score, 
% Age, Sex, Log(urinary biomarker) etc ) with n rows where n is the number 
% of participants. Each column represents a variable
% New_model:  Additional variables (eg one more urinary biomarker). 
% Risk:  Thresholds for a two groups (low, high).  eg 0.4
% bEvent:  Outcome vector (1=Event, 0=non event with 1 column).
% varargin:  Boolean vectors enabling choice of subcategories.  Eg Females
% and CKD patients.

% OUTPUTS
% A Risk Assessment Plot
% Structure containing
%   Reference=Model_Names(1);
%   New=Model_Names(2);
%   IDI_events=IDI for those with the event (95% CI)
%   IDI_nonevents=IDI for those without the event (95% CI)
%   IDI=IDI_events+IDI_nonevents (95% CI)
%   RelativeIDI=IDI releative to the reference discrimination curve;
%   IS_ref=Integrated Sensitivity of the reference model
%   IS_new=Integrated Sensitivity of the new model
%   IP_ref=Integrated 1-Specificity of the reference model
%   IP_new=Integrated 1-Specificity of the new model
%   cfNRI_events=category free/continuous NRI for those with the event (95% CI)
%   cfNRI_nonevents=category free/continuous NRI for those without the event (95% CI)
%   cfNRI=cfNRI_events+cfNRI_nonevents (95% CI)
%   twoCatNRI_events=two category NRI (<Risk, >=Risk)for those with the event (95% CI)
%   twoCatNRI_nonevents=two category NRI (<Risk, >=Risk) for those without the event (95% CI)
%   twoCatNRI=twoCatNRI_events+twoCatNRI_nonevents (95% CI)
%   twoCatThreshold=Risk
%   AUC: Areas under the receiver operator characteristic curve for the
%   reference and for the new models.  Difference in the areas. p value for
%   the difference using the DeLong method (Ref 3).
%   Nonevents_new_plotdata= Reference model Calculated risk (Col 1) v 1-Specificity (Col 2) 
%   Nonevents_ref_plotdata=New model Calculated risk (Col 1) v 1-Specificity (Col 2)
%   Events_new_plotdata=Reference model Calculated risk (Col 1) v Sensitivity (Col 2)
%   Events_ref_plotdata= New modelCalculated risk (Col 1) v Sensitivity (Col 2)
%   Risk_probabilities=Probability of event for each model (Col 1 = Ref, Col2 =New);
%   two_category_NRI_structure=output of multi_category_NRI for the two category model (see multi_category_NRI)


% Risk = sum(bEvent.*bchoice)/size(bEvent.*bchoice,1); % Option: Set the Risk to the prevalence

% Select the subcohort.
if isempty(varargin)
    bchoice=ones(size(BM,1),1);
else
bchoice=choice(varargin{:});
Reference=Reference(find(bchoice==1),:); BM=BM(find(bchoice==1)); bEvent=bEvent(find(bchoice==1)); bchoice=bchoice(find(bchoice==1));
end

% Calculate the basic NRI information including tables the AUCs and the Probabilities
[two_category_NRI_structure]=NRI(Model_Names, Reference, BM,  Risk, bEvent,  bchoice);

% Calculate the NRI for a two category model at the prevelance with a 95% CI (bootstrap)
[two_category_NRI_structure_ci]=multi_category_NRI_ci( Reference, BM,  Risk, bEvent,  bchoice);

% Calculate the NRI and IDI for a category free model at the prevelance with a 95% CI
[cfNRI_output_ci]=Category_Free_NRI_ci(Reference, BM,  bEvent,  bchoice);

% Prepare data for plotting
Pz=two_category_NRI_structure.Pz;  %Change for the appropriate NRI array

Events_ref_pz=Pz(find(bEvent),1);
Events_new_pz=Pz(find(bEvent),2);
Nonevents_new_pz=Pz(find(~bEvent),2);
Nonevents_ref_pz=Pz(find(~bEvent),1);
[f_ne_new, x_ne_new]=ecdf(Nonevents_new_pz);
[f_ne_ref, x_ne_ref]=ecdf(Nonevents_ref_pz);
[f_e_new, x_e_new]=ecdf(Events_new_pz);
[f_e_ref, x_e_ref]=ecdf(Events_ref_pz);

%Make sure that starts at 0,1
x_ne_new=vertcat(0, x_ne_new, 1);
x_e_new=vertcat(0, x_e_new, 1);
x_ne_ref=vertcat(0, x_ne_ref, 1);
x_e_ref=vertcat(0, x_e_ref, 1);
f_ne_new=vertcat(0, f_ne_new, 1);
f_e_new=vertcat(0, f_e_new, 1);
f_ne_ref=vertcat(0, f_ne_ref, 1);
f_e_ref=vertcat(0, f_e_ref, 1);


Nonevents_new_plotdata=[x_ne_new, 1-f_ne_new];  %ie Risk v 1- empirical cumulative distribution
Nonevents_ref_plotdata=[x_ne_ref, 1-f_ne_ref];
Events_new_plotdata=[x_e_new, 1-f_e_new];
Events_ref_plotdata=[x_e_ref, 1-f_e_ref];

% My Plot

stairs(Events_ref_plotdata(:,1), Events_ref_plotdata(:,2),'k--','markersize',4,'linewidth',1.5)
hold on
mycolors=colormap(lines(7));
stairs(Events_new_plotdata(:,1), Events_new_plotdata(:,2),'k-','markersize',4,'linewidth',1.5)
stairs(Nonevents_ref_plotdata(:,1), Nonevents_ref_plotdata(:,2),'r--','markersize',4,'linewidth',1.5)
stairs(Nonevents_new_plotdata(:,1), Nonevents_new_plotdata(:,2),'r-','markersize',4,'linewidth',1.5)

set(gca,'FontSize',14,'FontWeight','bold')
xlabel('Risk','FontSize',16,'FontWeight','bold')
ylabel('Sensitivity, 1-Specificity','FontSize',16,'FontWeight','bold')
legend('Reference: Events', 'Reference + Biomarker: Events','Reference: Non-events','Reference + Biomarker: Non-events', 'Location','NE')
set(gca,'Units','points');
rect= get(gca,'Position');
rect2=[rect(1), rect(2), 1.35*rect(3), 1.35*rect(4)];
set(gca,'Position',rect2);
fig_pos= get(gcf,'Position');
fig_pos_new=[fig_pos(1), fig_pos(2), 1.3*fig_pos(3), 1.3*fig_pos(4)];
set(gcf,'Position',fig_pos_new);
hold off

RAP_stats.Reference=Model_Names(1);
RAP_stats.New=Model_Names(2);
RAP_stats.IDI_events=cfNRI_output_ci.IDI_event_ci;
RAP_stats.IDI_nonevents=cfNRI_output_ci.IDI_nonevent_ci;
RAP_stats.IDI=cfNRI_output_ci.IDI_ci;
RAP_stats.RelativeIDI=cfNRI_output_ci.RelativeIDI_ci;
RAP_stats.IS_ref=cfNRI_output_ci.IS_ref_ci;
RAP_stats.IS_new=cfNRI_output_ci.IS_new_ci;
RAP_stats.IP_ref=cfNRI_output_ci.IP_ref_ci;
RAP_stats.IP_new=cfNRI_output_ci.IP_new_ci;
RAP_stats.cfNRI_event=cfNRI_output_ci.cfNRI_event_ci;
RAP_stats.cfNRI_nonevent=cfNRI_output_ci.cfNRI_nonevent_ci;
RAP_stats.cfNRI=cfNRI_output_ci.cfNRI_ci;
RAP_stats.twoCatNRI_event=two_category_NRI_structure_ci.NRI_event_ci;
RAP_stats.twoCatNRI_nonevent=two_category_NRI_structure_ci.NRI_nonevent_ci;
RAP_stats.twoCatNRI=two_category_NRI_structure_ci.NRI_ci;
RAP_stats.twoCatThreshold=Risk;
RAP_stats.AUC=two_category_NRI_structure.AUC;
RAP_stats.Nonevents_new_plotdata=[x_ne_new, 1-f_ne_new];  
RAP_stats.Nonevents_ref_plotdata=[x_ne_ref, 1-f_ne_ref];
RAP_stats.Events_new_plotdata=[x_e_new, 1-f_e_new];
RAP_stats.Events_ref_plotdata=[x_e_ref, 1-f_e_ref];
RAP_stats.Risk_probabilities=Pz;
RAP_stats.two_category_NRI_structure=two_category_NRI_structure;


end
