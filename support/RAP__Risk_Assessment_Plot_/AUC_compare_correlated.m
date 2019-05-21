function [AUC_compare]=AUC_compare_correlated(bm,  Disease, varargin)
% *********************  AUC_compare_correlated  ****************************
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
%
%  REFERENCE
%   DeLong E, DeLong D, Clarke-Pearson D. 
%   Comparing the areas under 2 or more correlated receiver operating 
%   characteristic curves - a nonparametric approach. 
%   Biometrics 1988;44(3):837?45. 
% *************************************************************************


% PURPOSE
% function to compare two correlated ROC curves. (ie when based on same
% sample of cases - eg comparing two analytes for predictive value of the
% same outcome (eg Sepsis)

% INPUTS
% bm = matrix of biomarkers two (ie two columns)
% bm_titles= Titles of biomarkers
% Disease  = defines diseases and not diseased
% varargin = defines cohort

% OUTPUS
% cell array
% AUCs with 95% confidence intervals for each column in bm
% The difference in AUCs with 95% confidence interval
% p value for the difference in AUCs


%excise any rows with NaN
n_bm=size(bm,2); % number of biomarkers

choice_N=~Disease;
choice_D=Disease;

% Choice
% choose normals
% choose the patients
if isempty(varargin)
    bchoice=ones(size(Disease,1),1);   
else
    bchoice=choice(varargin{:});
choice_D=choice_D(find(bchoice));
choice_N=choice_N(find(bchoice));
end

% Filter out the NaNs
choice_N(any(isnan(bm),2),:)=[];
choice_D(any(isnan(bm),2),:)=[];
bm=exciseRows(bm);

% FIND normals and diseases (vectors)
for i=1:n_bm
    [ N_temp, D_temp ] = FindNandD( bm(:,i), choice_N, choice_D );
    Normals(:,i)=N_temp;
    Diseased(:,i)=D_temp;
    
end

n_D=size(Diseased,1); % number of Diseased
n_N=size(Normals,1); % number of Normals

% compute the structural components (perhaps this can be done quicker, but,
% hey, it works!)
for p=1:n_bm
    for r=1:n_D
        s_temp=0;
        for q=1:n_N
            if Diseased(r,p)>Normals(q,p)
                s_temp=s_temp+1;
            end
            if  Diseased(r,p)==Normals(q,p)
                s_temp=s_temp+0.5;
            end
        end
        V_D(r,p)=s_temp/n_N;
    end
    for r=1:n_N
        s_temp=0;
        for q=1:n_D
            if Normals(r,p)<Diseased(q,p)
                s_temp=s_temp+1;
            end
            if Normals(r,p)==Diseased(q,p)
                s_temp=s_temp+0.5;
            end
        end
        V_N(r,p)=s_temp/n_D;
    end
end

%compute the estimates of the areas
for p=1:n_bm
    s_temp=0;
    for l=1:n_D
        s_temp=s_temp+V_D(l,p);
    end
    Area(p)=s_temp/n_D;
    [SE(p), CIMin95(p), CIPlus95(p), Area_notused]=CIAUC(n_N,n_D,Area(p));
end


%compute the estimated covariance matrix of the area vector
for k=1:n_bm
    for m=k:n_bm
        s_temp=0;
        for i=1:n_D
            s_temp=s_temp+((V_D(i,k)-Area(k))*(V_D(i,m)-Area(m)));
        end
        S10M(k,m)=s_temp/(n_D-1);
        S10M(m,k)=S10M(k,m);
    end
end
for k=1:n_bm
    for m=k:n_bm
        s_temp=0;
        for i=1:n_N
            s_temp=s_temp+((V_N(i,k)-Area(k))*(V_N(i,m)-Area(m)));
        end
        S01M(m,k)=s_temp/(n_N-1);
        S01M(k,m)=S01M(m,k);
    end
end
for k=1:n_bm
    for m=k:n_bm
        SM(m,k)=(S10M(m,k)/n_D)+(S01M(m,k)/n_N);
        SM(k,m)=SM(m,k);
    end
end

% summarise categories for output
%for k=1:n_bm
%    for h=1:n_N
%        s_N(k, Normals(k,h))=s_N(k, Normals(k,h))+1;
%    end
%    for h=1:n_D
%           s_D(k, Diseased(k,h))=s_D(k, Diseased(k,h))+1;
%    end
%end

%Calculations of chi2

L=[ 1, -1];
diff=L*Area';
temp2=Area*L';
LSLdash=L*SM*L';
x2=(diff/LSLdash)*temp2;
prob=1-chi2cdf(x2,1);

diff_plus95CI=diff+1.96*sqrt(LSLdash);
diff_minus95CI=diff-1.96*sqrt(LSLdash);
title={'AUC difference'};
AUC1string=  [num2str(Area(1),'%11.2g'),' (',num2str(CIMin95(1),'%11.2g'),' to ',num2str(CIPlus95(1),'%11.2g'),')'];
AUC2string=  [num2str(Area(2),'%11.2g'),' (',num2str(CIMin95(2),'%11.2g'),' to ',num2str(CIPlus95(2),'%11.2g'),')'];
AUCdiffstring=  [num2str(diff,'%11.2g'),' (',num2str(diff_minus95CI,'%11.2g'),' to ',num2str(diff_plus95CI,'%11.2g'),')'];
temp={title,AUC1string,AUC2string,AUCdiffstring,prob};
AUC_compare=temp;
end





