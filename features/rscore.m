function r = rscore(x, model, dflag)
% Compute risk score from LR model
% Sintax:
%     r = rscore(x, model)
%     r = rscore(x, model, dflag)
% Inputs:
%     x,      NxM feature matrix. Each row is a 1xM feature vector
%             corresponding to one image.
%     model,  structure risk model (as returned by rmodel).
%     dflag,  Optional boolean flag for displaying risk score, distribution
%             of risk among classes and confidence intervals.
%             Default is false.
%     
% Outputs,
%     r,      a scalar in [0,1] with the risk score. The higher
%             the better.
%             
% S. Pertuz
% Feb 05/2018
    
if nargin<3
    dflag = false;
end

x = (x-model.mu)./model.std;
r = glmval(model.weights, x, 'logit');

%Display risk score:
if dflag&&length(r)==1
    bandwidth   = 0.05;
    edgecolor   = [.5 .5 .5];
    hi_color    = [1 0 0];
    lo_color    = [0 1 0];
    alpha       = 0.3;
    
   
   %estimate density of low-risk:
    x0 = model.scores(~model.class);
    [density, value] = ksdensity(x0, 'bandwidth', bandwidth);
    density = density(value>=0&value<=1);
    density = 0.3*density/max(density);
    value = value(value>=0&value<=1);
    value(1) = 0;
    value(end) = 1;
    
    %plot density:
    xcor = [value, value(end:-1:1)];
    ycor = [density, -density(end:-1:1)];
    fill(xcor, ycor, [1 1 1], ...
        'facecolor', lo_color,...
        'edgecolor', edgecolor,...
        'facealpha', alpha);
    hold on
    
     %estimate density of high-risk:
    x1 = model.scores(model.class);
    [density, value] = ksdensity(x1, 'bandwidth', bandwidth);
    density = density(value>=0&value<=1);
    density = 0.3*density/max(density);
    value = value(value>=0&value<=1);
    value(1) = 0;
    value(end) = 1;
    
    %plot density:
    xcor = [value, value(end:-1:1)];
    ycor = [density, -density(end:-1:1)];
    fill(xcor, ycor, [1 1 1], ...
        'facecolor', hi_color,...
        'edgecolor', edgecolor,...
        'facealpha', alpha);
    
    legend({'High risk', 'Low risk'},...
        'Location','southwest',...
        'orientation','horizontal',...
        'autoupdate', 'off');
    
    
    %show mini box-plot
    Q = quantile(x1, [.25 .50 .75]);
    IQR = Q(3)-Q(1);
    fill(Q([1 1 3 3]), 0.01*[-1 1 1 -1], [1 1 1],...
        'facecolor', edgecolor, 'edgecolor', edgecolor);
    
    %show median:
    scatter(Q(2), 0, [], [1 1 1], 'filled')   
    
    %show mini box-plot
    Q = quantile(x0, [.25 .50 .75]);
    IQR = Q(3)-Q(1);
    fill(Q([1 1 3 3]), 0.01*[-1 1 1 -1], [1 1 1],...
        'facecolor', edgecolor, 'edgecolor', edgecolor);
    
    %show median:
    scatter(Q(2), 0, [], [1 1 1], 'filled')
    
%     %Notches for confidence intervals:
%     IQR = Q(3)-Q(1);
%     L_hi = Q(2) + 1.57*IQR/sqrt(length(x));
%     L_lo = Q(2) - 1.57*IQR/sqrt(length(x));
%     plot([L_lo, L_hi], [0 0], 'd')
    
%     %Show whiskers:
%     lo_wisk = max([0, Q(1)-1.5*IQR]);
%     hi_wisk = min([1, Q(3)+1.5*IQR]);
%     plot([Q(3), hi_wisk], [0 0], 'k', 'linewidth', 2)
%     plot([lo_wisk, Q(1)], [0 0], 'k', 'linewidth', 2)
    
    %show risk line:
    line([r, r], [-.4 .4], 'linewidth', 2, 'color', edgecolor)
    scatter(r, -0.4, [], edgecolor, 'filled', '^')
    scatter(r, +0.4, [], edgecolor, 'filled', 'v')
    
    str = sprintf('%1.2f', r);
    text(r-.05, -.42, str, 'fontsize', 12)
    set(gca,'xtick',[0 1], 'ytick',[])
    xlabel('Risk score')
end