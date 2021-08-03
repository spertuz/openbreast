function r = rscore(impath, model, dflag)
% Compute risk score of mammogram
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
tic
%Feature extraction
% parameters:
res     = 0.1; %resolution 0.1mm/pixel
flist = {'imin', 'imax', 'iavg', 'ient', 'istd', 'ip05', 'ip95', 'iba1', 'iba2', 'ip30', 'ip70', 'iske', 'ikur','iran',...
'cene', 'ccor', 'ccon', 'chom', 'cent',...
'rsre', 'rlre', 'rgln', 'rrpe', 'rrln', 'rlgr', 'rhgr',...
'sgra', 'slap', 'swas', 'swav', 'swar', 'stev', 'fdim'};

info = getinfo(impath);
[imn, im] = ffdmRead(impath, info);
mask  = segBreast(imresize(imn, .25), info.ismlo);
im = imresize(im, info.psize/res);
f = xfeatures(im, flist, mask, res);
        
%retrieve sensor info:
if strcmpi(info.target, 'MO')&& strcmpi(info.filter, 'RH')
    sensor = 0;
elseif strcmpi(info.target, 'TU')&& strcmpi(info.filter, 'AL')
    sensor = 1;
elseif strcmpi(info.target, 'RH')&& strcmpi(info.filter, 'RH')
    sensor = 2;
elseif strcmpi(info.target, 'RH')&& strcmpi(info.filter, 'SI')
    sensor = 3;
else
    error('target-filter not found!')
end
        
%add imaging parameters:
x = [f, info.KVP, info.H, info.cforce, sensor];

%Predict risk score
r = predict(model.lr, x);
t = toc;
[~, fname] = fileparts(impath);
fprintf('File  : %s\n', fname)
fprintf('Score : %1.3f\n', r)
fprintf('Age   : %s\n', info.age)
fprintf('Time  : %1.3f\n', t)

%Display risk score:
if dflag&&length(r)==1
    bandwidth   = 0.05;
    edgecolor   = [.5 .5 .5];
    hi_color    = [1 0 0];
    lo_color    = [0 1 0];
    alpha       = 0.3;
    
    figure, 
    subplot(121)
    im = mat2gray(imresize(im, .25));
    imshow(im);
    
    str{1} = sprintf('View    : %s%s',info.side, info.view);
    str{2} = sprintf('Age     : %s', info.age);
    str{3} = sprintf('R score : %1.2f', r);
    str{4} = 'Status: HIGH risk';
    text(30,size(im,1)-30,str(:), 'color','g')
    
    subplot(122)
    
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