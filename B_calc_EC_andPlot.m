
% -------------------------------------------------------------------------


% Created by Trevor Keenan, UC Berkeley and LBNL (trevorkeenan@berkeley.edu)
% December 2021

% this script will...

% calulate the constrained Beta
% based on inputs of:
% the modeled relationship
% and the observational constraint

% THIS SCRIPT IS CALLED BY A_varianceNormalization.m
% NB: results may differ slightly from published values due to bootstrapping

% TFK, December 2021 (trevorkeenan@berkeley.edu)

% -------------------------------------------------------------------------



% load the modeled data and constraint
load('./dataIntermediates/emergentConstraint_data')

nboot       = dataEmergent.nboot;                % size of bootstrap sample;

nmodels     = length(dataEmergent.modeledDeltaGPP); % number of models

% scatter plot and regression line
modeled_normSland = dataEmergent.normalizedSland;

modeledDeltaGPP = dataEmergent.modeledDeltaGPP;
modeledDeltaReco = dataEmergent.modeledDeltaReco;

s = rng;
% create a distribution of observed Sland values from which to sample
% UNC-1: uncertainty in GCP Sland
observedNEE_distribution = normrnd(dataEmergent.constraint,dataEmergent.constraintSD,[1,200]);
nboot_obs   = size(observedNEE_distribution, 2);    % size of bootstrap sample of observations

%%
% regression estimate of ECS PDF by bootstrap
ECSboot = zeros(nboot, 1);
for j=1:nboot
    ib      = randsample(nmodels, nmodels, 1);
    % UNC-2: uncertainty in relationship between models and normalized Sland
    % UNC-3: uncertainty in nonResp flux
    [B, BINT, R, RINT, STATS] = regress(modeledDeltaGPP(ib), [ones(nmodels, 1) modeled_normSland(ib)]);
    sigma      = sqrt(STATS(4));
    
    % predicted ECS (with normal errors)
    % this bootstraps the observed Sland, and adds error associated with
    % the relationship between bootstrapped models and bootstrappedNEE
    ECSboot(j) = B(1)+B(2)*observedNEE_distribution(1, randsample(nboot_obs, 1)) + sigma*randn(1);
end
ECSboot    = sort(ECSboot);

%%
% some statistics
ECS_med = median(ECSboot);
ECS_mean = mean(ECSboot);
ECS_std = std(ECSboot);

conf_int   = .66;
ECS_SDu = ECSboot(round((1+conf_int)/2*nboot));
ECS_SDl = ECSboot(round((1-conf_int)/2*nboot));

disp('Emergent constraint: -1 SD, median , + 1 SD:')
disp([ECSboot(round((1-conf_int)/2*nboot)), ECS_med, ECSboot(round((1+conf_int)/2*nboot))])
disp('Emergent constraint: -1 SD, mean, + 1 SD:')
disp([ECSboot(round((1-conf_int)/2*nboot)), ECS_mean, ECSboot(round((1+conf_int)/2*nboot))])


%%
% generate PDF for Beta GPP constraint

fig1=figure;
hold on
bandwidth=0.01;
x2= 0:bandwidth:1; % the range of the x-axis

% plot the distribution of modeled GPP % change
y1=modeledDeltaGPP;
y1(y1==0)=NaN;

y2 = normpdf(x2,nanmean(y1),nanstd(y1));
h1=histogram(y1,8,'Normalization','pdf','Visible','off');
b1 = bar(h1.BinEdges(1:end-1)+h1.BinWidth/2, max(y2*bandwidth)*h1.Values/max(h1.Values));
set(b1,'FaceColor',[0.8 0.8 0.8])

% plot the unconstrained PDF
p1=plot(x2,y2*bandwidth,'k-','LineWidth',2);
h = area(x2,y2*bandwidth);
h(1).FaceColor = [0.4 0.4 0.4];
h(1).FaceAlpha = 0.4;

indX= y2==max(y2);
disp(strcat('mean model Beta:', num2str(nanmean(y1))))
disp(strcat('median model Beta:', num2str(nanmedian(y1))))
disp(strcat('mean model Beta (norm dist):', num2str(x2(indX))))
plot([nanmean(y1) nanmean(y1)],[0 y2(indX)*bandwidth],'k--','LineWidth',1)

% now add the distribution of constrained GPP % change
colorConstrained = [0.9290 0.5940 0.1250];

ECS_unc1 = (ECS_SDu-ECS_SDl)/2; % joint uncertainty
ECS_unc2 = 0.012;  % uncertainty from normalization bootstrapping
ECS_unc3 = 0.0025;  % uncertainty from nonResp Flux
ECS_unc = sqrt(ECS_unc1^2 + ECS_unc2^2 + ECS_unc3^2);

constrainedMeanStd= [ECS_mean,ECS_unc];
% save for later use
save('./dataIntermediates/constrainedRangeBetaGPP','constrainedMeanStd')

y2 = normpdf(x2,ECS_mean, ECS_unc);

u95 = ECS_mean + 1.96*ECS_unc;
l95 = ECS_mean - 1.96*ECS_unc;


p1=plot(x2,y2*bandwidth,'-','LineWidth',2);
set(p1,'color',colorConstrained)
disp(strcat('mean constrained Beta:', num2str(ECS_mean)))
set(p1,'color',colorConstrained)
h = area(x2,y2*bandwidth);
h(1).FaceColor = colorConstrained;
h(1).FaceAlpha = 0.4;

indX= y2==max(y2);
disp(strcat('mean constrained Beta (norm dist):', num2str(x2(indX))))
p1=plot([x2(indX) x2(indX)],[0 y2(indX)*bandwidth],'--','LineWidth',1);
set(p1,'color',colorConstrained)

tmpx = y1;
tmpx(tmpx==0)=[];
[~,~] = ksdensity(tmpx,x2,'Bandwidth',0.5);

text(0.275, 0.025,'TBM prior','color','k','FontSize',18)
text(0.58, 0.1,{'Constrained';' distribution'},'color',colorConstrained,'FontSize',18)

ylabel('Probabilty density')
xlabel('\beta_R^{GPP} (%_{GPP} %_{CO2}^{-1})')

xlim([0.25 0.9])
ylim([0 0.125])
set(gca,'FontSize',20)

disp('Constraint % diff from unconstrained')
disp(100*(ECS_mean-nanmean(y1))/nanmean(y1))
disp('note the mean here is not that reflected in fig. 1 as a model was removed')
disp(strcat('true change =', num2str(100*(1-ECS_mean/0.589)),'%'))

disp('Constraint % diff from max unconstrained')
disp(100*(ECS_mean-nanmax(y1))/nanmax(y1))

disp('Reduction in uncertainty')
% disp(100*((ECS_SDu-ECS_SDl)/2-nanstd(y1))/nanstd(y1)) % this is just the
% joint uncertainty
disp(100*(ECS_unc-nanstd(y1))/nanstd(y1)) % the full uncertainty (unc1+unc2+unc3)


disp('implied constrained GPP change from 1981-2020')
CO2_1982 = 341.57;
CO2_2020 = 416;
constraintedBeta= ECS_mean;
constrainedUncertainty  = ECS_unc; %(ECS_SDu-ECS_SDl)/2; % this is one std.
constrainedUncertainty95  = 1.96*ECS_unc; %(ECS_SDu-ECS_SDl)/2; % this is one std.

impliedDeltaPgC = 118*((CO2_2020-CO2_1982)/CO2_1982)*constraintedBeta;
impliedDeltaPgC_std = 118*((CO2_2020-CO2_1982)/CO2_1982)*constrainedUncertainty;
impliedDeltaPgC_95 = 118*((CO2_2020-CO2_1982)/CO2_1982)*constrainedUncertainty95;

impliedDeltaPgC_stdUnconstrained = 118*((CO2_2020-CO2_1982)/CO2_1982)*nanstd(y1);
impliedDeltaPgC_stdUnconstrained_pct = 100*impliedDeltaPgC_stdUnconstrained;

disp(strcat('PgC increase to 2020: ', num2str(impliedDeltaPgC)))
disp(strcat('                  +/- (std): ', num2str(impliedDeltaPgC_std),' PgC'))
disp(strcat('                  +/- (95 CI): ', num2str(impliedDeltaPgC_95),' PgC'))
disp(strcat('= % increase to 2020: ', num2str(100*impliedDeltaPgC/118)))
disp(strcat('                  +/- (std): ', num2str(100*impliedDeltaPgC_std/118)))
disp(strcat('                  +/- (95 CI): ', num2str(100*impliedDeltaPgC_95/118)))


if saveFigures==1
    fname=strcat('./figures/emergent/Fig1c');
    eval(printStatement)
end


%%
% generate PDF for Beta GPP constraint, with each uncertainty contribution

fig1=figure;
hold on
bandwidth=0.01;
x2= 0:bandwidth:1; % the range of the x-axis

% plot the distribution of modeled GPP % change
y1=modeledDeltaGPP;
y1(y1==0)=NaN;

y2 = normpdf(x2,nanmean(y1),nanstd(y1));
h1=histogram(y1,8,'Normalization','pdf','Visible','off');
b1 = bar(h1.BinEdges(1:end-1)+h1.BinWidth/2, max(y2*bandwidth)*h1.Values/max(h1.Values));
set(b1,'FaceColor',[0.8 0.8 0.8])

% plot the unconstrained PDF
p1=plot(x2,y2*bandwidth,'k-','LineWidth',2);
h = area(x2,y2*bandwidth);
h(1).FaceColor = [0.4 0.4 0.4];
h(1).FaceAlpha = 0.4;

indX= y2==max(y2);
plot([nanmean(y1) nanmean(y1)],[0 y2(indX)*bandwidth],'k--','LineWidth',1)

% now add the distribution of constrained GPP % change
colorConstrained = [0.9290 0.5940 0.1250];

ECS_unc1 = (ECS_SDu-ECS_SDl)/2; % joint uncertainty
ECS_unc2 = 0.012;  % uncertainty from normalization bootstrapping
ECS_unc3 = 0.0025;  % uncertainty from nonResp Flux
ECS_unc = sqrt(ECS_unc1^2 + ECS_unc2^2 + ECS_unc3^2);

y2 = normpdf(x2,ECS_mean, ECS_unc);

p1=plot(x2,y2*bandwidth,'-','LineWidth',2);
% add dashed line for ECS_unc1
ECS_unc1_cont = ECS_unc*(ECS_unc1/(ECS_unc1+ECS_unc2+ECS_unc3));
y2a = normpdf(x2,ECS_mean, ECS_unc1_cont);
ydatax = y2a*bandwidth;
ydatax = ydatax*(max(y2*bandwidth)/max(ydatax));
px1=plot(x2,ydatax,'k--','LineWidth',2);

% add dashed line for ECS_unc2
ECS_unc2_cont =ECS_unc1_cont + ECS_unc*(ECS_unc2/(ECS_unc1+ECS_unc2+ECS_unc3));

y2a = normpdf(x2,ECS_mean, ECS_unc2_cont);
ydatax = y2a*bandwidth;
ydatax = ydatax*(max(y2*bandwidth)/max(ydatax));
px2=plot(x2,ydatax,'--','LineWidth',2,'color',[0.5 0.5 0.5]);

% add dashed line for ECS_unc3
ECS_unc3_cont =ECS_unc2_cont + ECS_unc*(ECS_unc3/(ECS_unc1+ECS_unc2+ECS_unc3));
y2a = normpdf(x2,ECS_mean, ECS_unc3_cont);
ydatax = y2a*bandwidth;
ydatax = ydatax*(max(y2*bandwidth)/max(ydatax));
px3=plot(x2,ydatax,'k:','LineWidth',2);

set(p1,'color',colorConstrained)
h = area(x2,y2*bandwidth);
h(1).FaceColor = colorConstrained;
h(1).FaceAlpha = 0.4;

indX= y2==max(y2);
p1=plot([x2(indX) x2(indX)],[0 y2(indX)*bandwidth],'--','LineWidth',1);
set(p1,'color',colorConstrained)

tmpx = y1;
tmpx(tmpx==0)=[];
[~,~] = ksdensity(tmpx,x2,'Bandwidth',0.5);

text(0.275, 0.025,'TBM prior','color','k','FontSize',18)
text(0.58, 0.1,{'Constrained';' distribution'},'color',colorConstrained,'FontSize',18)

ylabel('Probabilty density')
xlabel('\beta_R^{GPP} (%_{GPP} %_{CO2}^{-1})')

xlim([0.25 0.9])
ylim([0 0.125])

set(gca,'FontSize',20)

legend([px1 px2 px3],{'Joint uncertainty','Normalization uncertainty','\gamma uncertainty'},...
    'box','off','location','east','FontSize',14)

if saveFigures==1
    fname=strcat('./figures/emergent/ED_Fig6');
    eval(printStatement)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Do the same for RECO, using BetaGPP as the constraint

inferredBetaGPPdistribution = normrnd(ECS_mean,ECS_unc,[1,200]);
nboot_obs   = size(inferredBetaGPPdistribution, 2);    % size of bootstrap sample of observations

% regression estimate of ECS PDF by bootstrap
ECSboot_reco = zeros(nboot, 1);
for j=1:nboot
    ib      = randsample(nmodels, nmodels, 1);
    [B, BINT, R, RINT, STATS] = regress(modeledDeltaReco(ib), [ones(nmodels, 1) modeledDeltaGPP(ib)]);
    sigma      = sqrt(STATS(4));
    
    % predicted ECS (with normal errors)
    ECSboot_reco(j) = B(1)+B(2)*inferredBetaGPPdistribution(1, randsample(nboot_obs, 1)) + sigma*randn(1);
end
ECSboot_reco    = sort(ECSboot_reco);


%%
% some statistics
ECS_med_reco    = median(ECSboot_reco);
ECS_mean_reco    = mean(ECSboot_reco);
ECS_std_reco = std(ECSboot_reco);

conf_int   = .66;
ECS_SDu_reco = ECSboot_reco(round((1+conf_int)/2*nboot));
ECS_SDl_reco = ECSboot_reco(round((1-conf_int)/2*nboot));

disp('RECO Emergent constraint: -1 SD, median , + 1 SD:')
disp([ECSboot_reco(round((1-conf_int)/2*nboot)), ECS_med_reco, ECSboot_reco(round((1+conf_int)/2*nboot))])
disp('RECO Emergent constraint: -1 SD, mean, + 1 SD:')
disp([ECSboot_reco(round((1-conf_int)/2*nboot)), ECS_mean_reco, ECSboot_reco(round((1+conf_int)/2*nboot))])


%%
% generate PDF for Beta Reco constraint

fig1=figure;
hold on
bandwidth=0.01;
x2= 0:bandwidth:1; % the range of the x-axis

% plot the distribution of modeled RECO beta
y1=modeledDeltaReco;
y1(y1==0)=NaN;

% plot the bars
y2 = normpdf(x2,nanmean(y1),nanstd(y1));
h1=histogram(y1,10,'Normalization','pdf','Visible','off');
b1 = bar(h1.BinEdges(1:end-1)+h1.BinWidth/2, max(y2*bandwidth)*h1.Values/max(h1.Values));
set(b1,'FaceColor',[0.8 0.8 0.8])

% plot the distribution of unconstrained RECO beta
plot(x2,y2*bandwidth,'k-','LineWidth',2)
h = area(x2,y2*bandwidth);
h(1).FaceColor = [0.4 0.4 0.4];
h(1).FaceAlpha = 0.4;

indX= y2==max(y2);
disp(strcat('mean model Beta:', num2str(nanmean(y1))))
disp(strcat('median model Beta:', num2str(nanmedian(y1))))
disp(strcat('mean model Beta (norm dist):', num2str(x2(indX))))
plot([nanmean(y1) nanmean(y1)],[0 y2(indX)*bandwidth],'k--','LineWidth',1)

% now add the distribution of constrained RECO beta
colorConstrained = [0.9290 0.5940 0.1250];
unc1 = (ECS_SDu_reco-ECS_SDl_reco)/2;
unc2 = ECS_std_reco;
y2 = normpdf(x2,ECS_mean_reco, unc1);
p1=plot(x2,y2*bandwidth,'-','LineWidth',2);
set(p1,'color',colorConstrained)
h = area(x2,y2*bandwidth);
h(1).FaceColor = colorConstrained;
h(1).FaceAlpha = 0.4;

disp(strcat('mean constrained Beta:', num2str(ECS_mean_reco)))
indX= y2==max(y2);
disp(strcat('mean constrained Beta (norm dist):', num2str(x2(indX))))
p1=plot([x2(indX) x2(indX)],[0 y2(indX)*bandwidth],'--','LineWidth',1);
set(p1,'color',colorConstrained)

text(0.275, 0.025,'TBM prior','color','k','FontSize',18)
text(0.525, 0.08,{'Constrained';' distribution'},'color',colorConstrained,'FontSize',18)

tmpx = y1;
tmpx(tmpx==0)=[];
[f,xi] = ksdensity(tmpx,x2,'Bandwidth',0.5);

ylabel('Probabilty density')
xlabel('\beta_R^{Reco} (%_{Reco} %_{CO2}^{-1})')

xlim([0.25 0.9])
ylim([0 0.125])

set(gca,'FontSize',20)

disp('Constraint % diff from unconstrained')
disp(100*(ECS_mean_reco-nanmean(y1))/nanmean(y1))
disp('Constraint % diff from max unconstrained')
disp(100*(ECS_mean_reco-nanmax(y1))/nanmax(y1))

disp('Reduction in uncertainty')
disp(100*((ECS_SDu_reco-ECS_SDl_reco)/2-nanstd(y1))/nanstd(y1))


if saveFigures==1
    fname=strcat('./figures/emergent/Fig1d');
    eval(printStatement)
end


