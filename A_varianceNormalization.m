
% -------------------------------------------------------------------------


% Created by Trevor Keenan, UC Berkeley and LBNL (trevorkeenan@berkeley.edu)
% December 2021

% this script will...

% 1. load and plot trendy model quantities
% 2. use variance normalization to extract the partial relationship
%        between Beta_GPP and Sland
% 2. derive the emergent constraint (by calling B_calc_EC_andPlot.m)

% Figures produced:
% Figure 1a-d
% ED Figure 1
% ED Figure 2
% ED Figure 3
% ED Figure 6

% runtime ~ 40s, MBP M1 Max, 2021, Matlab R2021a
% NB: results may differ slightly from published values 
%       due to bootstrapping

% -------------------------------------------------------------------------


clearvars
close all

saveFigures =0;
addpath('./functions')
addpath('./functions/press')

printStatement='print(fig1,''-dpdf'', fname)';

nboot = 100000;
alpha=0.05;  % define the confidence interval for the estimation of beta (0.1 = 90%)
% later converted and presented as standard deviations

% 1. load the model derived data
mat = dir('./TRENDYv6_derived/*.mat');
for q = 1:length(mat)
    load(strcat('./TRENDYv6_derived/', mat(q).name));
end
clear mat q

%% 1. plot original Beta-GPP vs the Sink
fig1=figure;

hold on
for ii=1:width(T_R_GPP_S1)
    cModel = T_R_GPP_S1.Properties.VariableNames(ii);
    % get the Sland value
    xdata(ii) = T_normalizedSland.(cModel{:})(4);
    % get the Response Ratio (beta) of GPP
    ydata(ii) = T_R_GPP_S1.(cModel{:})(1);
    if isnan(xdata(ii))
        disp(strcat('Missing: ',cModel,': from normalized Sland'))
    end
end

xdataOrigSland = xdata;
plot(xdata,ydata,'.','MarkerSize',10)
[r p] =corrcoef(xdata,ydata,'rows','complete');
xlabel('\it{S}\rm_{LAND} GCP (PgC yr^{-1})')
ylabel('\beta_R^{GPP} (%_{GPP} %_{CO2}^{-1})')
set(gca,'FontSize',18)

% get the linear regression
x = xdata;
y = ydata;
indX=isnan(x);
x(indX)=[];y(indX)=[];
[poly, S] = polyfit(x,y,1);
slightExt = 0.05;
x2 = (min(x)-slightExt):0.01:(max(x)+slightExt);
xfit = x2;
alpha1 = 0.05;	% 68% CI Significance level
[Y,DELTA] = polyconf(poly,xfit,S,'alpha',alpha,'simopt','off','predopt','curve');
% convert CI to SE
DELTA = DELTA*2/3.92;
% plot the regression to the TRENDY models
hconf = plot(xfit,Y+DELTA,'r--');
plot(xfit,Y-DELTA,'r--')
shadedplot(xfit, Y+DELTA, Y-DELTA, [1 0.7 0.7],'none');
yfit = polyval(poly,xfit);
p(2)=plot(x2,yfit,'r-','LineWidth',1);

plot(x,y,'k.','MarkerSize',30);

ylim([0.3 0.9])

z_names = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N'};
nudgeX = 0.05*ones(size(xdata));
nudgeY = 0.01*ones(size(ydata));
xdata2=xdata+nudgeX;
ydata2=ydata+nudgeY;
text(xdata2,ydata2,z_names,'FontSize',20)

[r,p] =corrcoef(xdata,ydata,'rows','complete');
r_text1 = strcat('r^2= ',{' '},num2str(round(r(1,2)^2,2)));
r_text2 = strcat('p= ',{' '},num2str(round(p(1,2),2)));
r_text = {r_text1{:};r_text2{:}};
text(3.25,0.38,r_text,'FontSize',20)

if saveFigures==1
    fname='./figures/emergent/ED_Fig1';
    eval(printStatement)
end

%%
% Load the Global Carbon Project 2017 data (consistent with TRENDYv6)
GCPbudget = readtable('../data/GCP/GlobalCarbonBudgetSheet.xls');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create a table with the essentials
clear tmp
for ii = 1:width(T_R_GPP_S1)
    cModel = T_R_GPP_S1.Properties.VariableNames{ii};
    modelNames{ii}=cModel;
    
    tmp(1,ii) = T_normalizedSland.(cModel)(4); % b: Sland;
    
    tmp(2,ii) = T_R_GPP_S1.(cModel)(1); % c: Beta_GPP
    tmp(3,ii) = T_NonRespFlux.(cModel)(1);  % d: nonRespFlux
    try
        tmp(4,ii) = T_R_RECO_S1.(cModel)(1); % e: Beta_Reco
    catch
        tmp(4,ii) = NaN;
    end
    tmp(5,ii) = tmp(4,ii);    % f: Beta_Reco, not normalized
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% model Sland as a function of predictors

% generate the table for the linear model
T = array2table(tmp',...
    'RowNames',modelNames,...
    'VariableNames',{'b','c','d','e','f'});

% make a copy of T for manipulating
Tsub2 = T(:,{'b','c','d','e','f'});

indX=isnan(Tsub2.d);
Tsub2.d(indX)=nanmean(Tsub2.d);

% normalize d, e to zero mean
Tsub3=Tsub2;
Tsub3.d=(Tsub3.d-nanmean(Tsub3.d));
Tsub3.e=(Tsub3.e-nanmean(Tsub3.e));

% fit models
lmTest4 = fitlm(Tsub3,'b~c+f+d:e')
X=[Tsub3.c,Tsub3.f,Tsub3.d.*Tsub3.e];
y = Tsub3.b;
indX = isnan(y); y(indX)=[];X(indX,:)=[];
indX = isnan(X(:,3));y(indX)=[];X(indX,:)=[];
pr = press([X y]);
rPred = 1-pr/lmTest4.SST;
disp(strcat('Predicted r2:' ,num2str(rPred)))


%%
ypredOpt = predict(lmTest4,Tsub3);
fig1=figure;
hold on
xdata=Tsub3.b;
ydata =ypredOpt;

% get the linear regression
x = xdata;
y = ydata;
indX=isnan(x);
x(indX)=[];y(indX)=[];
indX=isnan(y);
x(indX)=[];y(indX)=[];

[poly, S] = polyfit(x,y,1);
slightExt = 0.05;
x2 = (min(x)-slightExt):0.01:(max(x)+slightExt);
xfit = x2;
alpha = 0.05;	% 95% CI Significanc e level
[Y,DELTA] = polyconf(poly,xfit,S,'alpha',alpha,'simopt','off','predopt','curve');

% plot the regression to the TRENDY models
hconf = plot(xfit,Y+DELTA,'r--');
plot(xfit,Y-DELTA,'r--')
shadedplot(xfit, Y+DELTA, Y-DELTA, [1 0.7 0.7],'none');
yfit = polyval(poly,xfit);
p(2)=plot(x2,yfit,'r-','LineWidth',1);

plot(xdata,ydata,'k.','MarkerSize',20);

z_names = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N'};
nudgeX = 0.05*ones(size(xdata));
nudgeX(strcmp(z_names,'C'))=nudgeX(strcmp(z_names,'C'))-0.2;
nudgeY = 0.01*ones(size(ydata));
nudgeY(strcmp(z_names,'K'))=nudgeY(strcmp(z_names,'K'))+0.1;
xdata2=xdata+nudgeX;
ydata2=ydata+nudgeY;
text(xdata2,ydata2,z_names,'FontSize',20)

[r,p] =corrcoef(xdata,ydata,'rows','complete');
r_text1 = strcat('r^2= ',{' '},num2str(round(r(1,2)^2,2)));
r_text2 = strcat('p < 0.01');
r_text = {r_text1{:};r_text2};
text(3.25,1.5,r_text,'FontSize',20)

xlabel('\it{S}\rm_{LAND} GCP (PgC yr^{-1})')
ylabel('Estimated \it{S}\rm_{LAND} (PgC yr^{-1})')
[r, ~] = corrcoef(Tsub3.b, ypredOpt,'rows','complete');

set(gca,'FontSize',18,'box','on')

xlim([0.75 4.25])
ylim([0.75 4.25])

if saveFigures==1
    fname='./figures/emergent/ED_Fig2';
    eval(printStatement)
end


%%
intercept = lmTest4.Coefficients.Estimate(1);
paramC= lmTest4.Coefficients.Estimate(2);
paramF= lmTest4.Coefficients.Estimate(3);
paramDE= lmTest4.Coefficients.Estimate(4);

predAll = intercept + paramC*Tsub3.c + paramF*Tsub3.f ...
    + paramDE*(Tsub3.d.*Tsub3.e);

% here normalizing by betaReco and interaction
predDE = paramF*Tsub3.f+...
    + paramDE*(Tsub3.d.*Tsub3.e);

normalizedSland = Tsub3.b - (predDE-nanmean(predDE)); % Sland - variance factors

% GCP uncertainty of decadal flux for past 6 decades 
% (0.9 PgC; Le Quere et al. 2017, Table 5 residual sink 1 std)
y = GCPbudget.ResidualSink';
indX1 = find(GCPbudget.Year==1982);
indX2 = find(GCPbudget.Year==2012);
stderror = 0.9/sqrt(6);
GCPbudgetBaseline(1,1:2) = [mean(y(indX1:indX2)), stderror];

% create a bootstrapped normalizedSland, sampling across values of
% nonRespFlux (j)
% this is used to propogate the uncertainty from nonRespFlux through to the final
% constraint
normalizedSland_NonRespBoot=nan(size(Tsub3.b,1),nboot);
for ii=1:nboot
    r = normrnd(0,std(Tsub3.d));
    r_uniform = rand*(nanmean(T.d));
    
    % assuming that the model ensemble is distributed around the mean,
    % propogate the associated uncertainty
    predDE_r = paramF*Tsub3.f + paramDE*((r+Tsub3.d).*Tsub3.e);
    normalizedSland_NonRespBoot(:,ii) =  Tsub3.b - (predDE_r-nanmean(predDE_r));
    
    % now get the uncertainty of the derived beta
    xdata = normalizedSland_NonRespBoot(:,ii);
    ydata = Tsub3.c;
    
    % get the linear regression
    x = xdata;
    y = ydata;
    indX=isnan(x);
    x(indX)=[];y(indX)=[];
    [poly, ~] = polyfit(x,y,1);
    
    % get the value of the fit line at x=GCP Sland
    emergentBeta_bootNonResp(ii)=polyval(poly,GCPbudgetBaseline(1,1));
    
end


%%
fig1=figure;
hold on;
shadedColor = [0.9 0.9 0.9];

xdata = normalizedSland;
ydata = Tsub3.c;

plot(xdata,ydata,'.','MarkerSize',20)

% get the linear regression
x = xdata;
y = ydata;
indX=isnan(x);
x(indX)=[];y(indX)=[];
[poly, S] = polyfit(x,y,1);
slightExt = 0.05;
x2 = (min(x)-slightExt):0.01:(max(x)+slightExt);
xfit = x2;
alpha = 0.05;	% 95% CI Significance level
[Y,DELTA] = polyconf(poly,xfit,S,'alpha',alpha,'simopt','off','predopt','curve');
% plot the regression to the TRENDY models
hconf = plot(xfit,Y+DELTA,'r--');
plot(xfit,Y-DELTA,'r--')
shadedplot(xfit, Y+DELTA, Y-DELTA, [1 0.7 0.7],'none');
yfit = polyval(poly,xfit);
p(2)=plot(x2,yfit,'r-','LineWidth',1);

% add the VERTICAL shaded area
emergent = GCPbudgetBaseline(1,1);
stde = GCPbudgetBaseline(1,2);
shaded_x1 = emergent -stde;
shaded_x2= emergent + stde;
shaded_maxy = 1;%polyval(poly,emergent);
ha = area([shaded_x1 shaded_x2], [shaded_maxy shaded_maxy] ,'FaceColor',shadedColor,'LineStyle','none');
ha.FaceAlpha = 0.6;
% add the vertical line
y1a=1;
p(1)=plot([emergent emergent],[0 y1a],'k--');

% add the horizontal emergent constraint line
x1a=[0 emergent];
x1a=[0 7];
y1a=polyval(poly,emergent);
p(1)=plot(x1a,[y1a y1a],'k--');

z_names = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N'};
nudgeX = 0.1*ones(size(xdata));
nudgeX(strcmp(z_names,'G'))=nudgeX(strcmp(z_names,'G'))-0.1;
nudgeX(strcmp(z_names,'I'))=nudgeX(strcmp(z_names,'I'))-0.05;
nudgeX(strcmp(z_names,'N'))=nudgeX(strcmp(z_names,'N'))-0.01;
nudgeY = 0.01*ones(size(ydata));
nudgeY(strcmp(z_names,'G'))=nudgeY(strcmp(z_names,'G'))+0.02;
nudgeY(strcmp(z_names,'I'))=nudgeY(strcmp(z_names,'I'))+0.01;
nudgeY(strcmp(z_names,'N'))=nudgeY(strcmp(z_names,'N'))-0.02;
xdata2=xdata+nudgeX;
ydata2=ydata+nudgeY;
text(xdata2,ydata2,z_names,'FontSize',20)

plot(xdata,ydata,'k.','MarkerSize',30)

[r,p] =corrcoef(xdata,ydata,'rows','complete');
r_text1 = strcat('r^2= ',{' '},num2str(round(r(1,2)^2,2)));
r_text2 = strcat('p < 0.01');
r_text = {r_text1{:};r_text2};
text(5,0.4,r_text,'FontSize',18)


xlabel('Normalized Terrestrial Sink (PgC yr^{-1})')
ylabel('\beta_R^{GPP} (%_{GPP} %_{CO2}^{-1})')
set(gca,'FontSize',18)

ylim([0.3 0.9])
xlim([0.25 6.25])

if saveFigures==1
    fname='./figures/emergent/Fig1a';
    eval(printStatement)
end


%%
fig1=figure;
hold on;
shadedColor = [0.9 0.9 0.9];

xdata = Tsub3.f;
ydata = Tsub3.c;

% remove VEGAS
xdata(end-1)=NaN;
ydata(end-1)=NaN;

plot(xdata,ydata,'.','MarkerSize',20)

[r,m,b] = regression(xdata',ydata');

plot(xdata,b+m*xdata,'r-')
r1=plot([0.3,1],[0.3,1],'r--');
set(r1,'LineStyle','--','color','r')

% get the linear regression
x = xdata;
y = ydata;
indX=isnan(x);
x(indX)=[];y(indX)=[];
[poly, S] = polyfit(x,y,1);
slightExt = 0.05;
x2 = (min(x)-slightExt):0.01:(max(x)+slightExt);
xfit = x2;
alpha = 0.05;	% 95% CI Significance level
[Y,DELTA] = polyconf(poly,xfit,S,'alpha',alpha,'simopt','off','predopt','curve');

% plot the regression to the TRENDY models
hconf = plot(xfit,Y+DELTA,'r--');
plot(xfit,Y-DELTA,'r--')
shadedplot(xfit, Y+DELTA, Y-DELTA, [1 0.7 0.7],'none');
yfit = polyval(poly,xfit);
p(2)=plot(x2,yfit,'r-','LineWidth',1);

p(1)=plot([0.485 0.485],[0 1],'k--');

% add the HORIZONTAL shaded area
x1 = [0 shaded_x2];
% add the horizontal emergent constraint line
p(1)=plot(x1a,[y1a y1a],'k--');

z_names = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N'};
nudgeX = 0.01*ones(size(xdata));
nudgeX(strcmp(z_names,'G'))=nudgeX(strcmp(z_names,'G'))-0.05;
nudgeX(strcmp(z_names,'I'))=nudgeX(strcmp(z_names,'I'))-0.005;
nudgeX(strcmp(z_names,'N'))=nudgeX(strcmp(z_names,'N'))-0.04;
nudgeY = 0.01*ones(size(ydata));
nudgeY(strcmp(z_names,'I'))=nudgeY(strcmp(z_names,'I'))+0.01;
nudgeY(strcmp(z_names,'N'))=nudgeY(strcmp(z_names,'N'))-0.02;
xdata2=xdata+nudgeX;
ydata2=ydata+nudgeY;
text(xdata2,ydata2,z_names,'FontSize',20)

plot(xdata,ydata,'k.','MarkerSize',30)

[r,p] =corrcoef(xdata,ydata,'rows','complete');
r_text1 = strcat('r^2= ',{' '},num2str(round(r(1,2)^2,2)));
r_text2 = strcat('p < 0.01');
r_text = {r_text1{:};r_text2};
text(0.775,0.4,r_text,'FontSize',18)


xlabel('\beta_R^{Reco} (%_{Reco} %_{CO2}^{-1})')
% ylabel('\beta_R^{GPP} (%_{GPP} %_{CO2}^{-1})')
set(gca,'FontSize',18)
set(gca,'YTickLabel','')

ylim([0.3 0.9])
xlim([0.3 0.9])

if saveFigures==1
    fname='./figures/emergent/Fig1b';
    eval(printStatement)
end

%%
fig1=figure;
hold on;
shadedColor = [0.9 0.9 0.9];

xdata = normalizedSland;
ydata = Tsub3.c;

plot(xdata,ydata,'.','MarkerSize',20)

% get the linear regression
x = xdata;
y = ydata;
indX=isnan(x);
x(indX)=[];y(indX)=[];
[poly, S] = polyfit(x,y,1);
slightExt = 0.05;
x2 = (min(x)-slightExt):0.01:(max(x)+slightExt);
xfit = x2;
alpha = 0.05;	% 95% CI Significance level
[Y,DELTA] = polyconf(poly,xfit,S,'alpha',alpha,'simopt','off','predopt','curve');
% plot the regression to the TRENDY models
hconf = plot(xfit,Y+DELTA,'r--');
plot(xfit,Y-DELTA,'r--')
shadedplot(xfit, Y+DELTA, Y-DELTA, [1 0.7 0.7],'none');
yfit = polyval(poly,xfit);
p(2)=plot(x2,yfit,'r-','LineWidth',1);

% add the VERTICAL shaded area
emergent = GCPbudgetBaseline(1,1);
stde = GCPbudgetBaseline(1,2);
shaded_x1 = emergent -stde;
shaded_x2= emergent + stde;
shaded_maxy = 1;
ha = area([shaded_x1 shaded_x2], [shaded_maxy shaded_maxy] ,'FaceColor',shadedColor,'LineStyle','none');
ha.FaceAlpha = 0.6;
% add the vertical line
y1a=1;
p(1)=plot([emergent emergent],[0 y1a],'k--');

% add the HORIZONTAL shaded area
x1 = [0 shaded_x2];
% add the horizontal emergent constraint line
x1a=[0 emergent];
y1a=polyval(poly,emergent);
p(1)=plot(x1a,[y1a y1a],'k--');

z_names = modelNames;
z_names{2} =strrep(z_names{2},'_','-');
z_names{8} =strrep(z_names{8},'_','-');
z_names{11} =strrep(z_names{11},'_','-');

nudgeX = 0.1*ones(size(xdata));
nudgeX(strcmp(z_names,'VISIT'))=nudgeX(strcmp(z_names,'VISIT'))-0.75;
nudgeX(strcmp(z_names,'JULES'))=nudgeX(strcmp(z_names,'JULES'))-0.5;
nudgeY = 0.01*ones(size(ydata));
nudgeY(strcmp(z_names,'LPJ'))=nudgeY(strcmp(z_names,'LPJ'))+0.03;
nudgeY(strcmp(z_names,'ORCHIDEE-MICT'))=nudgeY(strcmp(z_names,'ORCHIDEE-MICT'))-0.015;
nudgeY(strcmp(z_names,'JULES'))=nudgeY(strcmp(z_names,'JULES'))+0.02;
nudgeY(strcmp(z_names,'CLASS-CTEM'))=nudgeY(strcmp(z_names,'CLASS-CTEM'))+0.02;
xdata2=xdata+nudgeX;
ydata2=ydata+nudgeY;
text(xdata2,ydata2,z_names,'FontSize',20)

plot(xdata,ydata,'k.','MarkerSize',30);

xlabel('Normalized Terrestrial Sink (PgC yr^{-1})')
ylabel('\beta_R^{GPP} (%_{GPP} %_{CO2}^{-1})')
set(gca,'FontSize',18)

ylim([0.3 0.9])
xlim([0.25 6.25])

if saveFigures==1
    fname='./figures/emergent/ED_Fig3';
    eval(printStatement)
end


%% save a file with modeled data and emergent constraint
dataEmergent.modeledDeltaGPP = y;
x1 = Tsub3.f;
x1([8,13],:)=[];    % model 13 is empty for NEE
dataEmergent.modeledDeltaReco = x1;

dataEmergent.modeledNEE  = x;
x2=normalizedSland;
x2([8,13],:)=[];    % model 13 is empty
dataEmergent.normalizedSland = x2;

dataEmergent.constraint = GCPbudgetBaseline(1,1);
dataEmergent.constraintSD = GCPbudgetBaseline(1,2);

dataEmergent.nboot=nboot;

% save the intermediate file for the emergent constraint
save('./dataIntermediates/emergentConstraint_data','dataEmergent')

% run the emergent constraint code
run('./B_calc_EC_andPlot.m')
