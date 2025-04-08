%% coastline
[lat2,lon2,z] = read_kml('coast.kml');
% plot(lon2,lat2,'-k')


%% TG and SA load data

load('JA3_cor.mat')
load('S3A_cor.mat')
load('S3B_cor.mat')


load('Z:\data\processed_data\Tide Gauge Data\ProcessedTGdata.mat')
id=9;

[m1,~]=find(TGdate==dateshift(JA3_cor.Time(1),'start','hour','next')); % first date
[m2,~]=find(TGdate==dateshift(S3A_cor.Time(end),'start','hour','next')); % last date

lat_tg=TGinfo.Lat(id);
lon_tg=TGinfo.Lon(id);
dt_tg=TG(m1:m2,id);
date_tg=TGdate(m1:m2);

%% DTU MSS 18

latdtu18=[53.0000000000:0.0166666600:66.9999944000];
londtu18=[8.0000000000:0.0166666600:30.9999908000]';
fid=fopen('mss18.grd','r');
t=textscan(fid, '%f');
mssdtuu18=reshape(t{1,1},[1381,841]);  %cremove 6 first raw *******************


[lon_geo, lat_geo, z_geo]=grdread2('NKG2015.grd');
ndtu18=griddata(lon_geo, lat_geo, double(z_geo),londtu18,latdtu18); %m

mdtdtu18=(mssdtuu18'-ndtu18);

% mdt_tg=griddata(londtu18,latdtu18,mdtdtu18,lon_tg,lat_tg)*100; 

%% DTU MDT 19
latdtu19=[-80:0.12500000000000000:90];
londtu19=[0:0.12500000000000000:360];


fid=fopen('dtu19mdt.grd','r');
t=textscan(fid, '%f');
mdtdtu19=reshape(t{1,1},[2881,1361]); %m


mdt_tg=griddata(londtu19,latdtu19,mdtdtu19',lon_tg,lat_tg)*100;

%% DTU MDT 21  %dont run
latdtu21=ncread('DTU21MSS_1min.mss.nc','lat');
londtu21=ncread('DTU21MSS_1min.mss.nc','lon');
mssdtu21=ncread('DTU21MSS_1min.mss.nc','mss'); %m

mss_tg=griddata(londtu21,latdtu21,mssdtu21',lon_tg,lat_tg);

%% plot
figure(1)
contourf(londtu18,latdtu18,mssdtuu18','LineColor','none')
hold on
colorbar
[lat2,lon2,z] = read_kml('coast.kml');
plot(lon2,lat2,'-k')
caxis([14.1831 50.8624]);
ylim([53 67])
xlim([8 31])

figure(2)
contourf(lon_geo, lat_geo,z_geo,'LineColor','none')
hold on
colorbar
plot(lon2,lat2,'-k')
caxis([14.1831 50.8624]);
ylim([53 67])
xlim([8 31])

%MDT DTU18
figure(3)
contourf(londtu18,latdtu18,mdtdtu18*100,'LineColor','none')
hold on
colorbar
[lat2,lon2,z] = read_kml('coast.kml');
plot(lon2,lat2,'-k')
% caxis([14.1831 50.8624]);
ylim([53 67])
xlim([8 31])


%MDT DTU19
figure(4)
contourf(londtu19,latdtu19,mdtdtu19'*100,'LineColor','none')
hold on
colorbar
[lat2,lon2,z] = read_kml('coast.kml');
plot(lon2,lat2,'-k')
% caxis([14.1831 50.8624]);
ylim([53 67])
xlim([8 31])

%%  TG SLA and empty table SA SLA
% sla_tg=dt_tg-mdt_tg; % TG SLA from DTU model
sla_tg=dt_tg-mean(TG(:,id)); % TG SLA from TG data

sla = cell2table(cell(0,5));
sla.Properties.VariableNames{'Var1'}='date_sa';
sla.Properties.VariableNames{'Var2'}='sla_sa';
sla.Properties.VariableNames{'Var3'}='pass';
sla.Properties.VariableNames{'Var4'}='mission';
sla.Properties.VariableNames{'Var5'}='said';



clear m1
clear m2

%% 
%S3A
pass=283;
% pass=728;



sa=table(S3A_cor.id(S3A_cor.Pass==pass),S3A_cor.Time(S3A_cor.Pass==pass),S3A_cor.Lat(S3A_cor.Pass==pass),S3A_cor.Lon(S3A_cor.Pass==pass),S3A_cor.sacorr(S3A_cor.Pass==pass));
sa.Properties.VariableNames{1} = 'id';
sa.Properties.VariableNames{2} = 'Time';
sa.Properties.VariableNames{3} = 'Lat';
sa.Properties.VariableNames{4} = 'Lon';
sa.Properties.VariableNames{5} = 'DT';
sa= rmmissing(sa,'DataVariables',{'DT'}); %remove NaN data

for i=1:height(sa)
sa.dist(i)=distance([sa.Lat(i) sa.Lon(i)],[lat_tg,lon_tg],referenceEllipsoid('WGS84'))/1000; %measure TG to SA points disance
end

sa=sa(sa.dist<=30,:);

ID=unique(sa.id);
for i=1:length(ID)
mdt_sa(i,1)=mean((sa.DT(sa.id==ID(i))));
t=sa.Time(sa.id==ID(i));
date_sa(i,1)=dateshift(t(1),'start','hour','next');
clear t
end
sla_sa=mdt_sa-mean(TG(:,id));

clear ID
clear sa



sla1=table(date_sa,sla_sa);
sla1.pass(:,1)=pass;
mission = 'S3A';
sla1.mission(:,1) ={mission};
sla1.said(:,1) =2;



sla=[sla;sla1];
clear sla1
clear date_sa
clear sla_sa
clear mdt_sa
clear sa

%% %% 
%JA3
% pass=111;
pass=194;



sa=table(JA3_cor.id(JA3_cor.Pass==pass),JA3_cor.Time(JA3_cor.Pass==pass),JA3_cor.Lat(JA3_cor.Pass==pass),JA3_cor.Lon(JA3_cor.Pass==pass),JA3_cor.sacorr(JA3_cor.Pass==pass));
sa.Properties.VariableNames{1} = 'id';
sa.Properties.VariableNames{2} = 'Time';
sa.Properties.VariableNames{3} = 'Lat';
sa.Properties.VariableNames{4} = 'Lon';
sa.Properties.VariableNames{5} = 'DT';
sa= rmmissing(sa,'DataVariables',{'DT'}); %remove NaN data

for i=1:height(sa)
sa.dist(i)=distance([sa.Lat(i) sa.Lon(i)],[lat_tg,lon_tg],referenceEllipsoid('WGS84'))/1000; %measure TG to SA points disance
end

sa=sa(sa.dist<=30,:);

ID=unique(sa.id);
for i=1:length(ID)
mdt_sa(i,1)=mean((sa.DT(sa.id==ID(i))));
t=sa.Time(sa.id==ID(i));
date_sa(i,1)=dateshift(t(1),'start','hour','next');
clear t
end
sla_sa=mdt_sa-mean(TG(:,id));

clear ID
clear sa



sla1=table(date_sa,sla_sa);
sla1.pass(:,1)=pass;

mission = 'JA3';
sla1.mission(:,1) ={mission};
sla1.said(:,1) =1;

sla=[sla;sla1];
clear sla1
clear date_sa
clear sla_sa
clear mdt_sa
clear sa
%% %% 
%S3B
pass=283;
% pass=728;



sa=table(S3B_cor.id(S3B_cor.Pass==pass),S3B_cor.Time(S3B_cor.Pass==pass),S3B_cor.Lat(S3B_cor.Pass==pass),S3B_cor.Lon(S3B_cor.Pass==pass),S3B_cor.sacorr(S3B_cor.Pass==pass));
sa.Properties.VariableNames{1} = 'id';
sa.Properties.VariableNames{2} = 'Time';
sa.Properties.VariableNames{3} = 'Lat';
sa.Properties.VariableNames{4} = 'Lon';
sa.Properties.VariableNames{5} = 'DT';
sa= rmmissing(sa,'DataVariables',{'DT'}); %remove NaN data

for i=1:height(sa)
sa.dist(i)=distance([sa.Lat(i) sa.Lon(i)],[lat_tg,lon_tg],referenceEllipsoid('WGS84'))/1000; %measure TG to SA points disance
end

sa=sa(sa.dist<=30,:);

ID=unique(sa.id);
for i=1:length(ID)
mdt_sa(i,1)=mean((sa.DT(sa.id==ID(i))));
t=sa.Time(sa.id==ID(i));
date_sa(i,1)=dateshift(t(1),'start','hour','next');
clear t
end
sla_sa=mdt_sa-mean(TG(:,id));

clear ID
clear sa
clear mdt_sa

% 
sla1=table(date_sa,sla_sa);
sla1.pass(:,1)=pass;
mission = 'S3B';
sla1.mission(:,1) ={mission};
sla1.said(:,1) =3;


sla=[sla;sla1];
clear sla1
clear date_sa
clear sla_sa
clear mdt_sa
clear sa

%% 

sla = sortrows(sla,'date_sa','ascend');

%% plot
tg=table(date_tg,sla_tg);

plot(tg.date_tg,tg.sla_tg,'k')
hold on

plot(sla.date_sa(sla.said==3),sla.sla_sa(sla.said==3),'o','color',[0.4660 0.6740 0.1880],'MarkerSize',7,'MarkerFaceColor',[0.4660 0.6740 0.1880]) %S3B
plot(sla.date_sa(sla.said==2),sla.sla_sa(sla.said==2),'o','color',[0 0.4470 0.7410],'MarkerSize',7,'MarkerFaceColor',[0 0.4470 0.7410]) %S3A
plot(sla.date_sa(sla.said==1),sla.sla_sa(sla.said==1),'or','MarkerSize',7,'MarkerFaceColor','r') % JA3


ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);


%% Geoplot

geoplot(lat_tg,lon_tg,'^k','MarkerFaceColor','k','MarkerSize',10)
hold on
%S3B
geoplot(S3B_cor.Lat(S3B_cor.id==228),S3B_cor.Lon(S3B_cor.id==228),'-','color',[0.4660 0.6740 0.1880]) %728
geoplot(S3B_cor.Lat(S3B_cor.id==153),S3B_cor.Lon(S3B_cor.id==153),'-','color',[0.4660 0.6740 0.1880]) %283
%S3A
geoplot(S3A_cor.Lat(S3A_cor.id==462),S3A_cor.Lon(S3A_cor.id==462),'-','color',[0 0.4470 0.7410]) %728
geoplot(S3A_cor.Lat(S3A_cor.id==290),S3A_cor.Lon(S3A_cor.id==290),'-','color',[0 0.4470 0.7410]) %283
%JA3
geoplot(JA3_cor.Lat(JA3_cor.id==1125),JA3_cor.Lon(JA3_cor.id==1125),'-r') %111
geoplot(JA3_cor.Lat(JA3_cor.id==1237),JA3_cor.Lon(JA3_cor.id==1237),'-r') %194

geoplot(lat_tg,lon_tg,'ok','MarkerSize',400)




ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');



%% Autoregression (AR)

ar=10; % 10, 15, 20, 25, 30
sys = arima(ar,0,0);
Md1 = estimate(sys,sla.sla_sa);
residual1 = infer(Md1,sla_tg);
prediction1 = sla_tg + residual1;

R = corrcoef(sla_tg,prediction1)
MAE=mae(sla_tg-prediction1)
RMSE=rms(sla_tg-prediction1)


subplot(2,1,1)
plot(date_tg,sla_tg,'k')
hold on
plot(date_tg,prediction1,'b')
ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);
set(gca,'xticklabel',{[]})

subplot(2,1,2)
plot(date_tg,sla_tg-prediction1,'b')
yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
ylabel('\DeltaSLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);


%% Moving Average
ma=25; %10 15 20 25
sys = arima(0,0,ma);
Md1 = estimate(sys,sla.sla_sa);
residual2 = infer(Md1,sla_tg);
prediction2 = sla_tg + residual2;

R = corrcoef(sla_tg,prediction2)
MAE=mae(sla_tg-prediction2)
RMSE=rms(sla_tg-prediction2)

subplot(2,1,1)
plot(date_tg,sla_tg,'k')
hold on
plot(date_tg,prediction2,'r')
ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);
set(gca,'xticklabel',{[]})

subplot(2,1,2)
plot(date_tg,sla_tg-prediction2,'r')
yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
ylabel('\DeltaSLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);



%% Autoregressive Moving Average (ARMA)
ar=25; %10, 15, 20, 25
ma=5; %2, 5

sys = arima(ar,0,ma);
Md1 = estimate(sys,sla.sla_sa);
residual3 = infer(Md1,sla_tg);
prediction3 = sla_tg + residual3;

R = corrcoef(sla_tg,prediction3)
MAE=mae(sla_tg-prediction3)
RMSE=rms(sla_tg-prediction3)

subplot(2,1,1)
plot(date_tg,sla_tg,'k')
hold on
plot(date_tg,prediction3,'color',[0.4660 0.6740 0.1880])
ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);
set(gca,'xticklabel',{[]})

subplot(2,1,2)
plot(date_tg,sla_tg-prediction3,'color',[0.4660 0.6740 0.1880])
yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
ylabel('\DeltaSLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);




%% Autoregressive Integrated Moving Average (ARIMA)
ar=10; %10, 15, 20
I=5; 
ma=2; %2, 5


sys = arima(ar,I,ma);
Md1 = estimate(sys,sla.sla_sa);
residual4 = infer(Md1,sla_tg);
prediction4 = sla_tg + residual4;

R = corrcoef(sla_tg,prediction4)
MAE=mae(sla_tg-prediction4)
RMSE=rms(sla_tg-prediction4)

subplot(2,1,1)
plot(date_tg,sla_tg,'k')
hold on
plot(date_tg,prediction4,'color',[0.9290 0.6940 0.1250])
ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);
set(gca,'xticklabel',{[]})

subplot(2,1,2)
plot(date_tg,sla_tg-prediction4,'color',[0.9290 0.6940 0.1250])
yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
ylabel('\DeltaSLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);




%% Seasonal Autoregressive Integrated Moving-Average (SARIMA)
sys = arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:2,'SARLags',[12,24,36,48],'Seasonality',24,'SMALags',24,'Distribution','Gaussian');
Md1 = estimate(sys,sla.sla_sa);
residual5 = infer(Md1,sla_tg);
prediction5 = sla_tg + residual5;


R = corrcoef(sla_tg,prediction5)
MAE=mae(sla_tg-prediction5)
RMSE=rms(sla_tg-prediction5)

subplot(2,1,1)
plot(date_tg,sla_tg,'k')
hold on
plot(date_tg,prediction5,'color',[0.9290 0.6940 0.1250])
ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);
set(gca,'xticklabel',{[]})

subplot(2,1,2)
plot(date_tg,sla_tg-prediction5,'color',[0.9290 0.6940 0.1250])
yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
ylabel('\DeltaSLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);



%% GARCH Model
GARCH_X1 = garch('Offset',0,'GARCHLags',1:14,'ARCHLags',50,'Distribution','Gaussian');
GARCH_X1 = estimate(GARCH_X1,sla.sla_sa,'Display','off');
residual6 = infer(Md1,sla_tg);
prediction6 = sla_tg + residual6;



R = corrcoef(sla_tg,prediction6)
MAE=mae(sla_tg-prediction6)
RMSE=rms(sla_tg-prediction6)

subplot(2,1,1)
plot(date_tg,sla_tg,'k')
hold on
plot(date_tg,prediction6,'color',[0.9290 0.6940 0.1250])
ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);
set(gca,'xticklabel',{[]})

subplot(2,1,2)
plot(date_tg,sla_tg-prediction6,'color',[0.9290 0.6940 0.1250])
yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
ylabel('\DeltaSLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);


%% Glostan, Jagannathan and Runkle GARCH Model
% this give the same result as GRACH model in resampling
GJR_X1 = gjr('Offset',0,'GARCHLags',1:3,'ARCHLags',1,'LeverageLags',1,'Distribution','Gaussian');
GJR_X1 = estimate(GJR_X1,sla.sla_sa,'Display','off');
residual7 = infer(Md1,sla_tg);
prediction7 = sla_tg + residual7;


R = corrcoef(sla_tg,prediction7)
MAE=mae(sla_tg-prediction7)
RMSE=rms(sla_tg-prediction7)

subplot(2,1,1)
plot(date_tg,sla_tg,'k')
hold on
plot(date_tg,prediction7,'color',[0.9290 0.6940 0.1250])
ylabel('SLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);
set(gca,'xticklabel',{[]})

subplot(2,1,2)
plot(date_tg,sla_tg-prediction7,'color',[0.9290 0.6940 0.1250])
yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
ylabel('\DeltaSLA [cm]','FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
set(gca,'fontname','Times New Roman','FontSize',18);




%% 
hold on
subplot (2,3,1)
histogram(sla_tg-prediction1,'Normalization','probability')
xlim([-15 15])
subplot (2,3,2)
histogram(sla_tg-prediction2,'Normalization','probability')
xlim([-15 15])
subplot (2,3,3)
histogram(sla_tg-prediction3,'Normalization','probability')
xlim([-15 15])
subplot (2,3,4)
histogram(sla_tg-prediction4,'Normalization','probability')
xlim([-15 15])
subplot (2,3,5)
histogram(sla_tg-prediction5,'Normalization','probability')
xlim([-15 15])
subplot (2,3,6)
histogram(sla_tg-prediction6,'Normalization','probability')
xlim([-15 15])

%% Long Short-Term Memory Networks (LSTM)

numTimeStepsTrain = floor(0.9*numel(sla.sla_sa));
dataTrain = sla.sla_sa(1:numTimeStepsTrain+1);
dataTest = sla.sla_sa(numTimeStepsTrain+1:end);

% Standardize Data

mu = mean(dataTrain);
sig = std(dataTrain);
dataTrainStandardized = (dataTrain - mu) / sig;

% Prepare Predictors and Responses
XTrain = dataTrainStandardized(1:end-1)';
YTrain = dataTrainStandardized(2:end)';

% Define LSTM Network Architecture
numFeatures = 1;
numResponses = 1;
numHiddenUnits = 200;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',250, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');
% Train LSTM Network
net = trainNetwork(XTrain,YTrain,layers,options);

% Forecast Future Time Steps
dataTestStandardized = (dataTest - mu) / sig;
XTest = dataTestStandardized(1:end-1);
net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,YTrain(end));

numTimeStepsTest = numel(XTest);
for i = 2:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
end

% Unstandardize the predictions using the parameters calculated earlier.
YPred = sig*YPred + mu;

YTest = dataTest(2:end);
rmse = sqrt(mean((YPred-YTest').^2))

figure
plot(dataTrain(1:end-1))
hold on
idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
plot(idx,[data(numTimeStepsTrain) YPred],'.-')
hold off
ylabel("SLA")
title("Forecast")
legend(["Observed" "Forecast"])

figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Forecast"])
ylabel("Cases")
title("Forecast")

subplot(2,1,2)
stem(YPred - YTest')
xlabel("Month")
ylabel("Error")
title("RMSE = " + rmse)



% Update Network State with Observed Values
net = resetState(net);
net = predictAndUpdateState(net,XTrain);

YPred = [];
XTest=XTest';
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
end

% Unstandardize the predictions using the parameters calculated earlier.
YPred = sig*YPred + mu;
rmse = sqrt(mean((YPred-YTest').^2))

figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Predicted"])
ylabel("Cases")
title("Forecast with Updates")

subplot(2,1,2)
stem(YPred - YTest')
xlabel("Month")
ylabel("Error")
title("RMSE = " + rmse)

%% Extreme Learning Machine (ELM)
 
 % II. Training set/test set generation
 
 % training set
P_train = NIR(sla.sla_sa(1:numTimeStepsTrain+1),:)';
T_train = octane(sla.sla_sa(1:numTimeStepsTrain+1),:)';
 
 % test set
P_test = NIR(sla.sla_sa(numTimeStepsTrain+1:end),:)';
T_test = octane(sla.sla_sa(numTimeStepsTrain+1:end),:)';
N = size(P_test,2);
 
% III. Data normalization
% 1. Training set
[Pn_train,inputps] = mapminmax(P_train);
Pn_test = mapminmax('apply',P_test,inputps);
% 2. Test set
[Tn_train,outputps] = mapminmax(T_train);
Tn_test = mapminmax('apply',T_test,outputps);
 
% IV. ELM creation/training
[IW,B,LW,TF,TYPE] = elmtrain(Pn_train,Tn_train,30,'sig',0);
 
% V. ELM simulation test
tn_sim = elmpredict(Pn_test,IW,B,LW,TF,TYPE);
% 1. Anti-normalization
T_sim = mapminmax('reverse',tn_sim,outputps);
 
% VI. Comparison of results
result = [T_test' T_sim'];
% 1. Mean square error
E = mse(T_sim - T_test);

% 2. Determination coefficient
N = length(T_test);
R2=(N*sum(T_sim.*T_test)-sum(T_sim)*sum(T_test))^2/((N*sum((T_sim).^2)-(sum(T_sim))^2)*(N*sum((T_test).^2)-(sum(T_test))^2)); 

% VII. Drawing
figure(1)
plot(1:N,T_test,'r-*',1:N,T_sim,'b:o')
grid on
Legend('true value', 'predicted value')
Xlabel('sample number')
Ylabel('octane number')
String = {'ELM';['(mse = ' num2str(E) ' R^2 = ' num2str(R2) ')']};
title(string)
