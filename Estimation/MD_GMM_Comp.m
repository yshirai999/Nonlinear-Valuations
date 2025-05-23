%% Digital Moments
clear
clc
close all

%% Load Data
dataPath = getPath('Data');
Y = load(fullfile(dataPath, 'Y'));
Y = Y.Y;

SY = 10; startdate = strcat('20',num2str(SY,'%02.f'),'0104');
SD = str2double(startdate);

N = 250;

enforcegamma = 1;
enforceb = 1;
%lb = 1e-4;
lb = 0;

%options = optimset('MaxFunEval',1000,'MaxIter',100,'Display','iter','TolFun',1e-10,'TolX',1e-10,'OutputFcn', @outfun);
options = optimset('MaxFunEval',5000,'MaxIter',1000,'Display','off','TolFun',1e-10,'TolX',1e-10);
optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',1000,'Display','iter','FunctionTolerance',1e-10);

Nsim = 1000;
q = qrandstream('halton',1,'Skip',1e3,'Leap',1e2);
U = qrand(q,Nsim);

%% Load Data
varPath = getPath('VarArchive');
imin = find(Y(:,1)==SD);
imax = length(Y);
Delta = 5;
n = length(imin:Delta:imax);
theta = zeros(n,8);
theta(:,5:8) = Y(imin:Delta:imax,5:8);
d = Y(imin:Delta:imax,1);

ci = 1;
gammai = 0.2;
if enforceb == 1 
    bi = -log(1);
else
    bi = 1;
end
ai = 0.001;

%(c,gamma,b,a) = (3.642385e+02,9.094947e-13,1.191809e+02,1.445863e-03)

SPYMDGMM = load(fullfile(varPath, strcat('SPYMD',num2str(SY),'GMM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); SPYMDGMM = SPYMDGMM.SPYMD;
SPYMDDM = load(fullfile(varPath, strcat('SPYMD',num2str(SY),'DM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); SPYMDDM = SPYMDDM.SPYMD;
DeltaVec1 = zeros(2,length(imin+5:5:imax));
DeltaVec2 = zeros(2,length(imin+5:5:imax));
mu_u = zeros(2,length(imin+5:5:imax));
mu_l = zeros(2,length(imin+5:5:imax));
std_u = zeros(2,length(imin+5:5:imax));
std_l = zeros(2,length(imin+5:5:imax));
skew = zeros(2,length(imin+5:5:imax));
kur_l = zeros(2,length(imin+5:5:imax));
kur_u = zeros(2,length(imin+5:5:imax));


k = 1;
for i = imin+5:5:imax
    k = k+1;
    ind = Delta:Delta:N;
    NN = length(ind);
    YY = zeros(NN+1,3);
    YY(1,1) = Y(i-N,1);
    YY(1,2) = min(Y(i-N-Delta+1:i-N,2));
    YY(1,3) = max(Y(i-N-Delta+1:i-N,3));
    YY(1,4) = mean([YY(1,2),YY(1,3)]);
    for j=2:NN+1
        YY(j,1) = Y(i-N+ind(j-1),1);
        YY(j,2) = min(Y(i-N+ind(j-1)-Delta+1:i-N+ind(j-1),2));
        YY(j,3) = max(Y(i-N+ind(j-1)-Delta+1:i-N+ind(j-1),3));
        YY(j,4) = mean([YY(j,2),YY(j,3)]);
    end
    u = ( log(YY(2:NN+1,3)') - log(YY(1:NN,3)') )/Delta;
    l = ( log(YY(2:NN+1,2)') - log(YY(1:NN,2)') )/Delta;
    yp = -log(U)*Y(end,5);
    yn = -log(U)*Y(end,7);

    if i+5 > imax
        [DeltaVec2(2,k),DeltaVec1(2,k),mu_l(1,k),std_l(1,k),mu_l(2,k),std_l(2,k),kur_l(1,k),kur_l(2,k)] = ...
            plotL(l,[SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)],NN,yp,yn,1);
        [DeltaVec2(1,k),DeltaVec1(1,k),mu_u(1,k),std_u(1,k),mu_u(2,k),std_u(2,k),kur_u(1,k),kur_u(2,k)] = ...
            plotU(u,[SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)],NN,yp,yn,1);
    else
        [DeltaVec2(2,k),DeltaVec1(2,k),mu_l(1,k),std_l(1,k),mu_l(2,k),std_l(2,k),kur_l(1,k),kur_l(2,k)] = ...
            plotL(l,[SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)],NN,yp,yn,0);
        [DeltaVec2(1,k),DeltaVec1(1,k),mu_u(1,k),std_u(1,k),mu_u(2,k),std_u(2,k),kur_u(1,k),kur_u(2,k)] = ...
            plotU(u,[SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)],NN,yp,yn,0);
    end

    fprintf('i = %d, k = %d: dL = (%d,%d), dU = (%d,%d) \n', i, k, DeltaVec1(2,k),DeltaVec2(2,k), DeltaVec1(1,k), DeltaVec2(1,k))

end  

%% Visualization
close all

vizPath = getPath('Visualization');
DeltaVec1_out = rmoutliers(DeltaVec1,'percentiles',[5,95]);
DeltaVec2_out = rmoutliers(DeltaVec2,'percentiles',[5,95]);

figure
hold on
box on
grid on
plot([1:1:length(DeltaVec1_out(1,:))],DeltaVec1_out(1,:))
title('L1 Distance GMM-fDM Lower densities')
set(gca,'TickLabelInterpreter','latex')
str=strcat('L1_Norm_DeltaLowerDensities');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');

figure
hold on
box on
grid on
plot([1:1:length(DeltaVec1_out(2,:))],DeltaVec1_out(2,:))
title('L1 Distance GMM-fDM Upper densities')
set(gca,'TickLabelInterpreter','latex')
str=strcat('L1_Norm_DeltaUpperDensities');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');


figure
hold on
box on
grid on
scatter(DeltaVec1_out(1,:),DeltaVec1_out(2,:))
title('L1 Distance GMM-fDM Upper densities')
set(gca,'TickLabelInterpreter','latex')
str=strcat('L1_Norm_DeltaDensitiesScatter');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');


figure
hold on
box on
grid on
plot([1:1:length(DeltaVec2_out(1,:))],DeltaVec2_out(1,:))
title('L2 Distance GMM-fDM Lower densities')
set(gca,'TickLabelInterpreter','latex')
str=strcat('L2_Norm_DeltaLowerDensities');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');

figure
hold on
box on
grid on
plot([1:1:length(DeltaVec2_out(2,:))],DeltaVec2_out(2,:))
title('L2 Distance GMM-fDM Upper densities')
set(gca,'TickLabelInterpreter','latex')
str=strcat('L2_Norm_DeltaUpperDensities');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');


figure
hold on
box on
grid on
scatter(DeltaVec2_out(1,:),DeltaVec2_out(2,:))
title('L2 Distance GMM-fDM Upper densities')
set(gca,'TickLabelInterpreter','latex')
str=strcat('L2_Norm_DeltaDensitiesScatter');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');

%% Table quantiles DeltaVec

quant = [0.25;0.5;0.75];

QDVec = quantile([DeltaVec1(1,:)],quant);
fprintf('Upper Delta GMM-DM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f \\\\ \n', quant(i), QDVec(i,:))
end

QDVec = quantile([DeltaVec1(2,:)],quant);
fprintf('Lower Delta GMM-DM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f \\\\ \n', quant(i), QDVec(i,:))
end

% fprintf('\n')
% QDVec = quantile([DeltaVec2(1,:)],quant);
% for i=1:length(quant)
%     fprintf('%3.2f & %3.4f \\\\ \n', quant(i), QDVec(i,:))
% end
% 
% QDVec = quantile([DeltaVec2(2,:)],quant);
% fprintf('\n')
% for i=1:length(quant)
%     fprintf('%3.2f & %3.4f \\\\ \n', quant(i), QDVec(i,:))
% end

%% Visualize Mean
q = 0;
muGMM_u_out = rmoutliers(mu_u(1,2:end),'percentiles',[q,100-q]);
muGMM_l_out = rmoutliers(mu_l(1,2:end),'percentiles',[q,100-q]);

muDM_u_out = rmoutliers(mu_u(2,2:end),'percentiles',[q,100-q]);
muDM_l_out = rmoutliers(mu_l(2,2:end),'percentiles',[q,100-q]);

quant = [0.25;0.5;0.75];

QmuDMu = quantile([muDM_u_out],quant);
QmuDMl = quantile([muDM_l_out],quant);
QmuGMMu = quantile([muGMM_u_out],quant);
QmuGMMl = quantile([muGMM_l_out],quant);
fprintf('mu DM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f & %3.4f \\\\ \n', quant(i), QmuDMu(i,:), QmuDMl(i,:))
end
fprintf('mu GMM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f & %3.4f \\\\ \n', quant(i), QmuGMMu(i,:), QmuGMMl(i,:))
end
%%
[MDQ, ~, ~] = vqsplit([muDM_u_out;muDM_l_out;muGMM_u_out;muGMM_l_out],2^4);
MDQ = MDQ';
MDQ(1:5,:)
%%
dates = datenum(num2str(Y([imin+5:5:imax],1)),'yyyymmdd');
dates = datetime(dates,'ConvertFrom','datenum');

figure
hold on
box on
grid on
plot(dates,muGMM_u_out)
plot(dates,muGMM_l_out)
title('Mean GMM densities')
legend('Upper','Lower')
set(gca,'TickLabelInterpreter','latex')
str=strcat('muGMM_Upper');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');

figure
hold on
box on
grid on
plot(dates,muDM_u_out)
plot(dates,muDM_l_out)
title('Mean DM densities')
legend('Upper','Lower')
set(gca,'TickLabelInterpreter','latex')
str=strcat('muDM_Lower');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');


%% Visualize Dispersion
q = 0;
stdGMM_u_out = rmoutliers(std_u(1,2:end),'percentiles',[q,100-q]);
stdGMM_l_out = rmoutliers(std_l(1,2:end),'percentiles',[q,100-q]);

stdDM_u_out = rmoutliers(std_u(2,2:end),'percentiles',[q,100-q]);
stdDM_l_out = rmoutliers(std_l(2,2:end),'percentiles',[q,100-q]);

quant = [0.25;0.5;0.75];

QstdDMu = quantile([stdDM_u_out],quant);
QstdDMl = quantile([stdDM_l_out],quant);
QstdGMMu = quantile([stdGMM_u_out],quant);
QstdGMMl = quantile([stdGMM_l_out],quant);
fprintf('std DM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f & %3.4f \\\\ \n', quant(i), QstdDMu(i,:), QstdDMl(i,:))
end
fprintf('std GMM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f & %3.4f \\\\ \n', quant(i), QstdGMMu(i,:), QstdGMMl(i,:))
end

dates = datenum(num2str(Y([imin+5:5:imax],1)),'yyyymmdd');
dates = datetime(dates,'ConvertFrom','datenum');

figure
hold on
box on
grid on
plot(dates,stdGMM_u_out)
plot(dates,stdGMM_l_out)
title('Standard Deviation GMM densities')
legend('Upper','Lower')
set(gca,'TickLabelInterpreter','latex')
str=strcat('stdGMM_Upper');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');

figure
hold on
box on
grid on
plot(dates,stdDM_u_out)
plot(dates,stdDM_l_out)
title('Standard Deviation DM densities')
legend('Upper','Lower')
set(gca,'TickLabelInterpreter','latex')
str=strcat('stdDM_Lower');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');

%% Visualize Kurtosis
kurGMM_u_out = rmoutliers(kur_u(1,:),'percentiles',[10,90]);
kurGMM_l_out = rmoutliers(kur_l(1,:),'percentiles',[10,90]);

kurDM_u_out = rmoutliers(kur_u(2,:),'percentiles',[10,90]);
kurDM_l_out = rmoutliers(kur_l(2,:),'percentiles',[10,90]);

figure
hold on
box on
grid on
plot([1:1:length(kurGMM_u_out)],kurGMM_u_out)
plot([1:1:length(kurGMM_l_out)],kurGMM_l_out)
title('Kurtosis GMM densities')
legend('Upper','Lower')
set(gca,'TickLabelInterpreter','latex')
str=strcat('KurtosisGMM_Upper');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');

figure
hold on
box on
grid on
plot([1:1:length(kurDM_u_out)],kurDM_u_out)
plot([1:1:length(kurDM_l_out)],kurDM_l_out)
title('Kurtosis DM densities')
legend('Upper','Lower')
set(gca,'TickLabelInterpreter','latex')
str=strcat('KurtosisDM_Lower');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');

%% Visualize Specific Day

k = 510;
i = imin + 5*(k-1);

[DeltaVec2(2,k),DeltaVec1(2,k)] = plotL(l,[SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)],NN,yp,yn,1);
[DeltaVec2(1,k),DeltaVec1(1,k)] = plotU(u,[SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)],NN,yp,yn,1);

fprintf('i = %d, k = %d: dL = (%d,%d), dU = (%d,%d) \n', i, k, DeltaVec1(2,k),DeltaVec2(2,k), DeltaVec1(1,k), DeltaVec2(1,k))

%% Routines

    function [d1,d2,muGMM,sigmaGMM,muDM,sigmaDM,kGMM,kDM] = plotU(u,theta,N,yp,yn,pl)
    s = u;% + RCU(theta);
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,999,999)/1000;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(9); cp = theta(10); bn = theta(11); cn = theta(12);
    c = theta(1); gamma = theta(2); b = -theta(3); a = theta(4);
    Nsim = length(yp);
    B = 10000;
    NF=2^10;
    eta=B/NF;
    lambda = 2*pi/B;
    bb = lambda*NF/2;
    u=[0:NF-1]*eta;
    w = ones(1,NF); w(1)=1/2;
    x = -bb+lambda*[0:NF-1];
    phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
    
    phi = phi.* exp( sum(...
            (exp(1i*yp*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma)) )...
                .*  bp*cp .* exp(-c*cp*expint(yp/bp)) ./ yp...
            -(exp(-1i*yn*u)-1) .*  b*bn*cn .* exp(-c*cn*expint(yn/bn)) ./ yn) /Nsim);
    fGMM = max(real(fft((1/pi)*exp(1i*u*bb).*phi*eta)),0);
    
    muGMM = sum(x.*fGMM)*lambda;
    sigmaGMM = sqrt(sum( ((x-muGMM).^2) .* fGMM)*lambda);
    kGMM = sum( ( (x-muGMM)./sigmaGMM ).^4 )*lambda;
    
    Delta = 0.01;
    y = Delta:Delta:5;
    %Gammap = (b/c)*exp(-c*cp*expint(y/bp));
    Gammam = (b/c)*exp(-c*cn*expint(y/bn));
    %muGMMp = sum(Gammap)*Delta;
    muGMMm = sum(Gammam)*Delta;
    muGMM = muGMMm; % left tail measure distorted by Gamma_-

    c = theta(5); gamma = theta(6); b = -theta(7); a = theta(8);
    Nsim = length(yp);
    B = 10000;
    NF=2^10;
    eta=B/NF;
    lambda = 2*pi/B;
    bb = lambda*NF/2;
    u=[0:NF-1]*eta;
    w = ones(1,NF); w(1)=1/2;
    x = -bb+lambda*[0:NF-1];
    phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
    
    phi = phi.* exp( sum(...
            (exp(1i*yp*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma)) )...
                .*  bp*cp .* exp(-c*cp*expint(yp/bp)) ./ yp...
            -(exp(-1i*yn*u)-1) .*  b*bn*cn .* exp(-c*cn*expint(yn/bn)) ./ yn) /Nsim);
    fDM = max(real(fft((1/pi)*exp(1i*u*bb).*phi*eta)),0);
    
    muDM = sum(x.*fDM)*lambda;
    sigmaDM = sqrt(sum( ((x-muDM).^2) .* fGMM)*lambda);
    kDM = sum( ( (x-muDM)./sigmaDM ).^4 )*lambda;
    
    Delta = 0.01;
    y = Delta:Delta:5;
    %Gammap = (b/c)*exp(-c*cp*expint(y/bp));
    Gammam = (b/c)*exp(-c*cn*expint(y/bn));
    %muDMp = sum(Gammap)*Delta;
    muDMm = sum(Gammam)*Delta;
    muDM = muDMm; % left tail measure distorted by Gamma_-

    f = [0,max(interp1(x,fGMM,s),0)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)];
    PihatGMM = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
    PihatGMM(s>0) = 1-PihatGMM(s>0);

    f = [0,max(interp1(x,fDM,s),0)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)];
    PihatDM = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
    PihatDM(s>0) = 1-PihatDM(s>0);

    if pl == 1
        figure
        hold on
        grid on
        box on
        plot(x,fGMM,':','LineWidth',1);
        plot(x,fDM,'--','LineWidth',1)
        title('Upper Density')
        legend('GMM', 'DM')
        set(gca,'TickLabelInterpreter','latex')
        str=strcat('SPY_UpperDensity');
        fname=str;
        saveas(gcf, fullfile(vizPath, fname), 'epsc');

        figure
        hold on
        grid on
        box on
        plot(x,(fGMM-fDM),'LineWidth',1);
        %plot(x,fDM,':','LineWidth',1)
        title('Difference between GMM and DM Upper Density')
        %legend('GMM', 'DM')
        set(gca,'TickLabelInterpreter','latex')
        str=strcat('SPY_UpperDensityDelta');
        fname=str;
        saveas(gcf, fullfile(vizPath, fname), 'epsc');

        figure
        hold on
        grid on
        box on
        plot(s,max(PihatGMM,0),'--')
        plot(s,max(PihatDM,0),'--')
        plot(s,Pi,'-')
        legend('Estimated GMM', 'Estimated DM', 'Empirical','interpreter','latex')
        title('Upper Digital Moments Fitting')
        set(gca,'TickLabelInterpreter','latex')
        str=strcat('SPY_UpperGMMFitting');
        fname=str;
        saveas(gcf, fullfile(vizPath, fname), 'epsc');
    end
    
    d1 = norm(fGMM-fDM,1);
    d2 = norm(fGMM-fDM,2);
end

function [d1,d2,muGMM,sigmaGMM,muDM,sigmaDM,kGMM,kDM] = plotL(l,theta,N,yp,yn,pl)
    s = l;% - RCL(theta);
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,999,999)/1000;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(9); cp = theta(10); bn = theta(11); cn = theta(12);
    c = theta(1); gamma = theta(2); b = theta(3); a = theta(4);
    Nsim = length(yp);
    B = 10000;
    NF=2^10;
    eta=B/NF;
    lambda = 2*pi/B;
    bb = lambda*NF/2;
    u=[0:NF-1]*eta;
    w = ones(1,NF); w(1)=1/2;
    x = -bb+lambda*[0:NF-1];
    phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
    
    phi = phi .* exp( sum(...
            -(exp(1i*yp*u)-1) .*  b*bp*cp .* exp(-c*cp*expint(yp/bp)) ./ yp...
            +(exp(-1i*yn*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cn*expint(yn/bn))).^(-gamma/(1+gamma)) )...
                .*  bn*cn .* exp(-c*cn*expint(yn/bn)) ./ yn) /Nsim);
    fGMM = max(real(fft((1/pi)*exp(1i*u*bb).*phi*eta)),0);
    
    muGMM = sum(x.*fGMM)*lambda;
    sigmaGMM = sqrt(sum( ((x-muGMM).^2) .* fGMM)*lambda);
    kGMM = sum( ( (x-muGMM)./sigmaGMM ).^4 )*lambda;
    
    Delta = 0.01;
    y = Delta:Delta:5;
    Gammap = (b/c)*exp(-c*cp*expint(y/bp));
    %Gammam = (b/c)*exp(-c*cn*expint(y/bn));
    muGMMp = sum(Gammap)*Delta;
    %muGMMm = sum(Gammam)*Delta;
    muGMM = muGMMp; % right tail measure distorted by Gamma_-

    c = theta(5); gamma = theta(6); b = theta(7); a = theta(8);
    Nsim = length(yp);
    B = 10000;
    NF=2^10;
    eta=B/NF;
    lambda = 2*pi/B;
    bb = lambda*NF/2;
    u=[0:NF-1]*eta;
    w = ones(1,NF); w(1)=1/2;
    x = -bb+lambda*[0:NF-1];
    phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
    
    phi = phi .* exp( sum(...
            -(exp(1i*yp*u)-1) .*  b*bp*cp .* exp(-c*cp*expint(yp/bp)) ./ yp...
            +(exp(-1i*yn*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cn*expint(yn/bn))).^(-gamma/(1+gamma)) )...
                .*  bn*cn .* exp(-c*cn*expint(yn/bn)) ./ yn) /Nsim);
    fDM = max(real(fft((1/pi)*exp(1i*u*bb).*phi*eta)),0);
    
    muDM = sum(x.*fDM)*lambda;
    sigmaDM = sqrt(sum( ((x-muDM).^2) .* fDM)*lambda);
    kDM = sum( ( (x-muDM)./sigmaDM ).^4 )*lambda;
    
    Delta = 0.01;
    y = Delta:Delta:5;
    Gammap = (b/c)*exp(-c*cp*expint(y/bp));
    %Gammam = (b/c)*exp(-c*cn*expint(y/bn));
    muDMp = sum(Gammap)*Delta;
    %muDMm = sum(Gammam)*Delta;
    muDM = muDMp; % right tail measure distorted by Gamma_-


    f = [0,max(interp1(x,fGMM,s),0)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)];     
    PihatGMM = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
%     figure
%     plot(s,Pihat)
    PihatGMM(s>0) = 1-PihatGMM(s>0);

    f = [0,max(interp1(x,fDM,s),0)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)];     
    PihatDM = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
%     figure
%     plot(s,Pihat)
    PihatDM(s>0) = 1-PihatDM(s>0);

    if pl == 1
        figure
        hold on
        grid on
        box on
        plot(x,fGMM,':','LineWidth',1);
        plot(x,fDM,'--','LineWidth',1)
        title('Lower Density')
        legend('GMM', 'DM')
        set(gca,'TickLabelInterpreter','latex')
        str=strcat('SPY_LowerDensity');
        fname=str;
        saveas(gcf, fullfile(vizPath, fname), 'epsc');

        figure
        hold on
        grid on
        box on
        plot(x,(fGMM-fDM),'LineWidth',1);
        %plot(x,fDM,':','LineWidth',1)
        title('Difference between GMM and DM Lower Density')
        %legend('GMM', 'DM')
        set(gca,'TickLabelInterpreter','latex')
        str=strcat('SPY_LowerDensityDelta');
        fname=str;
        saveas(gcf, fullfile(vizPath, fname), 'epsc');

        figure
        hold on
        grid on
        box on
        plot(s,max(PihatGMM,0),'--')
        plot(s,max(PihatDM,0),'--')
        plot(s,Pi,'-')
        legend('Estimated GMM', 'Estimated DM', 'Empirical','interpreter','latex')
        title('Lower Digital Moments Fitting')
        set(gca,'TickLabelInterpreter','latex')
        str=strcat('SPY_LowerGMMFitting');
        fname=str;
        saveas(gcf, fullfile(vizPath, fname), 'epsc');    
    end

    d1 = norm(fGMM-fDM,1);
    d2 = norm(fGMM-fDM,2);
end