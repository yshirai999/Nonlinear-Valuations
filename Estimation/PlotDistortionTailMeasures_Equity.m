%% Digital Moments
clear
clc
close all

%% Load Data
DataPath = getPath('Data');

Y = load(fullfile(DataPath, 'Y'));
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

VarPath = getPath('VarArchive');

SPYMDGMM = load(fullfile(VarPath, strcat('SPYMD',num2str(SY),'GMM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); SPYMDGMM = SPYMDGMM.SPYMD;
SPYMDDM = load(fullfile(VarPath, strcat('SPYMD',num2str(SY),'DM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); SPYMDDM = SPYMDDM.SPYMD;
RCDM = load(fullfile(VarPath, strcat('RC',num2str(SY),'DM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); RCDM = RCDM.RC;
mu_mod_DM = load(fullfile(VarPath, strcat('mu_modDM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); mu_mod_DM = mu_mod_DM.mu_mod;
RCGMM = load(fullfile(VarPath, strcat('RC',num2str(SY),'GMM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); RCGMM = RCGMM.RC;
mu_mod_GMM = load(fullfile(VarPath, strcat('mu_modGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); mu_mod_GMM = mu_mod_GMM.mu_mod;
mu_mkt = load(fullfile(VarPath, strcat('mu_mkt',num2str(SY),num2str(Delta),num2str(N)))); mu_mkt = mu_mkt.mu_mkt;

DeltaVec1 = zeros(2,length(imin+5:5:imax));
DeltaVec2 = zeros(2,length(imin+5:5:imax));
mu_u = zeros(2,length(imin+5:5:imax));
mu = zeros(2,length(imin+5:5:imax));
mu_l = zeros(2,length(imin+5:5:imax));
std_u = zeros(2,length(imin+5:5:imax));
std_l = zeros(2,length(imin+5:5:imax));
skew = zeros(2,length(imin+5:5:imax));
kur_l = zeros(2,length(imin+5:5:imax));
kur_u = zeros(2,length(imin+5:5:imax));

theta = zeros(length(imin+5:5:imax),12);

Delta_y = 0.0001;
y = Delta_y:Delta_y:0.5;

for k=2:553%i = imin+5:5:imax
    i = imin+5*(k-1);
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
    
    theta(k,:) = [SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)];
    bp = theta(k,9); cp = theta(k,10); bn = theta(k,11); cn = theta(k,12);
    c = theta(k,1); gamma = theta(k,2); b = theta(k,3); a = theta(k,4);
    Gammam = exp(y).*(b/c).*(1-exp(-c*cp*expint(y/bp)))+exp(-y).*a.*(1-exp(-c*cn*expint(y/bn))).^(1/(1+gamma)); % Gamma_- Distorted Right Tail
    Gammap = exp(-y).*(b/c).*(1-exp(-c*cn*expint(y/bn)))+exp(y).*a.*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma));%+cp*expint(y/bp);%-a*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma)); % Gamma_- Distorted Left Tail
    %     mu_u(1,k) = sum(Gammap)*Delta_y;
    %     mu_u(1,k) = RCU([c,gamma,b,a,bp,cp,bn,cn]);
    mu_u(1,k) = mu_mod_GMM(k,1);
    mu(1,k) = cp*sum(exp(y).*expint(y/bp))*Delta_y-cn*sum(exp(-y).*expint(y/bn))*Delta_y;
    mu_l(1,k) = mu_mod_GMM(k,2);
    %     mu_l(1,k) = RCL([c,gamma,b,a,bp,cp,bn,cn]);
    %     mu_l(1,k) = sum(Gammam)*Delta_y;
    std_u(1,k) = sqrt(sum(y.*y.*Gammap)*Delta_y);
    std_l(1,k) = sqrt(sum(y.*y.*Gammam)*Delta_y);
    
    c = theta(k,5); gamma = theta(k,6); b = theta(k,7); a = theta(k,8);
    Gammam = exp(y).*(b/c).*(1-exp(-c*cp*expint(y/bp)))+exp(-y).*a.*(1-exp(-c*cn*expint(y/bn))).^(1/(1+gamma)); % Gamma_- Distorted Right Tail
    Gammap = exp(-y).*(b/c).*(1-exp(-c*cn*expint(y/bn)))+exp(y).*a.*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma));%+cp*expint(y/bp);%-a*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma)); % Gamma_- Distorted Left Tail
%     mu_u(2,k) = sum(Gammap)*Delta_y;
%     mu_u(2,k) = RCU([c,gamma,b,a,bp,cp,bn,cn]);
    mu_u(2,k) = mu_mod_DM(k,1);
    mu(2,k) = cp*sum(exp(y).*expint(y/bp))*Delta_y-cn*sum(exp(-y).*expint(y/bn))*Delta_y;
    mu_l(2,k) = mu_mod_DM(k,2);
%     mu_l(2,k) = RCL([c,gamma,b,a,bp,cp,bn,cn]);
%     mu_l(2,k) = sum(Gammam)*Delta_y;
    std_u(2,k) = sqrt(sum(y.*y.*Gammap)*Delta_y);
    std_l(2,k) = sqrt(sum(y.*y.*Gammam)*Delta_y);

%     theta = [SPYMDGMM(k,:),SPYMDDM(k,:),Y(i,5:8)];
%     bp = theta(9); cp = theta(10); bn = theta(11); cn = theta(12);
%     c = theta(1); gamma = theta(2); b = theta(3); a = theta(4);
%     Gammam = b*exp(-c*cp*expint(y/bp))...
%         -((a*c*exp(-c*cn*expint(y/bn)))/(1+gamma)).*(1-exp(-c*cn*expint(y/bn))).^(gamma/(1+gamma)); % Gamma_- Distorted Right Tail
%     Gammap = -((a*c*exp(-c*cp*expint(y/bp)))/(1+gamma)).*(1-exp(-c*cp*expint(y/bp))).^(gamma/(1+gamma))...
%         -b*(exp(y)-1).*exp(-c*cn*expint(y/bn)); % Gamma_- Distorted Left Tail
%     mu_u(1,k) = sum((exp(y)-1).*Gammap)*Delta_y;
%     mu(1,k) = ((1-bp).^(-cp))*((1+cp).^(-cn));
%     mu_l(1,k) = sum(exp(y).*Gammam)*Delta_y;
%     std_u(1,k) = sqrt(sum(y.*y.*Gammap)*Delta_y);
%     std_l(1,k) = sqrt(sum(y.*y.*Gammam)*Delta_y);
%     
%     c = theta(5); gamma = theta(6); b = theta(7); a = theta(8);
%     Gammam = b*exp(-c*cp*expint(y/bp))-((a*c*exp(-c*cn*expint(y/bn)))/(1+gamma)).*(1-exp(-c*cn*expint(y/bn))).^(1/(1+gamma)); % Gamma_- Distorted Right Tail
%     Gammap = ((a*c*exp(-c*cp*expint(y/bp)))/(1+gamma)).*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma))-b*exp(-c*cn*expint(y/bn)); % Gamma_- Distorted Left Tail
%     mu_u(2,k) = sum((exp(y)-1).*Gammap)*Delta_y;
%     mu(2,k) = -cn*sum((exp(y)-1).*expint(y/bn))*Delta_y + cp*sum((exp(y)-1).*expint(y/bp))*Delta_y;
%     mu_l(2,k) = sum((exp(y)-1).*Gammam)*Delta_y;
%     std_u(2,k) = sqrt(sum(y.*y.*Gammap)*Delta_y);
%     std_l(2,k) = sqrt(sum(y.*y.*Gammam)*Delta_y);
    
    mu_rel = mu(1,k)/mu(2,k)-1;
    fprintf('step = %d: relative error = %d, Upper = %d, mid = %d, Lower = %d\n',k, mu_rel, mu_u(2,k),mu(2,k),mu_l(2,k))
end  

%% Visualize Mean
q = 10;
muGMM_u_out = rmoutliers(mu_u(1,2:end),'percentiles',[q,100-q]);
muGMM_out = rmoutliers(mu(1,2:end),'percentiles',[q,100-q]);
muGMM_l_out = rmoutliers(mu_l(1,2:end),'percentiles',[q,100-q]);

muDM_u_out = rmoutliers(mu_u(2,2:end),'percentiles',[q,100-q]);
muDM_out = rmoutliers(mu(2,2:end),'percentiles',[q,100-q]);
muDM_l_out = rmoutliers(mu_l(2,2:end),'percentiles',[q,100-q]);

quant = [0.25;0.5;0.75];

QmuDMu = quantile(muDM_u_out,quant);
QmuDM = quantile(muDM_out,quant);
QmuDMl = quantile(muDM_l_out,quant);
QmuGMMu = quantile(muGMM_u_out,quant);
QmuGMM = quantile(muGMM_out,quant);
QmuGMMl = quantile(muGMM_l_out,quant);
fprintf('mu DM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f & %3.4f & %3.4f \\\\ \n', quant(i), QmuDMu(i,:), QmuDM(i,:), QmuDMl(i,:))
end
fprintf('mu GMM quantiles\n')
for i=1:length(quant)
    fprintf('%3.2f & %3.4f & %3.4f & %3.4f \\\\ \n', quant(i), QmuGMMu(i,:), QmuGMM(i,:), QmuGMMl(i,:))
end

[MDQDM, ~, ~] = NonlinearPricing.Functions.vqsplit([muDM_u_out;muDM_out;muDM_l_out],2^4);
MDQDM = 100*MDQDM'*252/Delta;

[MDQGMM, ~, ~] = NonlinearPricing.Functions.vqsplit([muGMM_u_out;muGMM_out;muGMM_l_out],2^4);
MDQGMM = 100*MDQGMM'*252/Delta;
fprintf('mu quantized \n')

fprintf('mu quantized \n')

for i=1:2^4
    fprintf('%3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f \\\\ \n', MDQDM(i,:), MDQGMM(i,:))
end

[MDQ, rep, ~] = vqsplit([muDM_u_out;muDM_out;muDM_l_out;muGMM_u_out;muGMM_out;muGMM_l_out],2^4);
MDQ = 100*MDQ'*252/Delta;

fprintf('\nmu quantized \n')

for i=1:2^4
    fprintf('%3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f \\\\ \n', MDQ(i,:))
end

%% Quantize Parameters
ParamQ = rmoutliers(theta(2:end,:),'percentiles',[q,100-q]);
size(ParamQ)
[ParamQ, ~, ~] = vqsplit(ParamQ',2^4);
ParamQ = ParamQ';
fprintf('mu quantized \n')
for i=1:2^4
    fprintf('%3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f \\\\ \n', ParamQ(i,:))
end

%% Visualization

VizPath = getPath('Visualization');

dates = datenum(num2str(Y(imin+5:5:imax,1)),'yyyymmdd');
dates = datetime(dates,'ConvertFrom','datenum');
dates = dates(1:end-1);

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
saveas(gcf, fullfile(VizPath, fname), 'epsc');

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
saveas(gcf, fullfile(VizPath, fname), 'epsc');

% %% Visualize Dispersion
% q = 0;
% stdGMM_u_out = rmoutliers(std_u(1,2:end),'percentiles',[q,100-q]);
% stdGMM_l_out = rmoutliers(std_l(1,2:end),'percentiles',[q,100-q]);
% 
% stdDM_u_out = rmoutliers(std_u(2,2:end),'percentiles',[q,100-q]);
% stdDM_l_out = rmoutliers(std_l(2,2:end),'percentiles',[q,100-q]);
% 
% quant = [0.25;0.5;0.75];
% 
% QstdDMu = quantile([stdDM_u_out],quant);
% QstdDMl = quantile([stdDM_l_out],quant);
% QstdGMMu = quantile([stdGMM_u_out],quant);
% QstdGMMl = quantile([stdGMM_l_out],quant);
% fprintf('std DM quantiles\n')
% for i=1:length(quant)
%     fprintf('%3.2f & %3.4f & %3.4f \\\\ \n', quant(i), QstdDMu(i,:), QstdDMl(i,:))
% end
% fprintf('std GMM quantiles\n')
% for i=1:length(quant)
%     fprintf('%3.2f & %3.4f & %3.4f \\\\ \n', quant(i), QstdGMMu(i,:), QstdGMMl(i,:))
% end
% [MDQ, ~, ~] = NonlinearPricing.Functions.vqsplit([stdDM_u_out;stdDM_l_out;stdGMM_u_out;stdGMM_l_out],2^4);
% MDQ = MDQ';
% fprintf('mu quantized \n')
% for i=1:10
%     fprintf('%3.4f & %3.4f & %3.4f & %3.4f \\\\ \n', MDQ(i,:))
% end

%% Visualization

fprintf('\n\n\n')
for k = 400:553
    Delta = 0.0001;
    y = Delta:Delta:0.05;
    c = theta(k,5); gamma = theta(k,6); b = theta(k,7); a = theta(k,8);
    bp = theta(k,9); cp = theta(k,10); bn = theta(k,11); cn = theta(k,12);
    Gammam = (b/c)*exp(y).*(1-exp(-c*cp*expint(y/bp)))+a*exp(-y).*(1-exp(-c*cn*expint(y/bn))).^(1/(1+gamma)); % Gamma_- Distorted Right Tail
    nup =  cp*exp(y).*expint(y/bp);
    
    nun =  exp(-y).*cn.*expint(y/bn);
    Gammap = (b/c)*exp(-y).*(1-exp(-c*cn*expint(y/bn)))+a*exp(y).*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma));%+cp*expint(y/bp);%-a*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma)); % Gamma_- Distorted Left Tail

%     figure
%     hold on
%     box on
%     plot(y,nup)
%     plot(y,Gammam)
%     legend('nu_p','Gammam')
%     
%     figure
%     hold on
%     box on
%     plot(y,nun)
%     plot(y,Gammap)
%     legend('nu_n','Gammap')
    
    %fprintf('Integrals = %d, %d, %d\n', -sum(nun)*Delta+sum(nup)*Delta, sum(Gammam)*Delta, sum(Gammap)*Delta)
end

%% Routines

function I = RCU(theta)
    c = abs(theta(1)); gamma = abs(theta(2));
    b = abs(theta(3));
    a = abs(theta(4));
    bp = abs(theta(5)); cp = abs(theta(6)); bn = abs(theta(7)); cn = abs(theta(8)); 
    
    L = @(w) cp*expint(log(1+w)/bp);
    K = @(w) cn*expint(-log(1-w)/bn);
    Gp = @(y) a*(1-exp(-c*y)).^(1/(1+gamma));
    Gm = @(y) (b/c)*(1-exp(-c*y));
    funp = @(w) Gp(L(w));
    funm = @(w) Gm(K(w));

    Ip = integral(funm,0,1);
    Im = integral(funp,0,Inf);
    I = Ip+Im;
end

function I = RCL(theta)
    c = abs(theta(1)); gamma =  abs(theta(2));
    b = abs(theta(3)); 
    a = abs(theta(4));
    bp = abs(theta(5)); cp = abs(theta(6)); bn = abs(theta(7)); cn = abs(theta(8)); 
    
    L = @(w) cp*expint(log(1+w)/bp);
    K = @(w) cn*expint(-log(1-w)/bn);
    Gp = @(y) a*(1-exp(-c*y)).^(1/(1+gamma));
    Gm = @(y) (b/c)*(1-exp(-c*y));
    funp = @(w) Gm(L(w));
    funm = @(w) Gp(K(w));

    Ip = integral(funm,0,1);
    Im = integral(funp,0,Inf);
    I = Ip+Im;
end