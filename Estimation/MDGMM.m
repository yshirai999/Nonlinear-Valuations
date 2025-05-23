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

%% Generalized Method of Moments
imin = find(Y(:,1)==SD);
imax = length(Y);
Delta = 5;
n = length(imin:Delta:imax);
theta = zeros(n,8);
mu_mod = zeros(n,3);
mu_mkt = zeros(n,3);
err = zeros(n,1);
RC = zeros(n,2);
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

varPath = getPath('VarArchive');
try
    SPYMD = load(fullfile(varPath,strcat('SPYMD',num2str(SY),'GMM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); SPYMD = SPYMD.SPYMD;
    RC = load(fullfile(varPath,strcat('RC',num2str(SY),'GMM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); RC = RC.RC;
    mu_mod = load(fullfile(varPath,strcat('mu_modGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); mu_mod = mu_mod.mu_mod;
    mu_mkt = load(fullfile(varPath,strcat('mu_mkt',num2str(SY),num2str(Delta),num2str(N)))); mu_mkt = mu_mkt.mu_mkt;
    err = load(fullfile(varPath,strcat('errGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N)))); err = err.err;

    for i = imax
        ind = [Delta:Delta:N];
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
        plotL(l,[SPYMD(end,:),Y(end,5:8)],NN,yp,yn)
        plotU(u,[SPYMD(end,:),Y(end,5:8)],NN,yp,yn)
    end    
catch
    for i = imin+Delta:Delta:imax
        tic
        theta0 = [ci,gammai,bi,ai];
        ind = [Delta:Delta:N];
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
        m = ( log(Y(i-N+1:i,4)') - log(Y(i-N:i-1,4)') );
    
        bp = Y(i,5);
        cp = Y(i,6);
        bn = Y(i,7);
        cn = Y(i,8);
        
        
        
        
        obj = @(theta) objU(u,abs(theta),bp,cp,bn,cn,enforcegamma,enforceb,lb)+objL(l,abs(theta),bp,cp,bn,cn,enforcegamma,enforceb,lb);
        
        [theta((i-imin)/Delta+1,1:4),err((i-imin)/Delta+1)] = fminsearch(obj,theta0(1:4),options);        
        gamma = theta((i-imin)/Delta+1,2);
        if enforcegamma == 1
            range = 1-lb;
            if gamma>1||gamma<lb
                n = floor((gamma-lb)/range);
                if mod(n,2)==1
                    gamma = gamma-n*range;
                else
                    gamma = 1+n*range-(gamma-lb);
                end
            end
            theta((i-imin)/Delta+1,2) = gamma;
        end
        
        if enforceb == 1
            theta((i-imin)/Delta+1,3) = exp(-abs(theta((i-imin)/Delta+1,3)));
        end
        
        theta = abs(theta);
    
        mu = ((1-bp)^(-cp))*(1+bn)^(-cn); 
        RC((i-imin)/Delta+1,:) = [RCU(theta((i-imin)/Delta+1,:)), RCL(theta((i-imin)/Delta+1,:))];
        mu_mod((i-imin)/Delta+1,:) = [log(mu)-RC((i-imin)/Delta+1,1),log(mu),log(mu)+RC((i-imin)/Delta+1,2)];
        mu_mkt((i-imin)/Delta+1,:) = [mean(u),mean(m),mean(l)];
        cputime = toc;
        fprintf('Date = %d: (c,gamma,b,a) = (%d,%d,%d,%d)\n', Y(i,1), theta((i-imin)/Delta+1,1:4))
        fprintf('Date = %d: (bp,cp,bn,cn) = (%d,%d,%d,%d)\n', Y(i,1), theta((i-imin)/Delta+1,5:8))
        fprintf('Date = %d: mu_u = %d, mu = %d, mu_l = %d\n', Y(i,1), mu_mod((i-imin)/Delta+1,1), mu_mod((i-imin)/Delta+1,2), mu_mod((i-imin)/Delta+1,3))
        fprintf('Date = %d: mu_u = %d, mu = %d, mu_l = %d\n', Y(i,1), mu_mkt((i-imin)/Delta+1,1), mu_mkt((i-imin)/Delta+1,2), mu_mkt((i-imin)/Delta+1,3))
        fprintf('Error = %d\n', err((i-imin)/Delta+1))
        fprintf('Estimated time to finish = %d\n\n', cputime*(imax-i)/Delta)
    end
    SPYMD = abs(theta(:,1:4));
end

%% Parameters Quantiles
quant = [0.05;0.25;0.5;0.75;0.95];
QSPY = quantile([SPYMD,SPYMD(:,3)./SPYMD(:,1)],quant);
SPY_out = rmoutliers(SPYMD,'percentiles',[5,95]);
QSPY_out = quantile([SPY_out,SPY_out(:,3)./SPY_out(:,1)],quant);

for i=1:length(quant)
    fprintf('%3.2f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f \\\\ \n', quant(i), QSPY(i,[1,2,4,3,5]))
end

fprintf('\n')

MDQ = quantize([SPYMD,SPYMD(:,3)./SPYMD(:,1)]);
for i=1:length(quant)
    fprintf('%3.4f & %3.4f & %3.4f & %3.4f & %3.4f \\\\ \n', MDQ(i,[1,2,4,3,5]))
end

%% ADF test on nonstationarity of Delta Drivers
RCDelta = RC(:,1)-RC(:,2);
RCDelta = RCDelta - mean(RCDelta);

fprintf('Augmented Dickey-Fueller Test Rejects Delta GMM Drivers Stationarity: %d\n', adftest(RCDelta))

%% Visualization
vizPath = getPath('Visualization');
d = datetime(d,'ConvertFrom','yyyymmdd'); 
if SY == 8 
    d = d(2:end);
end

figure
hold on
box on
grid on
plot(d,RC(:,1),'-')
plot(d,RC(:,2),'--')
legend('Upper','Lower','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('SPY_RCvsMktReturn_GMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

figure
hold on
box on
grid on
plot(d,RC(:,1)-RC(:,2),'-')
title('Difference between upper and lower GMM drivers','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('SPY_RCvsMktReturn_GMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off



xx = 0.00001:0.00001:0.1;
%c = SPYMD(end-1,1); g = SPYMD(end-1,2); b = SPYMD(end-1,3); a = SPYMD(end-1,4);
c = MDQ(1,1); g = MDQ(1,2); b = MDQ(1,3); a = MDQ(1,4);
Gp = a*(1-exp(-c*xx)).^(1/(1+g));
Gm = (b/c)*(1-exp(-c*xx));

% figure
% hold on
% box on
% grid on
% plot(xx,Gm)
% legend('$\Gamma_-$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% str=strcat('MeasureDistortionsGMMm',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb));
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off
% 
% figure
% hold on
% box on
% grid on
% plot(xx,Gp)
% legend('$\Gamma_+$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% str=strcat('MeasureDistortionsGMMp',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb));
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off

figure
hold on
box on
grid on
plot(xx,Gp,'-')
plot(xx,Gm,'--')
legend('$\Lambda_+$','$\Lambda_-$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex','yscale','log')
str=strcat('MeasureDistortionsGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

corrRC_mu_U = zeros(length(SPYMD)-50,1);
corrRC_mu_L = zeros(length(SPYMD)-50,1);
for i=1:length(SPYMD)-50
    corrRC_mu_U(i) = corr(mu_mkt(i:50+i,2),RC(i:50+i,1));
    corrRC_mu_L(i) = corr(mu_mkt(i:50+i,3),RC(i:50+i,2));
end

figure
hold on
box on
grid on
plot(d(50+1:end),corrRC_mu_U,'-')
plot(d(50+1:end),corrRC_mu_L,'--')
legend('Upper','Lower','interpreter','latex')
str=strcat('CorrGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

quant = [0,0.25,0.5,0.75,1];
CorrU= quantile(corrRC_mu_U,quant);
CorrL= quantile(corrRC_mu_L,quant);
for i=1:length(quant)
    fprintf('%3.2f & %3.2f & %3.2f \\\\ \n', quant(i), CorrU(i), CorrL(i))
end

figure
hold on
box on
grid on
plot(d(50+1:end),corrRC_mu_U-corrRC_mu_L,'-')
title('Difference between upper and lower GMM correlations','interpreter','latex')
str=strcat('CorrGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

%fprintf('corr_u = %d, corr_l = %d\n', corr(mu_mod(:,1),mu_mkt(:,1)), corr(mu_mod(:,3),mu_mkt(:,3)));
%fprintf('corr_u = %d, corr_l = %d\n', corr(RC(:,1),mu_mkt(:,2)), corr(RC(:,2),mu_mkt(:,3)));
fprintf('fraction of days RC negative correlated with mkt returns: %1.2f (upper) and %1.2f (lower)\n', sum(corrRC_mu_U<0)/length(corrRC_mu_U), sum(corrRC_mu_L<0)/length(corrRC_mu_L));
fprintf('Average (mu_u,mu,mu_l):(%d,%d,%d)\n', mean(mu_mkt))

% quant = [0,0.25,0.5,0.75,1];
% CorrSPY= quantile(corrRC_mu_U-corrRC_mu_L,quant);
% for i=1:length(quant)
%     fprintf('%3.2f %3.2f \\\\ \n', quant(i), CorrSPY(i))
% end

%% Save
prompt = 'Do you want to save results? Y/N: ';
s = input(prompt, 's');
%s = 'Y';
if strcmp(s,'Y')
    prompt = 'Warning: data will be saved and previous one overwritten. Input Y to continue: ';
    s = input(prompt, 's');
    if strcmp(s,'Y')
        save(fullfile(varPath,strcat('SPYMD',num2str(SY),'GMM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N))),'SPYMD');
        save(fullfile(varPath,strcat('RC',num2str(SY),'GMM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N))),'RC');
        save(fullfile(varPath,strcat('mu_modGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N))),'mu_mod');
        %save(fullfile(varPath,strcat('mu_mktGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N))),'mu_mkt');
        save(fullfile(varPath,strcat('ErrGMM',num2str(SY),num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N))),'err');
    end
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

function T = objU(u,theta,bp,cp,bn,cn,enforcegamma,enforceb,lb)

    if enforcegamma == 1
        gamma = theta(2);
        range = 1-lb;
        if gamma>1||gamma<lb
            n = floor((gamma-lb)/range);
            if mod(n,2)==1
                gamma = gamma-n*range;
            else
                gamma = 1+n*range-(gamma-lb);
            end
        end
        theta(2) = gamma;
    end
    b = theta(3);
    if enforceb == 1
        b = exp(-b);
    end
    theta(3) = b;
    
    mu = log ( ((1-bp)^(-cp))*((1+bn)^(-cn)) );
    RC = RCU([theta,bp,cp,bn,cn]);
    T = ( norm ( ( u - mu + RC ) * ( u.^([0;1;2;3;4;5;6]) )' ) )^2;
    %T =  ( u - mu + RC ) * ( u - mu + RC )';
end

function T = objL(l,theta,bp,cp,bn,cn,enforcegamma,enforceb,lb)
    if enforcegamma == 1
        gamma = theta(2);
        range = 1-lb;
        if gamma>1||gamma<lb
            n = floor((gamma-lb)/range);
            if mod(n,2)==1
                gamma = gamma-n*range;
            else
                gamma = 1+n*range-(gamma-lb);
            end
        end
        theta(2) = gamma;
    end
    b = theta(3);
    if enforceb == 1
        b = exp(-b);
    end
    theta(3) = b;
    
    mu = log ( ((1-bp)^(-cp))*((1+bn)^(-cn)) );
    RC = RCL([theta,bp,cp,bn,cn]);
    T = norm ( ( ( l - mu + RC ) * ( l.^([0;1]) )' ) )^2;
    %T = ( l - mu - RC ) * ( l - mu - RC )';
end

function [] = plotU(u,theta,N,yp,yn)
    s = u;% + RCU(theta);
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,999,999)/1000;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(5); cp = theta(6); bn = theta(7); cn = theta(8);
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
    f = max(real(fft((1/pi)*exp(1i*u*bb).*phi*eta)),0);
    
    figure
    grid on
    box on
    plot(x,f);
    title('Upper Density')
    set(gca,'TickLabelInterpreter','latex')
    str=strcat('SPY_UpperDensity');
    fname=str;
    %saveas(gcf, fullfile(vizPath, fname), 'epsc');

    f = [0,max(interp1(x,f,s),0)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)];
    Pihat = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
    Pihat(s>0) = 1-Pihat(s>0);
    
    figure
    hold on
    grid on
    box on
    plot(s,max(Pihat,0),'--')
    plot(s,Pi,'-')
    legend('Estimated', 'Empirical','interpreter','latex')
    title('Upper Digital Moments Fitting')
    set(gca,'TickLabelInterpreter','latex')
    str=strcat('SPY_UpperGMMFitting');
    fname=str;
    saveas(gcf, fullfile(vizPath, fname), 'epsc');
end

function [] = plotL(l,theta,N,yp,yn)
    s = l;% - RCL(theta);
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,999,999)/1000;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(5); cp = theta(6); bn = theta(7); cn = theta(8);
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
    f = max(real(fft((1/pi)*exp(1i*u*bb).*phi*eta)),0);
    
    figure
    grid on
    box on
    plot(x,f);
    title('Lower Density')
    set(gca,'TickLabelInterpreter','latex')
    str=strcat('SPY_LowerDensity');
    fname=str;
    %saveas(gcf, fullfile(vizPath, fname), 'epsc');

    f = [0,max(interp1(x,f,s),0)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)];
    Pihat = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
%     figure
%     plot(s,Pihat)
    Pihat(s>0) = 1-Pihat(s>0);
    
    figure
    hold on
    grid on
    box on
    plot(s,max(Pihat,0),'--')
    plot(s,Pi,'-')
    legend('Estimated', 'Empirical','interpreter','latex')
    title('Lower Digital Moments Fitting')
    set(gca,'TickLabelInterpreter','latex')
    str=strcat('SPY_LowerGMMFitting');
    fname=str;
    saveas(gcf, fullfile(vizPath, fname), 'pdf');
end