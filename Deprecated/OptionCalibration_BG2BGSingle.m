%% Data Calibration
clear
clc
close all

%% Load data
dat = load('ZC_2008_30d.mat');
R = dat.R;
dat = load('SPY_C_T1M_MONEY10_2020.mat');
C = dat.O;
dat = load('SPY_P_T1M_MONEY10_2020.mat');
P = dat.O;
dates = unique(C(:,2));
ndates = length(dates);
datesNum = datetime(dates,'ConvertFrom','yyyymmdd');

%% FFT parameters
N = 2^14;
eta = 0.15;
lambda = 2*pi/(eta*N);

alpha_C = 1;
alpha_P = 20;

param0 = [2.095450e-02,1.546497e+01,2.090907e-02,1.546497e+01];

%% Calibration 

X = zeros(ndates,4);
RC = zeros(ndates,2);
r = zeros(ndates,1);
ErrU = zeros(ndates,2);
ErrL = zeros(ndates,2);
Range = zeros(ndates,1);

optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',5000,'Display','iter');%,'FunctionTolerance',1e-10);
optionsBG = optimset('MaxFunEval',5000,'MaxIter',5000,'Display','iter','TolFun',1e-9,'TolX',1e-9);
optionsD = optimset('MaxFunEval',1000,'MaxIter',1000,'Display','iter','TolFun',1e-9,'TolX',1e-9);

nn = 0;

for d = 1%:ndates
    
    ind_R = (R(:,1) == dates(d));
    
    T = 30/360; r(d) = interp1(R(ind_R,2:3)/360,R(ind_R,4:5),T)/100;
    
    if isnan(r(d))
        fprintf('Warning: risk free rate is NAN')
    end
    
    ind_C = (C(:,2) == dates(d));
    S_C = C(ind_C,4); K_C = C(ind_C,5); B_C = C(ind_C,6); A_C = C(ind_C,7);
    ind_C = (K_C>S_C(1)); %if >, only OTM calls for calibration
    S_C = S_C(ind_C); K_C = K_C(ind_C); B_C = B_C(ind_C); A_C = A_C(ind_C);
    
    ind_P = (P(:,2) == dates(d));
    S_P = P(ind_P,4); K_P = C(ind_P,5); B_P = P(ind_P,6); A_P = P(ind_P,7);
    ind_P = (K_P<S_P(1)); %if <, only OTM puts for calibration
    S_P = S_P(ind_P); K_P = K_P(ind_P); B_P = B_P(ind_P); A_P = A_P(ind_P);
    
    S = [S_C;S_P]; K = [K_C;K_P]; L = [B_C;A_P]; U = [A_C;B_P]; 
    
    K_0 = S(1);
    
    if length(unique(S))>1
        fprintf('Warning: multiple stock prices for single date.')
    end
    
    f = @(x)(norm([BGC(abs(x(1)),abs(x(2)),abs(x(3)),abs(x(2)),T,r(d),alpha_C,K_0,K_C,lambda,N,eta);...
                   BGP(abs(x(1)),abs(x(2)),abs(x(3)),abs(x(2)),T,r(d),alpha_P,K_0,K_P,lambda,N,eta);
                   BGC(abs(x(3)),abs(x(2)),abs(x(1)),abs(x(2)),T,r(d),alpha_C,K_0,K_C,lambda,N,eta);
                   BGP(abs(x(3)),abs(x(2)),abs(x(1)),abs(x(2)),T,r(d),alpha_P,K_0,K_P,lambda,N,eta);]-[U;L]));
    bp = param0(1); cp = param0(2); bn = param0(3); cn = param0(2); exitflagBG = NaN;
    
    [x,~,exitflagBG] = fminsearch(f,[bp,cp,bn,cn],optionsBG); bp = abs(x(1)); cp = abs(x(2)); bn = abs(x(3)); cn = abs(x(2));
    
    tic;
    C_Mod_U = BGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta);
    C_Mod_L = BGC(bn,cp,bp,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta);
    
    P_Mod_L = BGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta);
    P_Mod_U = BGP(bn,cp,bp,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta);
    cputime = toc;

    [[C_Mod_U;P_Mod_U], [A_C;A_P], [C_Mod_L;P_Mod_L], [B_C;B_P], [C_Mod_U;P_Mod_U]-[C_Mod_L;P_Mod_L], [A_C;A_P]-[B_C;B_P]]
    
    ErrU(d,:) = [max(abs(C_Mod_U-A_C)),max(abs(P_Mod_U-A_P))];
    ErrL(d,:) = [max(abs(C_Mod_L-B_C)),max(abs(P_Mod_L-B_P))];
    Range(d) = length(U);
    
    %RC(d,:) = [RCU([c,gamma,b,a,bp,cp,bn,cn]),RCL([c,gamma,b,a,bp,cp,bn,cn])];
    
    X(d,:) = [bp,cp,bn,cn];
    
    %param0 = [bp,cp,bn,cn,c,gamma,a,b];
    
    fprintf('\nStep = %d: (bp,cp,bn,cn) = (%d,%d,%d,%d);\n', d, bp, cp, bn, cn)
    fprintf('Max Error U = %d; Max Error L = %d;\n', max(ErrU(d,2)), max(ErrL(d,2)))
    fprintf('# of Options Calibrated to = %d; flagBG = %d;\n', Range(d), exitflagBG)
    %fprintf('RC_U = %d; RC_L = %d \n', RC(d,1),RC(d,2))
    fprintf('cputime = %d\n\n', cputime)
    
    if d == 1
        figure
        hold on
        grid on
        box on
        plot([K_P;K_C]/K_0-1,([P_Mod_U;C_Mod_U]-[P_Mod_L;C_Mod_L])./[(A_P+B_P)/2;(A_C+B_C)/2])
        plot([K_P;K_C]/K_0-1,([A_P;A_C]-[B_P;B_C])./[(A_P+B_P)/2;(A_C+B_C)/2]);
        legend('Model Implied Bid-Ask spread','Observed Bid-Ask spread','Interpreter','latex')
        %xlabel('Moneyness')
        set(gca,'TickLabelInterpreter','latex')
        fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
        str=strcat('BidAskSpread_BG2BGSingle');
        fname=str;
        saveas(gcf, fullfile(fpath, fname), 'epsc');
        hold off
        
        [callU{1:length(K_C)}] = deal('Upper');
        [callL{1:length(K_C)}] = deal('Lower');
        callU = categorical(callU);
        callL = categorical(callL);

        figure
        hold on
        grid on
        box on
        scatter3(K_C/K_0-1,callU,A_C,'o','red')
        scatter3(K_C/K_0-1,callL,B_C,'o','black')
        scatter3(K_C/K_0-1,callU,C_Mod_U,'*','red')
        scatter3(K_C/K_0-1,callL,C_Mod_L,'*','black')
        view(171,18)
        xlabel('Moneyness')
        set(gca,'TickLabelInterpreter','latex')
        fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
        str=strcat('CallScatter_BG2BGSingle');
        fname=str;
        saveas(gcf, fullfile(fpath, fname), 'epsc');
        hold off
        
        [putU{1:length(K_P)}] = deal('Upper');
        [putL{1:length(K_P)}] = deal('Lower');
        putU = categorical(putU);
        putL = categorical(putL);

        figure
        hold on
        grid on
        box on
        scatter3(K_P/K_0-1,putU,A_P,'o','red')
        scatter3(K_P/K_0-1,putL,B_P,'o','black')
        scatter3(K_P/K_0-1,putU,P_Mod_U,'*','red')
        scatter3(K_P/K_0-1,putL,P_Mod_L,'*','black')
        view(171,18)
        xlabel('Moneyness')
        set(gca,'TickLabelInterpreter','latex')
        fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
        str=strcat('PutScatter_BG2BGSingle');
        fname=str;
        saveas(gcf, fullfile(fpath, fname), 'epsc');
        hold off
    end
    
    
end

%% Visualization

prompt = 'Do you want to visualize results? Y/N: ';
s = input(prompt, 's');
if s == 'Y'
    prompt = 'Warning: plots will be saved and previous ones overwritten. Input Y to continue: ';
    s = input(prompt, 's');
    if s == 'Y'
        s = true;
    else
        s = false;
    end
else
    s = false;
end
    
if s
    close all
    
    ErrMax = [max(ErrU,[],2),max(ErrL,[],2)];
    
    [ErrMaxout,ErrMaxoutind] = rmoutliers(ErrMax,'percentiles',[0,80]);
    figure
    hold on
    grid on
    box on
    plot(datesNum(~ErrMaxoutind),ErrMaxout(:,1))
    plot(datesNum(~ErrMaxoutind),ErrMaxout(:,2))
    legend('Upper Valuation Error','Lower Valuation Error','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
    str=strcat('CalibrationError_BG2BGSingle');
    fname=str;
    saveas(gcf, fullfile(fpath, fname), 'epsc');
    hold off
    
    fprintf('\nMedian Errors: Upper Call = %d, Upper Put = %d, Lower Call = %d, Lower Put = %d\n\n',...
            median(ErrU(:,1)),median(ErrU(:,2)),median(ErrL(:,2)),median(ErrL(:,2)))
    
    RCErr = RC(~ErrMaxoutind,:);
    DErr = datesNum(~ErrMaxoutind);
    [RCErrout,RCoutind] = rmoutliers(RCErr,'percentiles',[5,95]);
    figure
    hold on
    grid on
    box on
    plot(DErr(~RCoutind),RCErrout(:,1))
    plot(DErr(~RCoutind),RCErrout(:,2))
    plot(datesNum,r)
    legend('$RC^U$','$RC^L$','$r$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
    str=strcat('RC_RiskNeutral_BG2BGSingle');
    fname=str;
    saveas(gcf, fullfile(fpath, fname), 'epsc');
    hold off
end

%% save

prompt = 'Do you want to save results? Y/N: ';
s = input(prompt, 's');
if s == 'Y'
    prompt = 'Warning: previous results will be overwritten. Input Y to continue: ';
    s = input(prompt, 's');
    if s == 'Y'
        s = true;
    else
        fprintf('Results will not be saved\n')
        s = false;
    end
else
    fprintf('Results will not be saved\n')
    s = false;
end

if s
    save('BG_BG2BGSingle','X')
    %save('RC_BG2BGSingle','RC')
    save('ErrU_BG2BGSingle','ErrU')
    save('ErrL_BG2BGSingle','ErrL')
    save('Range_BG2BGSingle','Range')
    save('rf_2008_30d','r')
    fprintf('Results were saved, and previous ones overwritten')
end

%% Routines

function Phi = BGPhi(u,bp,cp,bn,cn,T,K_0,r)
d0 = 1i*u*(log(K_0)+T*(r+cp*log(1-bp)+cn*log(1+bn)));
Phi = exp(d0) .* ((1-1i*u*bp).^(-T*cp)).*((1+1i*u*bn).^(-T*cn));
end

function Psi = BGPsi(u,bp,cp,bn,cn,T,r,alpha,K_0)
Psi = (exp(-r*T)./((alpha+1i*u) .* (alpha+1i*u+1))) .* BGPhi(u-(alpha+1)*1i,bp,cp,bn,cn,T,K_0,r);
end

function C = BGC(bp,cp,bn,cn,T,r,alpha,K_0,K,lambda,N,eta)
I = cumsum([0,ones(1,N-1)]);
beta = log(K_0)-lambda*N/2;
k = beta+I*lambda;
u = I*eta;
w = ones(1,N)*eta; w(1)=w(1)/2;
x = exp(-1i*beta*u).*BGPsi(u,bp,cp,bn,cn,T,r,alpha,K_0).*w;
Call = (exp(-alpha*k)/pi).*max(real(fft(x,N)),0);
kk = log(K);
C = interp1(k,Call,kk);
end

function P = BGP(bp,cp,bn,cn,T,r,alpha,K_0,K,lambda,N,eta)
I = cumsum([0,ones(1,N-1)]);
beta = log(K_0)-lambda*N/2;
k = beta+I*lambda;
u = I*eta;
w = ones(1,N)*eta; w(1)=w(1)/2;
x = exp(-1i*beta*u).*BGPsi(u,bp,cp,bn,cn,T,r,-alpha,K_0).*w;
Put = (exp(alpha*k)/pi).*max(real(fft(x,N)),0);
kk = log(K);
P = interp1(k,Put,kk);
end

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