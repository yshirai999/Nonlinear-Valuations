%% Data Calibration
clear
clc
close all

%% Load data
DataPath = NonlinearPricing.Functions.getPath('Data');
SY = 20;
dat = load(fullfile(DataPath,'ZC_2008_30d.mat'));
R = dat.R;
dat = load(fullfile(DataPath,strcat('SPY_C_T1M_MONEY10_20',num2str(SY,'%02.f'),'.mat')));
C = dat.O;
dat = load(fullfile(DataPath,strcat('SPY_P_T1M_MONEY10_20',num2str(SY,'%02.f'),'.mat')));
P = dat.O;
dates = unique(C(:,2));
ndates = length(dates);
datesNum = datetime(dates,'ConvertFrom','yyyymmdd');

%% FFT parameters
N = 2^12;
eta = 0.15;
lambda = 2*pi/(eta*N);

alpha_C = 1;
alpha_P = 2;

prompt = 'Algo? 0/1/2: ';
algo = str2double(input(prompt,'s'));
while not(algo==0||algo==1||algo==2||algo==3)
    prompt = 'Algo? 0/1/2: ';
    algo = str2double(input(prompt,'s'));
end

prompt = 'Balanced? Y/N: ';
s1 = input(prompt, 's');
if s1 == 'Y'
    s1 = false;
else
    s1 = true;
end
    
prompt = 'Do you want to enforce gamma < 1? Y/N: ';
s2 = input(prompt, 's');
if s2 == 'Y'
    s2 = true;
else
    s2 = false;
end

if s1
    param0 = [0.03,1,0.05,1.2,0.01,0.25,1,1];
    %31 Dec
    %param0 = [3.876152e-03,6.145640e+02,9.791247e-02,3.717545e+00,1.438446e-02,3.404939e-01,1.067310e+00,9.235782e-06]; 
    param0 = [3.876152e-03,6.145640e+02,9.791247e-02,3.717545e+00,12,0.7,0.001,1]; 
else
    param0 = [0.03,1,0.05,1.2,1000,0.1,1000,0.1];
    param0 = [3.876152e-03,6.145640e+02,9.791247e-02,3.717545e+00,1,0.1,1,0.1]; 
    param0 = [3.876152e-03,6.145640e+02,9.791247e-02,3.717545e+00,1.829472e+03,5.892942e-01,2.059276e+03,3.174303e-01];
    %param0 = [0.03,1,0.05,1.2,1000,200,100,150];
    %param0 = [3.876142e-03,6.145676e+02,9.791253e-02,3.717541e+00,3.738839e+03,1.926281e+02,1.200680e+01,7.566878e+02];
end
% param0 = [6.779747e-03,2.877842e+01,5.210701e-02,6.110235e+00,...
%           1.735573e+03,2.427174e+02,6.948442e+01,1.400989e+02]; %Calibrated with Nsim = 100; MaxIter = 350; RC_U = 3.291709e-01; RC_L = 1.233392e-02
% param0 = [6.779747e-03,2.877842e+01,5.210701e-02,6.110235e+00,...
%           1.301409e+03,1.309773e-06,6.977208e+01,2.010158e+07]; %Calibrated with Nsim = 100; MaxIter = 1000; RC_U = 6.178159e-02; RC_L = 1.314772e-02
% param0 = [6.779746e-03,2.877843e+01,5.210701e-02,6.110235e+00,...
%           2.239072e+02,2.291395e-06,2.771185e+01,4.256353e+07]; %Calibrated with Nsim = 5000; MaxIter = 1000; RC_U = 5.048006e-02; RC_L = 3.268094e-02

% param0 = [9.761230e-02,1.267539e+00,5.263158e-02,5.999597e+00,1.121057e+00,9.066959e-02,1.725336e-02,6.942688e-04]; 
% param0 = [0.0039,614.5655,0.0979,3.7175,0.0146,0.2631,0.5104,0.0001];
% param0 = [0.0039,614.5640,0.0979,3.7175,0.0026,0.0026,2.4861,0.0005];
% param0 = [0.0039,614.5640,0.0979,3.7175,0.01,0.2,3.1483,0.01];
% param0 = [0.0039,614.5640,0.0979,3.7175,0.1,0.2,0,0.01];

param0 = [0.0039,614.5672,0.0979,3.7175, 0.0019,1.0328,0.0004,0.0065];
% Montecarlo Integration Parameters
Nsim = 5000;
q = qrandstream('halton',1,'Skip',1e3,'Leap',1e2);
U = qrand(q,Nsim);
Y = - log(U);
E = expint(-log(U));

%% Calibration 

X = zeros(ndates,8);
RC = zeros(ndates,2);
r = zeros(ndates,1);
ErrU = zeros(ndates,2);
ErrL = zeros(ndates,2);
Range = zeros(ndates,1);

optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',5000,'Display','iter');%,'FunctionTolerance',1e-10);
optionsBG = optimset('MaxFunEval',5000,'MaxIter',5000,'Display','iter','TolFun',1e-9,'TolX',1e-9);
optionsD = optimset('MaxFunEval',245,'MaxIter',245,'Display','iter');%,'TolFun',1e-9,'TolX',1e-9);

nn = 0;

for d = ndates
    
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
    
    S = [S_C;S_P]; K = [K_C;K_P]; B = [B_C;B_P]; A = [A_C;A_P]; 
    
    K_0 = S(1);
    
    if length(unique(S))>1
        fprintf('Warning: multiple stock prices for single date.')
    end
    
    f = @(x)(norm(([BGC(abs(x(1)),abs(x(2)),abs(x(3)),abs(x(4)),T,r(d),alpha_C,K_0,K_C,lambda,N,eta);...
                   BGP(abs(x(1)),abs(x(2)),abs(x(3)),abs(x(4)),T,r(d),alpha_P,K_0,K_P,lambda,N,eta)]-(A+B)/2)));
    bp = param0(1); cp = param0(2); bn = param0(3); cn = param0(4); exitflagBG = NaN;
    
    if algo == 0
        [x,~,exitflagBG] = fminsearch(f,[bp,cp,bn,cn],optionsBG); bp = abs(x(1)); cp = abs(x(2)); bn = abs(x(3)); cn = abs(x(4));
    elseif algo == 1
        [x,~,exitflagBG] = fminunc(f,[bp,cp,bn,cn],optionsBG); bp = abs(x(1)); cp = abs(x(2)); bn = abs(x(3)); cn = abs(x(4));
    elseif algo == 2
        [x,~,exitflag] = patternsearch(f,[bp,cp,bn,cn],[],[],[],[],[],[],[],optionspsrch);
    end
    
    C_Mod_M = BGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta);
    P_Mod_M = BGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta);
    
    if s1
        f = @(x)(norm( ([UBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,abs(x(1)),x(2),abs(x(3)),abs(x(4)),Nsim,Y,E,s2)-...
                         LBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,abs(x(1)),x(2),abs(x(3)),abs(x(4)),Nsim,Y,E,s2);...
                         UBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,abs(x(1)),x(2),abs(x(3)),abs(x(4)),Nsim,Y,E,s2)-...
                         LBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,abs(x(1)),x(2),abs(x(3)),abs(x(4)),Nsim,Y,E,s2)]...
                         -(A-B) ) ./ (A+B) ) );
        c = param0(5); gamma = param0(6); a = param0(7); b = param0(8); exitflagD = NaN;
        if algo == 0
            [x,~,exitflagD] = fminsearch(f,[c,gamma,a,b],optionsD); c = abs(x(1)); gamma = x(2); a = abs(x(3)); b = abs(x(4));
        elseif algo == 1
            [x,~,exitflagD] = fminunc(f,[c,gamma,a,b],optionsD); c = abs(x(1)); gamma = x(2); a = abs(x(3)); b = abs(x(4));
        elseif algo == 2
            [x,~,exitflag] = patternsearch(f,[c,gamma,a,b],[],[],[],[],[],[],[],optionspsrch); c = abs(x(1)); gamma = x(2); a = abs(x(3)); b = abs(x(4));
        else
            c = param0(5); gamma = param0(6); a = param0(7); b = param0(8);
        end
        if s2
            range = 1;
            if gamma>1||gamma<0
                n = floor(gamma);
                if mod(n,2)==1
                    gamma = gamma-n;
                else
                    gamma = 1+n-gamma;
                end
            end
        end
        C_Mod_U = UBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2);
        C_Mod_L = LBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2);
    
        P_Mod_L = LBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2);
        P_Mod_U = UBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2);
        
        RC(d,:) = [RCU([c,gamma,a,b,bp,cp,bn,cn]),RCL([c,gamma,b,a,bp,cp,bn,cn])];
        
        X(d,:) = [bp,cp,bn,cn,c,gamma,a,b];
    else
        f = @(x)(norm( ([UBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,abs(x(1)),x(2),1/abs(x(1)),1,Nsim,Y,E,s2)-...
                         LBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,abs(x(3)),x(4),1/abs(x(3)),1,Nsim,Y,E,s2);...
                         UBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,abs(x(3)),x(4),1/abs(x(3)),1,Nsim,Y,E,s2)-...
                         LBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,abs(x(1)),x(2),1/abs(x(1)),1,Nsim,Y,E,s2)]./...
                         [C_Mod_M;P_Mod_M])-(A-B)./((A+B)/2) ) );
        cU = param0(5); gammaU = param0(6);
        cL = param0(7); gammaL = param0(8); exitflagD = NaN;
    
        if algo == 0
            [x,~,exitflagD] = fminsearch(f,[cU,gammaU,cL,gammaL],optionsD); 
        elseif algo == 1
            [x,~,exitflagD] = fminunc(f,[cU,gammaU,cL,gammaL],optionsD);
        elseif algo == 2
            [x,~,exitflag] = patternsearch(f,[cU,gammaU,cL,gammaL],[],[],[],[],[],[],[],optionspsrch);
        end
        
        cU = abs(x(1)); gammaU = x(2); aU = 1/cU; bU = 1;
        cL = abs(x(3)); gammaL = x(4); aL = 1/cL; bL = 1;
        if s2
            if gammaU>1||gammaU<0
                n = floor(gammaU);
                if mod(n,2)==1
                    gammaU = gammaU-n;
                else
                    gammaU = 1+n-gammaU;
                end
            end
            if gammaL>1||gammaL<0
                n = floor(gammaL);
                if mod(n,2)==1
                    gammaL = gammaL-n;
                else
                    gammaL = 1+n-gammaL;
                end
            end
        else
            gammaU = abs(gammaU);
            gammaL = abs(gammaL);
        end    
    
        C_Mod_U = UBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,cU,gammaU,aU,bU,Nsim,Y,E,s2);
        C_Mod_L = LBGC(bp,cp,bn,cn,T,r(d),alpha_C,K_0,K_C,lambda,N,eta,cL,gammaL,aL,bL,Nsim,Y,E,s2);
    
        P_Mod_L = LBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,cU,gammaU,aU,bU,Nsim,Y,E,s2);
        P_Mod_U = UBGP(bp,cp,bn,cn,T,r(d),alpha_P,K_0,K_P,lambda,N,eta,cL,gammaL,aL,bL,Nsim,Y,E,s2);
        
        RC(d,:) = [RCU([cU,gammaU,aU,bU,bp,cp,bn,cn]),RCL([cL,gammaL,bL,aL,bp,cp,bn,cn])];
        
        X(d,:) = [bp,cp,bn,cn,cU,gammaU,cL,gammaL];
    end
    

    %[[C_Mod_U;P_Mod_U], A, [C_Mod_M;P_Mod_M], (A+B)/2, [C_Mod_L;P_Mod_L], B, [C_Mod_U;P_Mod_U]-[C_Mod_L;P_Mod_L], A-B];
    
    ErrU(d,:) = [max(abs(C_Mod_U-A_C)),max(abs(P_Mod_U-A_P))];
    ErrL(d,:) = [max(abs(C_Mod_L-B_C)),...
               max(abs(P_Mod_L-B_P))];
    Range(d) = length(A);
    
    %param0 = [bp,cp,bn,cn,c,gamma,a,b];
    
    if s1
        fprintf('\nStep = %d: (bp,cp,bn,cn) = (%3.4f,%3.4f,%3.4f,%3.4f); (c,gamma,a,b) = (%3.4f,%3.4f,%3.4f,%3.4f);\n',...
            d, bp, cp, bn, cn, c, gamma, a, b)
    else
        fprintf('\nStep = %d: (bp,cp,bn,cn) = (%3.4f,%3.4f,%3.4f,%3.4f); (cU,gammaU,aU,bU) = (%3.4f,%3.4f,%3.4f,%3.4f); (cL,gammaL,aL,bL) = (%3.4f,%3.4f,%3.4f,%3.4f);\n',...
            d, bp, cp, bn, cn, cU, gammaU, aU, bU, cL, gammaL, aL, bL)
    end
    fprintf('Max Error U = %d; Max Error L = %d; \n', max(ErrU(d,2)), max(ErrL(d,2)))
    fprintf('# of Options Calibrated to = %d; flagBG = %d; flagD = %d \n', Range(d), exitflagBG, exitflagD)
    fprintf('RC_U = %d; RC_L = %d \n', RC(d,1),RC(d,2))
    fprintf('cputime = %d\n\n', cputime)
    
end

if s1
    fprintf('param0 = [%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f];\n', bp,cp,bn,cn,c,gamma,a,b);
else
    fprintf('param0 = [%d,%d,%d,%d,%d,%d,%d,%d];\n', bp,cp,bn,cn,cU,gammaU,cL,gammaL);
end

%% Visualization
vizPath = NonlinearPricing.Functions.getPath('Visualization');

figure
hold on
grid on
box on
scatter([K_P;K_C]/K_0-1,([P_Mod_U;C_Mod_U]-[P_Mod_L;C_Mod_L])./[P_Mod_M;C_Mod_M],'*') %plot(K/K_0-1,[C_Mod_U;P_Mod_U]-[C_Mod_L;P_Mod_L])
scatter([K_P;K_C]/K_0-1,([A_P;A_C]-[B_P;B_C])./[(A_P+B_P)/2;(A_C+B_C)/2],'o'); %plot(K/K_0-1,A-B)
legend('Model Implied Bid-Ask spread','Observed Bid-Ask spread','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('BidAskSpread_GcgammaDouble',num2str(s1),num2str(s2));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off
        
[callU{1:length(K_C)}] = deal('Upper');
[callM{1:length(K_C)}] = deal('Mid');
[callL{1:length(K_C)}] = deal('Lower');
callU = categorical(callU);
callM = categorical(callM);
callL = categorical(callL);

figure
hold on
grid on
box on
scatter3(K_C/K_0-1,callU,A_C,'o','red')
scatter3(K_C/K_0-1,callM,(A_C+B_C)/2,'o','blue')
scatter3(K_C/K_0-1,callL,B_C,'o','black')
scatter3(K_C/K_0-1,callU,C_Mod_U,'*','red')
scatter3(K_C/K_0-1,callM,C_Mod_M,'*','blue')
scatter3(K_C/K_0-1,callL,C_Mod_L,'*','black')
view(171,18)
%xlabel('Moneyness')
set(gca,'TickLabelInterpreter','latex')
str=strcat('CallScatter_GcgammaDouble',num2str(s1),num2str(s2));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

[putU{1:length(K_P)}] = deal('Upper');
[putM{1:length(K_P)}] = deal('Mid');
[putL{1:length(K_P)}] = deal('Lower');
putU = categorical(putU);
putM = categorical(putM);
putL = categorical(putL);

figure
hold on
grid on
box on
scatter3(K_C/K_0-1,callU,A_C,'o','red')
scatter3(K_C/K_0-1,callM,(A_C+B_C)/2,'o','blue')
scatter3(K_C/K_0-1,callL,B_C,'o','black')
scatter3(K_C/K_0-1,callU,C_Mod_U,'*','red')
scatter3(K_C/K_0-1,callM,C_Mod_M,'*','blue')
scatter3(K_C/K_0-1,callL,C_Mod_L,'*','black')
scatter3(K_P/K_0-1,putU,A_P,'o','red')
scatter3(K_P/K_0-1,putM,(A_P+B_P)/2,'o','blue')
scatter3(K_P/K_0-1,putL,B_P,'o','black')
scatter3(K_P/K_0-1,putU,P_Mod_U,'*','red')
scatter3(K_P/K_0-1,putM,P_Mod_M,'*','blue')
scatter3(K_P/K_0-1,putL,P_Mod_L,'*','black')
view(189,9)
%xlabel('Moneyness')
set(gca,'TickLabelInterpreter','latex')
str=strcat('Scatter_GcgammaDouble',num2str(s1),num2str(s2));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc'); % Updated to use vizPath
hold off    

figure
hold on
grid on
box on
scatter3(K_P/K_0-1,putU,A_P,'o','red')
scatter3(K_P/K_0-1,putM,(A_P+B_P)/2,'o','blue')
scatter3(K_P/K_0-1,putL,B_P,'o','black')
scatter3(K_P/K_0-1,putU,P_Mod_U,'*','red')
scatter3(K_P/K_0-1,putM,P_Mod_M,'*','blue')
scatter3(K_P/K_0-1,putL,P_Mod_L,'*','black')
view(189,9)
%xlabel('Moneyness')
set(gca,'TickLabelInterpreter','latex')
str=strcat('PutScatter_GcgammaDouble',num2str(s1),num2str(s2));
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc'); % Updated to use vizPath
hold off    

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
    str=strcat('CalibrationError_GcgammaDouble',num2str(s1),num2str(s2));
    fname=str;
    saveas(gcf, fullfile(vizPath, fname), 'epsc'); % Updated to use vizPath
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
    str=strcat('RC_RiskNeutral_GcgammaDouble',num2str(s1),num2str(s2));
    fname=str;
    saveas(gcf, fullfile(vizPath, fname), 'epsc'); % Updated to use vizPath
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
    varArchivePath = NonlinearPricing.Functions.getVarArchivePath();
    save(fullfile(varArchivePath, strcat('BG_GcgammaDouble',num2str(s1),num2str(s2))),'X')
    save(fullfile(varArchivePath, strcat('RC_GcgammaDouble',num2str(s1),num2str(s2))),'RC')
    save(fullfile(varArchivePath, strcat('ErrU_GcgammaDouble',num2str(s1),num2str(s2))),'ErrU')
    save(fullfile(varArchivePath, strcat('ErrL_GcgammaDouble',num2str(s1),num2str(s2))),'ErrL')
    save(fullfile(varArchivePath, strcat('Range_GcgammaDouble',num2str(s1),num2str(s2))),'Range')
    save(fullfile(varArchivePath, strcat('rf_20',num2str(SY,'%02.f'),'_30d',num2str(s1),num2str(s2))),'r')
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

function Phi = UBGPhi(u,bp,cp,bn,cn,T,K_0,r,c,gamma,a,b,Nsim,Y,E)
yp = Y*bp;
yn = Y*bn;
d0 = 1i*u*(log(K_0)+T*(r+cp*log(1-bp)+cn*log(1+bn)));
Phi = exp(d0) .* ((1-1i*u*bp).^(-T*cp)).*((1+1i*u*bn).^(-T*cn)) .* exp( T* sum(...
            (exp(1i*yp*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cp*E)).^(-gamma/(1+gamma)) )...
                .*  bp*cp .* exp(-c*cp*E) ./ yp...
            -(exp(-1i*yn*u)-1) .*  b*bn*cn .* exp(-c*cn*E) ./ yn) /Nsim);
end

function Psi = UBGPsi(u,bp,cp,bn,cn,T,r,alpha,K_0,c,gamma,a,b,Nsim,Y,E)
Psi = (exp(-r*T)./((alpha+1i*u) .* (alpha+1i*u+1))) .* UBGPhi(u-(alpha+1)*1i,bp,cp,bn,cn,T,K_0,r,c,gamma,a,b,Nsim,Y,E);
end

function Phi = LBGPhi(u,bp,cp,bn,cn,T,K_0,r,c,gamma,a,b,Nsim,Y,E)
yp = Y*bp;
yn = Y*bn;
d0 = 1i*u*(log(K_0)+T*(r+cp*log(1-bp)+cn*log(1+bn)));
Phi = exp(d0) .* ((1-1i*u*bp).^(-T*cp)).*((1+1i*u*bn).^(-T*cn)) .* exp( T * sum(...
            -(exp(1i*yp*u)-1) .*  b*bp*cp .* exp(-c*cp*E) ./ yp...
            +(exp(-1i*yn*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cn*E)).^(-gamma/(1+gamma)) )...
                .*  bn*cn .* exp(-c*cn*E) ./ yn) /Nsim);
end

function Psi = LBGPsi(u,bp,cp,bn,cn,T,r,alpha,K_0,c,gamma,a,b,Nsim,Y,E)
Psi = (exp(-r*T)./((alpha+1i*u) .* (alpha+1i*u+1))) .* LBGPhi(u-(alpha+1)*1i,bp,cp,bn,cn,T,K_0,r,c,gamma,a,b,Nsim,Y,E);
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

function C = UBGC(bp,cp,bn,cn,T,r,alpha,K_0,K,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2)
if s2
    if gamma>1||gamma<0
        n = floor(gamma);
        if mod(n,2)==1
            gamma = gamma-n;
        else
            gamma = 1+n-gamma;
        end
    end
else
    gamma = abs(gamma);
end
I = cumsum([0,ones(1,N-1)]);
beta = log(K_0)-lambda*N/2;
k = beta+I*lambda;
u = I*eta;
w = ones(1,N)*eta; w(1)=w(1)/2;
x = exp(-1i*beta*u).*UBGPsi(u,bp,cp,bn,cn,T,r,alpha,K_0,c,gamma,a,b,Nsim,Y,E).*w;
Call = (exp(-alpha*k)/pi).*max(real(fft(x,N)),0);
kk = log(K);
C = interp1(k,Call,kk);
end

function C = LBGC(bp,cp,bn,cn,T,r,alpha,K_0,K,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2)
if s2
    if gamma>1||gamma<0
        n = floor(gamma);
        if mod(n,2)==1
            gamma = gamma-n;
        else
            gamma = 1+n-gamma;
        end
    end
else
    gamma = abs(gamma);
end
I = cumsum([0,ones(1,N-1)]);
beta = log(K_0)-lambda*N/2;
k = beta+I*lambda;
u = I*eta;
w = ones(1,N)*eta; w(1)=w(1)/2;
x = exp(-1i*beta*u).*LBGPsi(u,bp,cp,bn,cn,T,r,alpha,K_0,c,gamma,a,b,Nsim,Y,E).*w;
Call = (exp(-alpha*k)/pi).*max(real(fft(x,N)),0);
kk = log(K);
C = interp1(k,Call,kk);
end

function P = UBGP(bp,cp,bn,cn,T,r,alpha,K_0,K,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2)
if s2
    if gamma>1||gamma<0
        n = floor(gamma);
        if mod(n,2)==1
            gamma = gamma-n;
        else
            gamma = 1+n-gamma;
        end
    end
else
    gamma = abs(gamma);
end
I = cumsum([0,ones(1,N-1)]);
beta = log(K_0)-lambda*N/2;
k = beta+I*lambda;
u = I*eta;
w = ones(1,N)*eta; w(1)=w(1)/2;
x = exp(-1i*beta*u).*LBGPsi(u,bp,cp,bn,cn,T,r,-alpha,K_0,c,gamma,a,b,Nsim,Y,E).*w;
Put = (exp(alpha*k)/pi).*max(real(fft(x,N)),0);
kk = log(K);
P = interp1(k,Put,kk);
end

function P = LBGP(bp,cp,bn,cn,T,r,alpha,K_0,K,lambda,N,eta,c,gamma,a,b,Nsim,Y,E,s2)
if s2
    if gamma>1||gamma<0
        n = floor(gamma);
        if mod(n,2)==1
            gamma = gamma-n;
        else
            gamma = 1+n-gamma;
        end
    end
else
    gamma = abs(gamma);
end
I = cumsum([0,ones(1,N-1)]);
beta = log(K_0)-lambda*N/2;
k = beta+I*lambda;
u = I*eta;
w = ones(1,N)*eta; w(1)=w(1)/2;
x = exp(-1i*beta*u).*UBGPsi(u,bp,cp,bn,cn,T,r,-alpha,K_0,c,gamma,a,b,Nsim,Y,E).*w;
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