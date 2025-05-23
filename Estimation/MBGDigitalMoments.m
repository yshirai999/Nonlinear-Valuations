%% Digital Moments MBG
clear
clc
close all

%% Load Data
ticker = {'SPY', 'VIX', 'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY'};
ETF = matlab.lang.makeValidName(ticker);
D = length(ticker);

for d=1:D

Y = load(strcat('Y',ticker{d},'.mat'));
Y = Y.Y;
ind = (Y(:,1)>20151231);
eval([ETF{d} ' = Y(ind,2:end);']);

end

%SPY(:,1:2) = SPY(:,2:-1:1);

dates = datetime(Y(ind,1),'ConvertFrom','yyyymmdd');

N = length(dates);

lowvec = zeros(N,D);
highvec = zeros(N,D);
midvec = zeros(N,D);
bpvec = zeros(N,D);
cpvec = zeros(N,D);
bnvec = zeros(N,D);
cnvec = zeros(N,D);
Y = zeros(N,D);

for d = 1:D
    eval(['lowvec(:,d) = ' ETF{d} '(:,1);']);
    eval(['highvec(:,d) = ' ETF{d} '(:,2);']);
    eval(['midvec(:,d) = ' ETF{d} '(:,3);']);
    eval(['bpvec(:,d) = ' ETF{d} '(:,4);']);
    eval(['cpvec(:,d) = ' ETF{d} '(:,5);']);
    eval(['bnvec(:,d) = ' ETF{d} '(:,6);']);
    eval(['cnvec(:,d) = ' ETF{d} '(:,7);']);
    eval(['Y(:,d) = ' ETF{d} '(:,2);']);
end

vis = 1;

options = optimset('MaxFunEval',10000,'MaxIter',10000,'Display','iter','TolFun',1e-10,'TolX',1e-10);
optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',25,'Display','off','FunctionTolerance',1e-10);

%% Sample Generation

varPath = getPath('VarArchive');

Delta = 252; %lookback period for each estimate
Nmin = Delta; %date at which estimation can start
kmax = 252; %sample size at each date

try
    x = load(fullfile(varPath, strcat('xMBGSample',num2str(kmax))));
    x = x.x;
    TT = load(fullfile(varPath, strcat('TTMBGTailProb',num2str(kmax))));
    TT = TT.TT;
catch
    x = zeros(kmax,D); %sample to be generated
    TT = zeros(N-Nmin,kmax); %observed tail for each date and x in sample
    
    for n=Nmin+1:N
        k = 1;
        s = log(midvec(n-Nmin+1:n,:)./midvec(n-Nmin:n-1,:));
        alpha = max(-min(min(s)),max(max(s))); %max abs jump size in the sample to be generated
        while k<=kmax
            j = randperm(D,1); ij = randperm(D,j); xij = 2*alpha*rand(1,j)-alpha;
            x(k,ij,n-Nmin) = xij;
            TT(n-Nmin,k) = T(x(k,:,n-Nmin),s);
            if TT(n-Nmin,k)>0.025
                k=k+1;
                %fprintf('k=%d, T_k = %d, j = %d\n', k-1, TT(n-Nmin,k-1), j)
            else
                x(k,ij,n-Nmin) = zeros(1,j);
            end
        end
        fprintf('n=%d\n',n-Nmin+1)
    end

    save(fullfile(varPath, strcat('xMBGSample',num2str(kmax))), 'x')
    save(fullfile(varPath, strcat('TTMBGTailProb',num2str(kmax))), 'TT')
end

%% Tails 
M = 100;
rng('default'); %freeze the seed
rng(1);
Z = randn(M,D);
eps = 1e-10;
g = zeros(N-Nmin,D*(D-1)/2);
zeta = zeros(N-Nmin,1);

g0 = [60.3465,651.3773,-0.0038,-205.6770,4.8418,24.2426,-0.0664,0.0613,0.0595,-290.6120,161.5743,...
      -0.0038,0.0027,421.5420,535.0741,-0.0009,0.0145,0.0159,0.0747,0.0021,-0.0016,-0.0722,...
      315.5769,-0.0080,0.0150,-0.0045,-0.0013,-0.2185,-0.0011,-0.0009,-0.0043,-0.0007,-0.0084,...
      0.0105,-0.0185,-0.0012,282.6299,0.0075,0.0007,0.0169,-0.0043,-6.9368,0.1568,0.0002,...
      0.0024,-594.3617,-5.4250,-0.0153,-0.0027,253.6199,0.0012,0.0096,0.0323,0.0007,-0.1906];
zeta0 = 1.8209;

for n = Nmin+1:Nmin+1
    bp = bpvec(n,:);
    cp = cpvec(n,:);
    bn = bnvec(n,:);
    cn = cnvec(n,:);
    
    zeta_min = 1/(min(min([cp;cn])));
    
%     g0 = zeros(1,D*(D-1)/2);
%     zeta0 = zeta_min;
    
    f = @(param)obj(param,zeta_min,Z,eps,bp,cp,bn,cn,x(:,:,n-Nmin),TT(n-Nmin,:)');
    
    [param] = patternsearch(f,[g0,zeta0],[],[],[],[],[],[],[],optionspsrch);
%     [param] = fminsearch(f,[g0,zeta0],options);
    g(n-Nmin,:) = param(1:end-1);
    zeta(n-Nmin) = param(end);
    zeta = max(zeta,zeta_min);
    
end

%%
C = C_AH(param(1:end-1),D,eps); C = C*C;
fprintf('g0 = [%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f];\nzeta0 = %.4f;\n', param)
%% Visualization

vizPath = getPath('Visualization');
% figure
% hold on
% box on
% grid on
% plot(d,RC(:,1))
% plot(d,RC(:,2))
% plot(d,mu_mod(:,2))
% legend('$RC^U$','$RC^L$','$\mu$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% str=strcat('SPY_RCvsMktReturn_DigitalMoments');
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off



%% Routines

function Tstat = T(x,s)
    Tstat = 0;
    [N,~] = size(s);
    for n=1:N
        tp = (s(n,:)>x).*(s(n,:)>0);
        tn = (s(n,:)<x).*(s(n,:)<0);
        Tstat = Tstat+prod(tp+tn)/N;
    end
end

function [C,A] = C_AH(g,D,eps)
%return SQRT of unique correlation matrix exp(A[x]) for A symmetric.
    x = ones(D,1);
    A = diag(x,0);
    imax = 0;
    for i=1:D-1 
        imin = imax+1;
        imax = imax+D-i;
        A = A + diag(g(imin:imax),-i)+diag(g(imin:imax),i);
    end
    xnew = x-log(diag(expm(A)));
    A = A-diag(x,0)+diag(xnew,0);
    while norm(xnew-x,2)>eps
        x = xnew;
        xnew = x-log(diag(expm(A)));
        A = A-diag(x,0)+diag(xnew,0);
    end
    C = expm(A/2); %WARNING: this is the SQRT of the correlation matrix!
end

function T = hatT(g,zeta,Z,eps,bp,cp,bn,cn,x)
    [M,D] = size(Z);
    rng('default'); %freeze the seed
    rng(1);
    gam = gamrnd(1/zeta,zeta,M,1);
    sigma = sqrt(2*bp.*bn/zeta);
    vartheta = (bp-bn)/zeta;
    C = C_AH(g,D,eps);
    S = diag(sigma)*C;
    X = gam*vartheta+sqrt(gam')*Z*S;
    [K,~] = size(x);
    T = zeros(K,1);
    if sum(sum(isnan(C)))>0
        T = 1e5*ones(K,1);
    else     
        for k=1:K
%           for m=1:M
%           Pr_p = (cp+1/zeta).*expint( (x(k,:)+X(m,:))./bp );
%           Pr_n = (cn+1/zeta).*expint( (x(k,:)+X(m,:))./bn );
%           Pr_p(isnan(Pr_p)) = 0; Pr_p(~isfinite(Pr_p)) = 0;
%           Pr_n(isnan(Pr_n)) = 0; Pr_n(~isfinite(Pr_n)) = 0;
%           T(k) = T(k)+prod(Pr_p+Pr_n)/M;
%           end
        
            Pr_p = zeros(M,D);
            Pr_n = zeros(M,D);
%           for m=1:M
%           Pr_p(m,:) = (cp+1/zeta).*expint( (x(k,:)+X(m,:))./bp );
%           end
            for d=1:D
                Pr_p(:,d) = (cp(d)+1/zeta).*expint( (x(k,d)+X(:,d))./bp(d) );
                Pr_n(:,d) = (cn(d)+1/zeta).*expint( (x(k,d)+X(:,d))./bn(d) );
            end
%           Pr_p = (cp+1/zeta) .* expint( (x(k,:)+X)./bp );
%             Pr_n = (cn+1/zeta) .* expint( (x(k,:)+X)./bn );
%           fprintf('pos = (%d,%d), neg = (%d,%d)\n',size(Pr_p),size(Pr_n))
        
            Pr_p=Pr_p(imag(Pr_p)==0); Pr_n=Pr_n(imag(Pr_n)==0);
            T(k) = sum(prod(Pr_p,2)+prod(Pr_n,2))/M;
        end
    end
end

function T = obj(param,zeta_min,Z,eps,bp,cp,bn,cn,x,TT)
    fprintf('g0 = [%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,...\n      %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f];\nzeta0 = %.4f;\n', param)
    g = param(1:end-1);
    zeta = param(end);
    zeta = max(zeta,zeta_min);
    hatTT = hatT(g,zeta,Z,eps,bp,cp,bn,cn,x);
    T = sum( ((TT-hatTT).^2)./(TT.*(1-TT)) );
    %[TT, hatTT]
    fprintf('T = %.4f\n\n', T)
end

%%

% g0 = [-256.0000,319.0000,0.0000,192.0000,0.0000,512.0000,0.0000,0.0000,0.0000,-512.0000,0.0000,...
%       0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,...
%       0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,...
%       0.0000,0.0000,0.0000,512.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,...
%       0.0000,0.0000,0.0000,0.0000,0.0000,512.0000,0.0000,0.0000,0.0000,0.0000,0.0000];
% zeta0 = 1.5831;
% 
% g0 = [-138.9872,248.8192,0.0010,175.3614,0.0016,496.2071,0.0019,-0.0011,0.0061,-330.8070,-0.0005,...
%       0.0031,0.0052,-0.0030,0.0009,-0.0006,-0.0021,0.0020,-0.0005,-0.0037,-0.0007,0.0008,...
%       -0.0009,0.0059,-0.0106,0.0023,-0.0032,-0.0020,0.0013,0.0017,0.0016,-0.0014,-0.0092,...
%       0.0102,-0.0073,-0.0031,273.4575,-0.0027,0.0029,0.0008,-0.0018,-0.0021,-0.0012,0.0005,...
%       0.0010,-0.0013,-0.0024,0.0017,-0.0014,3.4210,0.0011,0.0054,0.0146,0.0024,0.0011];
% zeta0 = 1.7610;
% 
% g0 = [-138.9872,248.8192,0.0010,175.3614,252.0016,496.2071,0.0331,0.0311,0.0686,-330.8070,-127.0288,...
%       0.0031,0.0052,-0.0030,128.0009,-0.0006,-0.0021,0.0020,-0.0005,-0.0037,-0.0007,0.0008,...
%       255.9991,0.0059,-0.0106,0.0023,-0.0032,-0.0020,0.0013,0.0017,0.0016,-0.0014,-0.0092,...
%       0.0102,-0.0073,-0.0031,273.4575,-0.0027,0.0029,0.0008,-0.0018,-0.0021,-0.0012,0.0005,...
%       0.0010,-512.0013,-0.0024,0.0017,-0.0014,3.4210,0.0011,0.0054,0.0146,0.0024,0.0011];
% zeta0 = 1.5831;
% 
% g0 = [-208.7952,164.5707,0.0012,86.1283,3.0301,54.9208,-0.0231,0.0170,0.0943,-204.5011,0.1618,...
%       -0.0015,0.0033,-0.0019,536.4822,-0.0008,-0.0084,0.0013,-0.0005,-0.0024,-0.0007,0.0006,...
%       390.8303,0.0012,0.0117,0.0023,-0.0013,-0.0038,-0.0008,0.0030,0.0033,-0.0091,-0.0071,...
%       -0.0021,-0.0233,0.0028,303.3318,-0.0016,0.0039,0.0010,-0.0027,-0.0017,-0.0012,0.0003,...
%       0.0007,-569.1018,-0.0038,-0.0002,-0.0010,3.3468,0.0008,0.0055,0.0467,-0.0013,0.0010];
% zeta0 = 1.7179;
% 
% g0 = [47.2048,164.5707,0.0012,86.1283,3.0301,54.9208,-0.0231,0.1108,0.0943,-204.5011,15.1618,...
%       -0.0015,0.0033,1.7481,536.4822,-0.0008,-0.0084,0.0013,-0.0005,-0.0024,-0.0007,0.0006,...
%       390.8303,0.0012,0.0117,0.0023,-0.0013,-0.0038,-0.0008,0.0030,0.0033,-0.0091,-0.0071,...
%       -0.0021,-0.0233,0.0028,303.3318,-0.0016,0.0039,0.0010,-0.0027,-4.0017,-0.1262,0.0003,...
%       0.0007,-569.1018,-0.0038,-0.0002,-0.0010,3.3468,0.0008,0.0055,0.0467,-0.0013,0.0010];
% zeta0 = 1.7179;
% 
% g0 = [66.5368,159.3966,0.0012,82.5168,2.8912,54.6468,-0.0204,0.0984,0.0811,-219.5003,14.5076,...
%       -0.0030,0.0036,0.3805,506.3880,-0.0010,-0.0085,0.0012,-0.0001,-0.0028,-0.0012,0.0006,...
%       396.3557,0.0018,0.0216,-0.0011,-0.0010,-0.0030,-0.0009,0.0021,0.0019,-0.0092,-0.0061,...
%       -0.0028,-0.0253,0.0036,307.0303,-0.0017,0.0030,0.0011,-0.0026,-5.3580,-0.1483,0.0004,...
%       0.0008,-600.3056,0.0003,-0.0001,-0.0015,3.3770,0.0006,0.0066,0.0440,-0.0009,0.0018];
% zeta0 = 1.7236;
% 
% g0 = [66.5368,159.3966,0.0012,-173.4832,2.8912,54.6468,-0.0204,0.0984,0.0811,-219.5003,142.5076,...
%       -0.0030,0.0036,505.8805,506.3880,-0.0010,-0.0085,0.0012,0.2343,-0.0028,-0.0012,0.0006,...
%       396.3557,0.0018,0.0216,-0.0011,-0.0010,-0.0030,-0.0009,0.0021,0.0019,-0.0092,-0.0061,...
%       -0.0028,-0.0253,0.0036,307.0303,-0.0017,0.0030,0.0011,-0.0026,-5.3580,-0.1483,0.0004,...
%       0.0008,-600.3056,-0.2497,-0.0157,-0.0015,3.3770,0.0006,0.0066,0.0440,-0.0009,-0.9982];
% zeta0 = 1.7236;
% 
% g0 = [63.3282,183.6605,0.0022,-206.3141,4.7604,24.7824,-0.0233,0.0938,0.0787,-207.0250,159.0497,...
%       -0.0019,0.0016,546.5750,534.9688,-0.0011,-0.0052,0.0009,0.2641,-0.0046,-0.0006,0.0005,...
%       316.3509,0.0016,0.0223,-0.0030,-0.0014,-0.0002,-0.0010,-0.0009,0.0020,-0.0091,-0.0055,...
%       -0.0045,-0.0142,0.0071,317.2489,-0.0037,0.0024,0.0002,-0.0017,-5.8065,-0.1115,0.0004,...
%       0.0011,-527.0614,-0.3437,-0.0151,-0.0022,1.9091,0.0006,0.0081,0.0417,-0.0005,-0.1628];
% zeta0 = 1.7572;
% 
% g0 = [63.3282,439.6605,0.0022,-206.3141,4.7604,24.7824,-0.0233,0.0938,0.0787,-207.0250,159.0497,...
%       -0.0019,0.0016,546.5750,534.9688,-0.0011,-0.0052,0.0302,0.2641,-0.0046,-0.0006,-0.0620,...
%       316.3509,0.0641,0.0223,-0.0030,-0.0014,-0.2190,-0.0010,-0.0009,0.0020,-0.0091,-0.0055,...
%       -0.0045,-0.0142,0.0071,317.2489,-0.0037,0.0024,0.0236,-0.0017,-5.8065,-0.1115,0.0004,...
%       0.0011,-527.0614,-0.3437,-0.0151,-0.0022,65.9091,0.0006,0.0081,0.0417,-0.0005,-0.1628];
% zeta0 = 1.7572;
