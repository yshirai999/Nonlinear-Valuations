%% Convex Portfolio Choice and Investment Amount Optimization
clear
clc
close all

thetai0 = 'unif'; %% 'unif' for uniform initial condition, 'prev' for previous
LongShort = '';

% DeltaBG = load('DeltaBG');
% DeltaBG = DeltaBG.DeltaS;

%% Testing Parameters.
% MBG parameters
ticker = {'SPY', 'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY'}; % materials, energy, industrial, financials, tech, consumer staples, utilities, healthcare, consumer discretionary
ETF = matlab.lang.makeValidName(ticker);
D = length(ticker);

for d=1:D

Y = load(strcat('Y',ticker{d},'.mat'));
Y = Y.Y;
ind = logical((Y(:,1)>20200131).*(Y(:,1)<20200501));
eval([ETF{d} ' = Y(ind,2:end);']);

end

datesnum = Y(ind,1);
dates = datetime(Y(ind,1),'ConvertFrom','yyyymmdd');

N = length(dates);

bpvec = zeros(N,D);
cpvec = zeros(N,D);
bnvec = zeros(N,D);
cnvec = zeros(N,D);
Y = zeros(N,D);

for d = 1:D
    eval(['bpvec(:,d) = ' ETF{d} '(:,4);']);
    eval(['cpvec(:,d) = ' ETF{d} '(:,5);']);
    eval(['bnvec(:,d) = ' ETF{d} '(:,6);']);
    eval(['cnvec(:,d) = ' ETF{d} '(:,7);']);
    eval(['Y(:,d) = ' ETF{d} '(:,2);']);
end

vis = 1;

mup = zeros(N,D);
sigmap = zeros(N,D);
mun = zeros(N,D);
sigman = zeros(N,D);

for i = 1 : N
    mup(i,:) = ( bpvec(i,:).*cpvec(i,:) )';
    sigmap(i,:) = ( bpvec(i,:).*sqrt(cpvec(i,:)) )';
    mun(i,:) = ( bnvec(i,:).*cnvec(i,:) )';
    sigman(i,:) = ( bnvec(i,:).*sqrt(cnvec(i,:)) )';
end


g0 = [60.3465,651.3773,-0.0038,-205.6770,4.8418,24.2426,-0.0664,0.0613,0.0595,-290.6120,161.5743,...
      -0.0038,0.0027,421.5420,535.0741,-0.0009,0.0145,0.0159,0.0747,0.0021,-0.0016,-0.0722,...
      315.5769,-0.0080,0.0150,-0.0045,-0.0013,-0.2185,-0.0011,-0.0009,-0.0043,-0.0007,-0.0084,...
      0.0105,-0.0185,-0.0012,282.6299,0.0075,0.0007,0.0169,-0.0043,-6.9368,0.1568,0.0002,...
      0.0024,-594.3617,-5.4250,-0.0153,-0.0027,253.6199,0.0012,0.0096,0.0323,0.0007,-0.1906];
zeta0 = 1.8209;

eps = 1e-10;
C = C_AH(g0,D+1,eps);
C(:,2)=[]; %remove VIX
C(2,:)=[]; %remove VIX
gamma_sim = 1;
M = 10000; %sample size for montecarlo integration

% Rebate and Distortion Parameters
cu = 100;
cl = 2;
gam = 0.01;
aw = 1000;
bw = 0.1;

chi = [250];
chi2 = 2*[1];

% Term and Investment Amounts
T = 1;

%% Maximization

theta = zeros(N,D);
P = zeros(N,1);
r_p = zeros(N,1);
r_SPY = zeros(N,1);

mu_p = zeros(N,1);
mu_SPY = zeros(N,1);
sigma_p = zeros(N,1);
sigma_SPY = zeros(N,1);
Sharpe_p = zeros(N,1);
Sharpe_SPY = zeros(N,1);

thetai = [ones(D-1,1)/D];
Pi = 2500000;

try
    r = load(strcat('ConvexMBG_FebApr2020',thetai0,'_r',LongShort,'.mat'));
    theta = load(strcat('ConvexMBG_FebApr2020',thetai0,'_theta',LongShort,'.mat'));
    P = load(strcat('ConvexMBG_FebApr2020',thetai0,'_P',LongShort,'.mat'));
    r = r.r;
    theta = theta.theta;
    P = P.P;
    r_p = r(:,1);
    r_SPY = r(:,2);
catch
    alpha = [chi,chi2];
    for i = 1:N-1
        % MBG Parameters
        tic;
    
        bp = bpvec(i,:)';
        cp = cpvec(i,:)';
        bn = bnvec(i,:)';
        cn = cnvec(i,:)';
        
        sigma = sqrt(2*bp.*bn/zeta0);
        vartheta = (bp-bn)/zeta0;
        S = diag(sigma)*C*C*diag(sigma);
        d = sqrt(zeta0)*eye(D);
        Q = pinv(d*S*d); 
        [V_Q,L_Q] = eig(Q);
        [V_Qinv, L_Qinv] = eig(d*S*d);
        L_Qinv = max(real(diag(L_Qinv)),0);
        V_Q = real(V_Q);
        
        delta = vartheta'*Q*vartheta+2;
    
        a = log( ((1./(1-bp)).^cp).*((1./(1+bn)).^cn) );
        
        rng('default'); %freeze the seed
        rng(1);
        Unif = rand(M,D-1)*2*pi; %angles of the ellipsoid {z'Qz=1}
        R = gamrnd(gamma_sim,1/sqrt(delta),M,1); %ray of the ellipsoid {z'Qz=1}
        U = ones(M,D).*L_Qinv';
        
        for dd = 1:D-1
            U(:,dd) = U(:,dd).*cos(Unif(:,dd));
            U(:,dd+1:D) = U(:,dd+1:D).*sin(Unif(:,dd));
        end
        
        U = (V_Q*U')';
        Z = R.*U;

        % Maximization
        
        rho = rhoVG(delta,vartheta,zeta0,Z,gamma_sim);
        
        options = optimset('MaxFunEval',100,'MaxIter',100,'Display','iter','TolFun',1e-9,'TolX',1e-9);
        
        if strcmp(LongShort,'_Short')
            f = @(x)-H([x(1:end-1);1-sum(x(1:end-1))],x(end),gam,aw,bw,bp,cp,bn,cn,a,Z,rho,cu,cl,alpha,T);
        else
            f = @(x)-H([x(1:end-1);1-sum(x(1:end-1))],abs(x(end)),gam,aw,bw,bp,cp,bn,cn,a,Z,rho,cu,cl,alpha,T);
        end
        
        if strcmp(thetai0,'unif')
            thetai = [ones(D-1,1)/D];
            Pi = 2500000;
        elseif strcmp(thetai0,'prev')
        else
            fprintf('Warning: initial conditions not correctly set')
        end

        [thetap,Sopt,exitflagBG] = fminsearch(f,[thetai;Pi],options); 
        
        thetai(1:end) = thetap(1:end-1)-floor(thetap(1:end-1)); %Short ONLY
        if strcmp(LongShort,'_Short')
            Pi = thetap(end);
        else
            Pi = abs(thetap(end));
        end
        
        theta(i,:) = [thetai',1-sum(thetai)];
        P(i) = Pi;
    
        r_p(i) = P(i)*theta(i,:)*(Y(i+1,:)'-Y(i,:)');
        r_SPY(i) = (Y(i+1,1)'-Y(i,1)');
        
        imin = max(i-252,1);
        mu_p(i) = mean(r_p(1:i));
        sigma_p(i) = std(r_p(1:i));
    
        mu_SPY(i) = mean(r_SPY(1:i));
        sigma_SPY(i) = std(r_SPY(1:i));
    
        Sharpe_p(i) = mu_p(i)./sigma_p(i);
        Sharpe_SPY(i) = mu_SPY(i)./sigma_SPY(i);

        cputime = toc;
    
        fprintf('i = %d, theta_i(1,2) = (%d,%d), Sharpe = (%d,%d), P(i) = %d, time = %d\n',...
                    i, theta(i,1), theta(i,2), Sharpe_p(i), Sharpe_SPY(i), P(i), cputime)
        
    end    
end


%% Visualization

r_p = zeros(N,1);
r_SPY = zeros(N,1);

r_pmon = zeros(N,1);
r_SPYmon = zeros(N,1);

mu_p = zeros(N,1);
mu_SPY = zeros(N,1);
sigma_p = zeros(N,1);
sigma_SPY = zeros(N,1);

Sharpe_p = zeros(N,1);
Sharpe_SPY = zeros(N,1);

Delta = 1;

for i=1:N-1
    r_p(i) = P(i) * theta(i,:) * ( (Y(i+1,:)'-Y(i,:)') ./ Y(i,:)' );
    r_SPY(i) = 10 * 1e6 * (Y(i+1,1)'-Y(i,1)') / Y(i,1)';
    
    r_pmon(i) = sum(r_p(max(i-Delta+1,1):i));
    r_SPYmon(i) = sum(r_SPY(max(i-Delta+1,1):i));
    
%     r_pmon(i) = r_p(i);
%     r_SPYmon(i) = r_SPY(i);
end
% 
per = 1; %Lookback period to estimate mean and stdev of daily returns
for i=1:N-1
     imin = max(i-per,1); 
     
     mu_p(i) = mean(r_p(imin:i));
     sigma_p(i) = std(r_p(imin:i));
     
     mu_SPY(i) = mean(r_SPY(imin:i));
     sigma_SPY(i) = std(r_SPY(imin:i));
     
     Sharpe_p(i) = sqrt(per)*mu_p(i)./sigma_p(i);
     Sharpe_SPY(i) = sqrt(per)*mu_SPY(i)./sigma_SPY(i);     
     
end

% fprintf('Portfolio Sharpe > SPY Sharpe %d, and Sharpe correlation is %d\n',...
%     mean(Sharpe_p(per:N-1)>Sharpe_SPY(per:N-1)), corr(Sharpe_p(per:N-1),Sharpe_SPY(per:N-1)))

% for i = 1:N
%     mu_p(i) = mean(Pr_p(1:i));
%     sigma_p(i) = std(Pr_p(1:i));
%     
%     mu_SPY(i) = mean(5*1e6*r_SPY(1:i));
%     sigma_SPY(i) = std(5*1e6*r_SPY(1:i));
%     
%     Sharpe_p(i) = mu_p(i)./sigma_p(i);
%     Sharpe_SPY(i) = mu_SPY(i)./sigma_SPY(i);
% end

figure
hold on
box on
grid on
plot(dates,P)
datetick('x','mmm yyyy')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('ConvexMBGFebApr2020_P',LongShort);
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
fpath=('C:\Users\Yoshihiro Shirai\Desktop\Spectral Martingale Measures');
saveas(gcf, fullfile(fpath, fname), 'epsc');
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Dissertation\Dissertation');
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

figure
hold on
box on
grid on
plot(dates,cumsum(r_p),'-')
plot(dates,cumsum(r_SPY),'-.')
datetick('x','mmm yyyy')
legend('Portfolio','SPY','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('ConvexMBGFebApr2020_r',LongShort);
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
fpath=('C:\Users\Yoshihiro Shirai\Desktop\Spectral Martingale Measures');
saveas(gcf, fullfile(fpath, fname), 'epsc');
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Dissertation\Dissertation');
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

figure
hold on
box on
grid on
plot(dates,Y(:,1))
legend('SPY','interpreter','latex')
datetick('x','mmm yyyy')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('SPY');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
fpath=('C:\Users\Yoshihiro Shirai\Desktop\Spectral Martingale Measures');
saveas(gcf, fullfile(fpath, fname), 'epsc');
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Dissertation\Dissertation');
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

Nmin = 1;
Nmax = N;

[quantile(Sharpe_p(Nmin:Nmax),[0.05;0.25;0.5;0.75;0.95]),quantile(Sharpe_SPY(Nmin:Nmax),[0.05;0.25;0.5;0.75;0.95])]

Sharpe_p_out = rmoutliers(Sharpe_p,'percentiles',[2.5,97.5]);
Sharpe_SPY_out = rmoutliers(Sharpe_SPY,'percentiles',[2.5,97.5]);
[quantile(Sharpe_p_out(2:end),[0;0.25;0.5;0.75;1]),quantile(Sharpe_SPY_out(2:end),[0;0.25;0.5;0.75;1])]

%% Save

prompt = 'Do you want to save results? Y/N: ';
s = input(prompt, 's');
%s = 'Y';
if strcmp(s,'Y')
    prompt = 'Warning: data will be saved and previous one overwritten. Input Y to continue: ';
    s = input(prompt, 's');
    if strcmp(s,'Y')
        r = [r_p,r_SPY];
        save(strcat('ConvexMBG_FebApr2020',thetai0,'_theta',LongShort),'theta')
        save(strcat('ConvexMBG_FebApr2020',thetai0,'_P',LongShort),'P')
        save(strcat('ConvexMBG_FebApr2020',thetai0,'_r',LongShort),'r')
    end
end

% 
% DeltaS = Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1);
% save(strcat('DeltaMBGConvex',thetai0),'DeltaS')

%% Routines

function n = nup(p,theta,bp,cp,bn,cn)
    n = zeros(size(p));
    ind = (p>0);
    for k=1:length(theta)
        if theta(k)>0
            n(ind) = n(ind)+cp(k)*expint(log(p(ind)/theta(k)+1)/bp(k));
            n(isnan(n)) = 0;
        elseif theta(k)<0
            indn = (theta(k)<-p);
            ind = logical(indn.*ind);
            n(ind) = n(ind)+cn(k)*expint(-log(p(ind)/theta(k)+1)/bn(k));
        end
    end
end

function n = nun(p,theta,bp,cp,bn,cn)
    n = zeros(size(p));
    ind = (p<0);
    for k=1:length(theta)
        if theta(k)<0
            n(ind) = n(ind)+cp(k)*expint(log(p(ind)/theta(k)+1)/bp(k));
            n(isnan(n)) = 0;
        elseif theta(k)>0
            indn = (theta(k)>-p);
            indk = logical(indn.*ind);
            n(indk) = n(indk) + cn(k)*expint(-log(p(indk)/theta(k)+1)/bn(k));
        end
    end
end

function d = dpsin(x,c,bw)
    %d = bw*exp( -c * x );
    d = exp( -c * x );
end
 
function d = dpsip(x,c,gamma,aw)
%     d = (aw*c/(1+gamma))...
%         .* ( ( 1-exp( - c * x ) ) .^ (-gamma/(1+gamma)) )...
%         .* exp( -c * x );
    d = (1/(1+gamma))...
        .* ( ( 1-exp( - c * x ) ) .^ (-gamma/(1+gamma)) )...
        .* exp( -c * x ).*(x>0);
    d(isnan(d))=0;
    d(~isfinite(d)) = 0;
end

function k = kappa(x,theta,bp,cp,bn,cn)
    k = theta.*(exp(x)-1).*...
        ( cp * ( exp(-abs(x)/bp) ./ abs(x) ) .* (x>0)...
        + cn * ( exp(-abs(x)/bn) ./ abs(x) ) .* (x<0));
    k(isnan(k))=0;
    k(~isfinite(k)) = 0;
end

function t = Hc(theta,p,c,gamma,aw,bw,bp,cp,bn,cn,cu,cl,alpha,T,linvar,D)
    theta = theta-floor(theta);
    %theta = theta-sign(theta).*floor(abs(theta));
    t = zeros(D,1);
    for j = 1:D
        funp = @(yj) dpsin( nup(theta(j)*(exp(yj)-1),theta,bp,cp,bn,cn) , c, bw ).*...
                    kappa(yj,theta(j),bp(j),cp(j),bn(j),cn(j));
                
        funn = @(yj) dpsip( nun(theta(j)*(exp(yj)-1),theta,bp,cp,bn,cn) , c , gamma, aw).*...
                    kappa(yj,theta(j),bp(j),cp(j),bn(j),cn(j));
        
        Ip = integral(funp,0,Inf);
        In = integral(funn,-Inf,0);
        t(j) = Ip - In;
    end
    
    t = p*linvar-p*sum(t)+beta(c,cu,cl,alpha)/(T);
end

function t = H(theta,p,gamma,aw,bw,bp,cp,bn,cn,a,Z,rho,cu,cl,alpha,T)
    [M,D] = size(Z);
    linvar = theta'*a + sum ( (theta'*(exp(Z')-1))' .* rho )/M;
    f = @(c)Hc(theta,p,c,gamma,aw,bw,bp,cp,bn,cn,cu,cl,alpha,T,linvar,D);
    options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
    Nc = 10;
    c = linspace(cl,cu,Nc);
    fi = zeros(Nc,1);
    for i = 1:Nc
        fi(i) = f(c(i));
    end
    [~,cind] = min(fi);
    c = c(cind);
    c = fminsearch(f,c,options);
    t = f(c);
end

function rho = rhoVG(delta,vartheta,zeta,Z,gamma_sim)
    [~,D] = size(Z);
    rho = 2*pi * sqrt(delta) * (zeta^(D/2-2/3)) * exp(vartheta'*Z')'...
          * gamma(gamma_sim) ./ ( gamma(D/2) * (delta^(gamma_sim/2)) .* (diag(Z*Z')).^(gamma_sim/2) );
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

function b = beta(c,cu,cl,alpha)
    %b = alpha(1)*exp( 1./(c-cl).^alpha(2) ).*exp(-1./(cu-c).^alpha(2)).*(c<cu);
    
    b = ( exp(alpha(1)*(c-cu*1.5))-1-alpha(1)*(c-cu*1.5) ) ./ (c-cl*0.7).^alpha(2);
    
    b(~isfinite(b)) = 0;
    b(isnan(b)) = 0;
end