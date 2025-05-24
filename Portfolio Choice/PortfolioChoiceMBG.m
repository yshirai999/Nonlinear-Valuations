%% Portfolio Choice
clear
clc
close all

thetai0 = 'prev'; %% 'unif' for uniform initial condition, 'prev' for previous

DeltaBG = load('DeltaBG');
DeltaBG = DeltaBG.DeltaS;

%% Testing Parameters.
% MBG parameters
ticker = {'SPY', 'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY'}; % materials, energy, industrial, financials, tech, consumer staples, utilities, healthcare, consumer discretionary
ETF = matlab.lang.makeValidName(ticker);
D = length(ticker);

dataPath = getPath('Data');

for d=1:D
    Y = load(fullfile(dataPath, strcat('Y',ticker{d},'.mat')));
    Y = Y.Y;
    ind = (Y(:,1)>20151231);
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

% Distortion Parameters
c = 0.01;
gam = 0.25;
aw = 1000;
bw = 0.1;

%% Maximization

varPath = getPath('VarArchive');

theta = zeros(N,D);
r_p = zeros(N,1);
r_SPY = zeros(N,1);

mu_p = zeros(N,1);
mu_SPY = zeros(N,1);
sigma_p = zeros(N,1);
sigma_SPY = zeros(N,1);
Sharpe_p = zeros(N,1);
Sharpe_SPY = zeros(N,1);

thetai = ones(D-1,1)/D;

try
    theta = load(strcat('thetaMBG',thetai0,'.mat'));
    theta = theta.theta;
catch
    for i = 1:N-1
        % MBG Parameters
    
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
        
        %draw sample of size M from ellipsoid {z'Qz=1}
        Unif = rand(M,D-1)*2*pi; %angles of the ellipsoid {z'Qz=1}
        R = gamrnd(gamma_sim,1/sqrt(delta),M,1); %ray of the ellipsoid {z'Qz=1}
        U = ones(M,D).*L_Qinv';
%         U(:,1:D-1) = cos(Unif(:,1)).*ones(M,D-1);
%         U(:,D) = sin(Unif(:,1));
%         
%         for dd = 2:D-1
%             try
%                 U(:,dd:D) = U(:,dd:D).*sin(Unif(:,dd));
%             catch
%             end
%         end
        
        for dd = 1:D-1
            U(:,dd) = U(:,dd).*cos(Unif(:,dd));
            U(:,dd+1:D) = U(:,dd+1:D).*sin(Unif(:,dd));
        end
        


%         figure
%         try
%             hold on
%             scatter3(U(:,1),U(:,2),U(:,3))
%             hold off
%         catch
%             figure
%             hold on
%             scatter(U(:,1),U(:,2))
%             hold off
%         end
%         
        U = (V_Q*U')';
        Z = R.*U;
%         try
%             figure
%             hold on
%             scatter3(Z(:,1),Z(:,2),Z(:,3))
%             hold off
%         catch          
%             figure
%             hold on
%             scatter(Z(:,1),Z(:,2))
%             hold off
%         end
    
        % Maximization
        
        rho = rhoVG(delta,vartheta,zeta0,Z,gamma_sim);
        
        options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
        f = @(theta)-H([theta;1-sum(theta)],c,gam,aw,bw,bp,cp,bn,cn,a,Z,rho);
    
    %     thetaplot = 0.01:0.01:0.99;
    %     ff = zeros(size(thetaplot));
    %     for ii = 1:length(thetaplot)
    %         ff(ii) = f(thetaplot(ii));
    %     end
    %     figure
    %     plot(thetaplot,-ff);
    %     
    %     [~,ii] = min(ff);
    %     thetai = thetaplot(ii);
    
        if strcmp(thetai0,'unif')
            thetai = ones(D-1,1)/D;
        elseif strcmp(thetai0,'prev')
        else
            fprintf('Warning: initial conditions not correctly set')
        end
    
        [thetai,Sopt,exitflagBG] = fminsearch(f,thetai,options); 
        
        %theta = theta-sign(theta).*floor(abs(theta)); 
        thetai = thetai-floor(thetai); %Short ONLY
        
        theta(i,:) = [thetai',1-sum(thetai)];
    
        r_p(i) = theta(i,:)*(Y(i+1,:)'-Y(i,:)');
        r_SPY(i) = (Y(i+1,1)'-Y(i,1)');
        
        imin = max(i-252,1);
        mu_p(i) = mean(r_p(1:i));
        sigma_p(i) = std(r_p(1:i));
    
        mu_SPY(i) = mean(r_SPY(1:i));
        sigma_SPY(i) = std(r_SPY(1:i));
    
        Sharpe_p(i) = mu_p(i)./sigma_p(i);
        Sharpe_SPY(i) = mu_SPY(i)./sigma_SPY(i);
    
        fprintf('i = %d, theta_i(1,2) = (%d,%d), Sharpe = (%d,%d)\n',...
                    i, theta(i,1), theta(i,2), Sharpe_p(i), Sharpe_SPY(i))
        
    end    
end
%% Visualization

vizPath = getPath('VizArchive');

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
    r_p(i) = theta(i,:)*(Y(i+1,:)'-Y(i,:)');
    r_SPY(i) = (Y(i+1,1)'-Y(i,1)');
    
    r_pmon(i) = sum(r_p(max(i-Delta+1,1):i));
    r_SPYmon(i) = sum(r_SPY(max(i-Delta+1,1):i));
    
%     r_pmon(i) = r_p(i);
%     r_SPYmon(i) = r_SPY(i);
end


per = 252; %Lookback period to estimate mean and stdev of daily returns
for i=1:N-1
     imin = max(i-252,1); 
     
     mu_p(i) = mean(r_pmon(imin:i));
     sigma_p(i) = std(r_pmon(imin:i));
     
     mu_SPY(i) = mean(r_SPYmon(imin:i));
     sigma_SPY(i) = std(r_SPYmon(imin:i));
     
     Sharpe_p(i) = sqrt(252)*mu_p(i)./sigma_p(i);
     Sharpe_SPY(i) = sqrt(252)*mu_SPY(i)./sigma_SPY(i);     
     
end

fprintf('Portfolio Sharpe > SPY Sharpe %d, and Sharpe correlation is %d\n',...
    mean(Sharpe_p(per:N-1)>Sharpe_SPY(per:N-1)), corr(Sharpe_p(per:N-1),Sharpe_SPY(per:N-1)))

[quantile(Sharpe_p,[0.05;0.25;0.5;0.75;0.95]),quantile(Sharpe_SPY,[0.05;0.25;0.5;0.75;0.95])]

% Sharpe_p_out = rmoutliers(Sharpe_p,'percentiles',[2.5,97.5]);
% Sharpe_SPY_out = rmoutliers(Sharpe_SPY,'percentiles',[2.5,97.5]);
%[quantile(Sharpe_p_out(2:end),[0;0.25;0.5;0.75;1]),quantile(Sharpe_SPY_out(2:end),[0;0.25;0.5;0.75;1])]

figure 
hold on
grid on
box on
plot(dates(per:N-1),Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1))
%plot(dates(2*252:N-1),(mup(2*252:N-1,1)-mun(2*252:N-1,1))./sqrt(sigmap(2*252:N-1,1).^2+sigman(2*252:N-1,1).^2))
str=strcat('AnnualSharpeMBGDelta');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

figure 
hold on
grid on
box on
plot(dates(per:N-1),Sharpe_p(per:N-1,1))
plot(dates(per:N-1),Sharpe_SPY(per:N-1,1))
%plot(dates(2*252:N-1),(mup(2*252:N-1,1)-mun(2*252:N-1,1))./sqrt(sigmap(2*252:N-1,1).^2+sigman(2*252:N-1,1).^2))
legend('Portfolio','SPY')
str=strcat('AnnualSharpeMBG');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

figure 
hold on
grid on
box on
plot(dates(per:N-1),Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1))
plot(dates(per:N-1),DeltaBG)
legend('MBG','BG')
str=strcat('AnnualSharpeMBGBG');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

figure 
hold on
grid on
box on
plot(dates(per:N-1),r_pmon(per:N-1),'LineWidth', 2)
plot(dates(per:N-1),r_SPYmon(per:N-1),'-.','LineWidth', 0.5)
legend('Portfolio','SPY')
title('Daily Returns')
str=strcat('AnnualReturnMBG');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

%% Save

save(fullfile(varPath, strcat('thetaMBG',thetai0)), 'theta')

DeltaS = Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1);
save(fullfile(varPath, strcat('DeltaBG',thetai0)), 'DeltaS')

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
    d = bw*exp( -c * x );
end
 
function d = dpsip(x,c,gamma,aw)
    d = (aw*c/(1+gamma))...
        .* ( ( 1-exp( - c * x ) ) .^ (-gamma/(1+gamma)) )...
        .* exp( -c * x );
    d(~isfinite(d)) = 0;
end

function k = kappa(x,theta,bp,cp,bn,cn)
    k = theta.*(exp(x)-1).*...
        ( cp * ( exp(-abs(x)/bp) ./ abs(x) ) .* (x>0)...
        + cn * ( exp(-abs(x)/bn) ./ abs(x) ) .* (x<0));
    k(isnan(k))=0;
end

function t = H(theta,c,gamma,aw,bw,bp,cp,bn,cn,a,Z,rho)
    [M,D] = size(Z);
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
    
    linvar = theta'*a + sum ( (theta'*(exp(Z')-1))' .* rho )/M;
    t = linvar-sum(t);
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