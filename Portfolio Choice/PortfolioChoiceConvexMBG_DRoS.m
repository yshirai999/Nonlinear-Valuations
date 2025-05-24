%% Convex Portfolio Choice and Investment Amount Optimization
clear
clc
close all

thetai0 = 'prev'; %% 'unif' for uniform initial condition, 'prev' for previous

% DeltaBG = load('DeltaBG');
% DeltaBG = DeltaBG.DeltaS;

%% Testing Parameters.
% MBG parameters
ticker = {'SPY', 'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY'}; % materials, energy, financials, industrial, tech, consumer staples, utilities, healthcare, consumer discretionary
%ticker = {'SPY', 'XLI'};
ETF = matlab.lang.makeValidName(ticker);
D = length(ticker);

dataPath = getPath('Data');

for d=1:D
    Y = load(fullfile(dataPath, strcat('Y',ticker{d},'.mat')));
    Y = Y.Y;
    ind = logical((Y(:,1)>20200227).*(Y(:,1)<20200303));
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
C = C_AH(g0,11,eps);
%C(:,[2:5,7:11])=[]; %remove VIX & others
%C([2:5,7:11],:)=[]; %remove VIX & others
C(:,2)=[]; %remove VIX & others
C(2,:)=[]; %remove VIX & others
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

varPath = getPath('VarArchive');

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

thetai = ones(D-1,1)/D;
Np = 100;
Pmin = 500;
%Pmax = 100000000;
Pmax = 10000000;
Pi = linspace(Pmin,Pmax,Np)';

try
    theta = load(fullfile(varPath, strcat('thetaMBG_Convex',thetai0,'.mat')));
    theta = theta.theta;
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
        thetaik = zeros(Np,D-1);
        fi = zeros(Np,1);
        for k=1:Np
            f = @(theta)-H([theta;1-sum(theta)],Pi(k),gam,aw,bw,bp,cp,bn,cn,a,Z,rho,cu,cl,alpha,T);            
            %[thetaik(k,:),fi(k)] = fminsearch(f,thetai,options);
            fi(k) = -f(thetai);
            fprintf('P = %d, theta(1,-1)=(%d,%d),fik = %d\n',Pi(k),thetai(1),1-thetai(1),fi(k))
        end
        
%         for k=1:Np
%             f = @(theta)-H([theta;1-sum(theta)],Pi(k),gam,aw,bw,bp,cp,bn,cn,a,Z,rho,cu,cl,alpha,T);            
%             %[thetaik(k,:),fi(k)] = fminsearch(f,thetai,options);
%             %fi(k) = -f([0.2,0.05,0.05,0.1,0.1,0.1,0.1,0.1,0.1]');
%             fi(k) = -f([0.5,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04]'); %becomes negative with cu=100 and cl = 20;
%             fprintf('P = %d, theta(1,-1)=(%d,%d),fik = %d\n',Pi(k),thetai(1),1-thetai(1),fi(k))
%         end
%         plot(Pi,fi);
        grid on
        box on
        hold off
     end    
end

%% Visualization
vizPath = getPath('Visualization');
figure 
hold on
grid on
box on
plot(Pi,fi);
xlabel('$\varpi$', 'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('MBG_thetaunif_02282020');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

%% Visualization

% r_p = zeros(N,1);
% r_SPY = zeros(N,1);
% 
% r_pmon = zeros(N,1);
% r_SPYmon = zeros(N,1);
% 
% mu_p = zeros(N,1);
% mu_SPY = zeros(N,1);
% sigma_p = zeros(N,1);
% sigma_SPY = zeros(N,1);
% 
% Sharpe_p = zeros(N,1);
% Sharpe_SPY = zeros(N,1);
% 
% Delta = 1;
% 
% for i=1:N-1
%     r_p(i) = theta(i,:)*(Y(i+1,:)'-Y(i,:)');
%     r_SPY(i) = (Y(i+1,1)'-Y(i,1)');
%     
%     r_pmon(i) = sum(r_p(max(i-Delta+1,1):i));
%     r_SPYmon(i) = sum(r_SPY(max(i-Delta+1,1):i));
%     
% %     r_pmon(i) = r_p(i);
% %     r_SPYmon(i) = r_SPY(i);
% end
% 
% per = 252; %Lookback period to estimate mean and stdev of daily returns
% for i=1:N-1
%      imin = max(i-252,1); 
%      
%      mu_p(i) = mean(r_pmon(imin:i));
%      sigma_p(i) = std(r_pmon(imin:i));
%      
%      mu_SPY(i) = mean(r_SPYmon(imin:i));
%      sigma_SPY(i) = std(r_SPYmon(imin:i));
%      
%      Sharpe_p(i) = sqrt(252)*mu_p(i)./sigma_p(i);
%      Sharpe_SPY(i) = sqrt(252)*mu_SPY(i)./sigma_SPY(i);     
%      
% end
% 
% fprintf('Portfolio Sharpe > SPY Sharpe %d, and Sharpe correlation is %d\n',...
%     mean(Sharpe_p(per:N-1)>Sharpe_SPY(per:N-1)), corr(Sharpe_p(per:N-1),Sharpe_SPY(per:N-1)))
% 
% [quantile(Sharpe_p,[0.05;0.25;0.5;0.75;0.95]),quantile(Sharpe_SPY,[0.05;0.25;0.5;0.75;0.95])]
% 
% % Sharpe_p_out = rmoutliers(Sharpe_p,'percentiles',[2.5,97.5]);
% % Sharpe_SPY_out = rmoutliers(Sharpe_SPY,'percentiles',[2.5,97.5]);
% %[quantile(Sharpe_p_out(2:end),[0;0.25;0.5;0.75;1]),quantile(Sharpe_SPY_out(2:end),[0;0.25;0.5;0.75;1])]
% 
% figure 
% hold on
% grid on
% box on
% plot(dates(per:N-1),Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1))
% %plot(dates(2*252:N-1),(mup(2*252:N-1,1)-mun(2*252:N-1,1))./sqrt(sigmap(2*252:N-1,1).^2+sigman(2*252:N-1,1).^2))
% str=strcat('AnnualSharpeMBGDelta');
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off
% 
% figure 
% hold on
% grid on
% box on
% plot(dates(per:N-1),Sharpe_p(per:N-1,1))
% plot(dates(per:N-1),Sharpe_SPY(per:N-1,1))
% %plot(dates(2*252:N-1),(mup(2*252:N-1,1)-mun(2*252:N-1,1))./sqrt(sigmap(2*252:N-1,1).^2+sigman(2*252:N-1,1).^2))
% legend('Portfolio','SPY')
% str=strcat('AnnualSharpeMBG');
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off
% 
% figure 
% hold on
% grid on
% box on
% plot(dates(per:N-1),Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1))
% plot(dates(per:N-1),DeltaBG)
% legend('MBG','BG')
% str=strcat('AnnualSharpeMBGBG');
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off
% 
% figure 
% hold on
% grid on
% box on
% plot(dates(per:N-1),r_pmon(per:N-1),'LineWidth', 2)
% plot(dates(per:N-1),r_SPYmon(per:N-1),'-.','LineWidth', 0.5)
% legend('Portfolio','SPY')
% title('Daily Returns')
% str=strcat('AnnualReturnMBG');
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off

%% Save

% save(fullfile(varPath, strcat('thetaMBG_Convex',thetai0),'theta')
%
% DeltaS = Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1);
% save(fullfile(varPath, strcat('DeltaMBGConvex',thetai0),'DeltaS')

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
    indp = (x>0);
    indn = (x<0);
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