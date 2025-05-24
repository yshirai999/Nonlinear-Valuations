%% Portfolio Choice
clear
clc
close all

thetai0 = 'prev'; %% 'unif' for uniform initial condition, 'prev' for previous

%% Testing Parameters.
% MBG parameters
ticker = {'SPY', 'VIX', 'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY'};
ETF = matlab.lang.makeValidName(ticker);
D = length(ticker);

dataPath = getPath('Data');

for d=1:D
    Y = load(fullfile(dataPath,strcat('Y',ticker{d},'.mat')));
    Y = Y.Y;
    ind = (Y(:,1)>20151231);
    eval([ETF{d} ' = Y(ind,2:end);']);
end

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
    
% Distortion Parameters
c = 0.01;
gamma = 0.25;
aw = 100;
bw = 0.1;

%% Maximization

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
    theta = load(strcat('thetaBG',thetai0,'.mat'));
    theta = theta.theta;
catch
    for i = 1:N-1
        % MBG Parameters
    
        bp = bpvec(i,:)';
        cp = cpvec(i,:)';
        bn = bnvec(i,:)';
        cn = cnvec(i,:)';
    
        a = log( ((1./(1-bp)).^cp).*((1./(1+bn)).^cn) );
    
        % Maximization
    
        options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
        f = @(theta)-H([theta;1-sum(theta)],c,gamma,aw,bw,bp,cp,bn,cn,a,D);
    
    %      thetaplot = 0.01:0.01:0.99;
    %      ff = zeros(size(thetaplot));
    %      for ii = 1:length(thetaplot)
    %          ff(ii) = f(thetaplot(ii));
    %      end
    %     figure
    %     plot(thetaplot,-ff);
    %     
    %      [~,ii] = min(ff);
    %      thetai = thetaplot(ii);
    
        if strcmp(thetai0,'unif')
            thetai = ones(D-1,1)/D;
        elseif strcmp(thetai0,'prev')
        else
            fprintf('Warning: initial conditions not correctly set')
        end
    
        [thetai,Sopt,exitflagBG] = fminsearch(f,thetai,options); thetai = thetai-floor(thetai);
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

vizPath = getPath('Visualization');

close all

r_p = zeros(N,1);
r_SPY = zeros(N,1);

r_pmon = zeros(N,1);
r_SPYmon = zeros(N,1);

r_pmon_down = zeros(N,1);
r_SPYmon_down = zeros(N,1);

mu_p = zeros(N,1);
mu_SPY = zeros(N,1);
sigma_p = zeros(N,1);
sigma_SPY = zeros(N,1);

sigma_p_down = zeros(N,1);
sigma_SPY_down = zeros(N,1);

Sharpe_p = zeros(N,1);
Sharpe_SPY = zeros(N,1);

Sortino_p = zeros(N,1);
Sortino_SPY = zeros(N,1);

Delta = 1;

for i=1:N-1
    r_p(i) = theta(i,:)*(Y(i+1,:)'-Y(i,:)');
    r_SPY(i) = (Y(i+1,1)'-Y(i,1)');
    
    r_pmon(i) = sum(r_p(max(i-Delta+1,1):i));
    r_SPYmon(i) = sum(r_SPY(max(i-Delta+1,1):i));

    r_pmon_down(i) = min(r_pmon(i),0);
    r_SPYmon_down(i) = min(r_SPYmon(i),0);
    %     
%     r_pmon(i) = r_p(i);
%     r_SPYmon(i) = r_SPY(i);
end


per = 252; %Lookback period to estimate mean and stdev of daily returns
for i=1:N-1
     imin = max(i-252,1); 
     
     mu_p(i) = mean(r_pmon(imin:i));
     sigma_p(i) = std(r_pmon(imin:i));
     sigma_p_down(i) = std(r_pmon_down(imin:i));
     
     mu_SPY(i) = mean(r_SPYmon(imin:i));
     sigma_SPY(i) = std(r_SPYmon(imin:i));
     sigma_SPY_down(i) = std(r_SPYmon_down(imin:i));
     
     Sharpe_p(i) = sqrt(252)*mu_p(i)./sigma_p(i);
     Sharpe_SPY(i) = sqrt(252)*mu_SPY(i)./sigma_SPY(i);     
     
     Sortino_p(i) = sqrt(252)*mu_p(i)./sigma_p_down(i);
     Sortino_SPY(i) = sqrt(252)*mu_SPY(i)./sigma_SPY_down(i);     
end

fprintf('Portfolio Sharpe > SPY Sharpe %d, and Sharpe correlation is %d\n',...
    mean(Sharpe_p(per:N-1)>Sharpe_SPY(per:N-1)), corr(Sharpe_p(per:N-1),Sharpe_SPY(per:N-1)))

fprintf('Portfolio Sortino > SPY Sortino %d, and Sortino correlation is %d\n',...
    mean(Sortino_p(per:N-1)>Sortino_SPY(per:N-1)), corr(Sortino_p(per:N-1),Sortino_SPY(per:N-1)))

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
title('Annualized Daily Sharpe Ratio')
str=strcat('AnnualSharpeBG');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');
hold off

figure 
hold on
grid on
box on
plot(dates(per:N-1),Sortino_p(per:N-1,1))
plot(dates(per:N-1),Sortino_SPY(per:N-1,1))
%plot(dates(2*252:N-1),(mup(2*252:N-1,1)-mun(2*252:N-1,1))./sqrt(sigmap(2*252:N-1,1).^2+sigman(2*252:N-1,1).^2))
legend('Portfolio','SPY')
title('Annualized Daily Sortino Ratio')
str=strcat('AnnualSortinoBG');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');
hold off

figure 
hold on
grid on
box on
plot(dates(per:N-1),r_pmon(per:N-1))
plot(dates(per:N-1),r_SPYmon(per:N-1))
legend('Portfolio','SPY')
title('Daily Returns')
str=strcat('AnnualReturnBG');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');
hold off

%% Save

varPath = getPath('VarArchive');

save(fullfile(varPath,strcat('thetaBG',thetai0)), 'theta')

DeltaS = Sharpe_p(per:N-1,1)-Sharpe_SPY(per:N-1,1);
save(fullfile(varPath,'DeltaBG'),'DeltaS');

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

function t = H(theta,c,gamma,aw,bw,bp,cp,bn,cn,a,D)
    theta = theta-floor(theta);
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
    t = theta'*a-sum(t);
end