%% Portfolio Choice
clear
clc
close all

%% Data and Parameters specification.
% Risk free rate
r = 0.0001;
dt = 252;
% BG parameters
ticker = {'SPY'};
ETF = matlab.lang.makeValidName(ticker);
D = length(ticker);

for d=1:D

Y = load(strcat('Y',ticker{d},'.mat'));
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

% Rebate and Distortion Parameters
cu = 2;
cl = 0.25;
gamma = 0.01;

chi = [1];
chi2 = 2*[1];
M = length(chi);

%% Maximization

cplot = linspace(cl*1.001,cu*0.999,50);
% b = beta(cplot,cu,cl,alpha);
% plot(cplot,b)

P = 100;
p_l = 10;
p_u = 1000;
p = linspace(p_l,p_u,P)';

theta = zeros(N,M);

short = zeros(N,M);
tmin = 1000;
tmax = 1259;

figure
hold on
for j=1:M
    alpha = [chi(j),chi2(j)];
    for t = tmin:tmax
        bp = bpvec(t);
        cp = cpvec(t)*dt;
        bn = bnvec(t);
        cn = cnvec(t)*dt;
        a = log( ( (1-bp).^(-cp) ) .* ( (1+bn).^(-cn) ) );
        if a < 0
            [bn,bp]=deal(bp,bn);
            [cn,cp]=deal(cp,cn);
            a = log( ( (1-bp).^(-cp) ) .* ( (1+bn).^(-cn) ) );
            short(t,j) = 1;
        end
        L = zeros(P,1);
        for i = 1:P            
            % Maximization
    
            options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
            f = @(c)H(c,gamma,bp,cp,bn,cn,cu,cl,alpha,p(i));
    
            fcplot = zeros(size(cplot));
    
            for cind = 1:length(cplot)
                fcplot(cind) = f(cplot(cind));
            end
    
            [~,cind] = min(fcplot);
            c = cplot(cind);
            L(i) = f(c);
    
            %[c,L] = fminsearch(f,c,options);
        end
    
        %plot(p/p(end),(p(end)-p)*dt*r+p*a+L(:))
        [~,theta_t] = max((p(end)-p)*dt*r+p*a+L(:));
        theta(t,j) = theta_t/P;
        fprintf('t = %d, j = %d, theta = %d, bp = %d, a = %d\n',t,j,theta(t,j),bp,a)
    end
    plot(dates(tmin:tmax),theta(tmin:tmax,j))
    leg{j} = strcat('$\chi = $', num2str(chi(j)));
end

legend = legend(leg,'interpreter','latex','location','best');
grid on
box on
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('ConvexPortfolioChoiceMD_2020');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

%% Save

save('ConvexPCthetaMD_2020','theta')

%% Routines

function d = dpsin(x,c)
    d = 0.1*exp( -c * x );
end
 
function d = dpsip(x,c,gamma)
    d = (1/(1+gamma))...
        .* ( ( 1-exp( - c * x ) ) .^ (-gamma/(1+gamma)) )...
        .* exp( -c * x );
    d(~isfinite(d)) = 0;
end

function k = kappa(x,bp,cp,bn,cn)
    k = (exp(x)-1).*...
        ( cp * ( exp(-abs(x)/bp) ./ abs(x) ) .* (x>0)...
        + cn * ( exp(-abs(x)/bn) ./ abs(x) ) .* (x<0));
    k(isnan(k))=0;
end

function t = H(c,gamma,bp,cp,bn,cn,cu,cl,alpha,P)
    funp = @(y) dpsin( cp*expint(y/bp), c ).*...
                    kappa(y,bp,cp,bn,cn);
                
    funn = @(y) dpsip( cn*expint(y/bn), c, gamma).*...
                    kappa(-y,bp,cp,bn,cn);
        
    Ip = integral(funp,0,Inf);
    In = integral(funn,0,Inf);
    t = Ip - In;
    t = -P*t+beta(c,cu,cl,alpha);
end

function b = beta(c,cu,cl,alpha)
    b = alpha(1)*exp( 1./(c-cl).^alpha(2) ).*exp(-1./(cu-c).^alpha(2)).*(c<cu);
    
    b = ( exp(alpha(1)*(c-cu))-1-alpha(1)*(c-cu) ) ./ (c-cl).^alpha(2);
    
    b(~isfinite(b)) = 0;
    b(isnan(b)) = 0;
end