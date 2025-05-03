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

% Coefficient of Absolute Risk Aversion
eps = [5];
M = length(eps);

%% Maximization

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
    for t = tmin:tmax
        bp = bpvec(t);
        cp = cpvec(t)*dt;
        bn = bnvec(t);
        cn = cnvec(t)*dt;
        [fBG,FBG,x] = BG(bp,cp,bn,cn); %Fourier inversion
        a = ( (1-bp)^(-cp) )*( (1+bn)^(-cn) );
        if a < 1
            [bn,bp]=deal(bp,bn);
            [cn,cp]=deal(cp,cn);
            [fBG,FBG,x] = BG(bp,cp,bn,cn); %Fourier inversion
            short(t,j) = 1;
        end
        L = zeros(P,1);
        for i = 1:P            
            % Maximization
    
            options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
            
            f = @(theta)H(eps(j),p(end),theta,r,dt,fBG,x);
    
            L(i) = f(p(i)/p(end));
        end
        [~,theta_t] = max(L(:));
        theta(t,j) = theta_t/P;
        fprintf('t = %d, j = %d, theta = %d, bp = %d, a = %d\n',t,j,theta(t,j),bp,a)
    end
    plot(dates(tmin:tmax),theta(tmin:tmax,j))
    leg{j} = strcat('$\eps = $', num2str(eps(j)));
end

legend = legend(leg,'interpreter','latex','location','best');
grid on
box on
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('ConvexPortfolioChoiceCRRA_2020');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

%% Save

save('ConvexPCthetaCRRA_2020','theta')

%% Routines

function [f,F,x] = BG(bp,cp,bn,cn)
N=2^18;
eta=0.5;
lambda = 2*pi/(N*eta);
beta = lambda*N/2;
u=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
x = -beta+lambda*[0:N-1];

phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
f = max(real(fft((1/pi)*exp(1i*u*beta).*phi*eta)),0);
f = f/(sum(f)*lambda);
F = cumsum(f)*lambda;
end

function t = H(eps,P,theta,r,T,f,x)   
    lambda = x(2)-x(1);
    t = sum( ( ( (1-theta)*(exp(r*T))+theta*(exp(x)) ).^(1-eps) ).*f ) * lambda;
    t = P * ( t.^(1/(1-eps)) );
end
