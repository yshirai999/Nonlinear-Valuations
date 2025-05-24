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
gammau = 0.75;
gammal = 0.05;

chi = [0.3];
chi2 = 2;
M = length(chi);

%% Maximization

vizPath = getPath('Visualization');

gammaplot = linspace(gammal*1.001,gammau*0.999,50);

P = 100;
p_l = 10;
p_u = 1000;
p = linspace(p_l,p_u,P)';
p = [0;p]; P=P+1;

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
        
        [fBG,FBG,x] = BG(bp,cp,bn,cn); %Fourier inversion
        a = sum(exp(x).*fBG);
        if a < 1
            [bn,bp]=deal(bp,bn);
            [cn,cp]=deal(cp,cn);
            [fBG,FBG,x] = BG(bp,cp,bn,cn); %Fourier inversion
            a = sum(exp(x).*fBG);
            short(t,j) = 1;
        end
        L = zeros(P,1);
        for i = 1:P            
            % Maximization
    
            options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
            
            f = @(gamma)H(gamma,gammal,gammau,alpha,p(end),p(i)/p(end),r,dt,fBG,FBG,x);
    
            fgammaplot = zeros(size(gammaplot));
    
            for gammaind = 1:length(gammaplot)
                fgammaplot(gammaind) = f(gammaplot(gammaind));
            end
            
            [~,gammaind] = min(fgammaplot);
            gamma = gammaplot(gammaind);
            
            L(i) = f(gamma);
    
        end
    
        [~,theta_t] = max(L(:));
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
str=strcat('ConvexPortfolioChoicePD_2020');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

%% Save
varPath = getPath('VarArchive');
save(fullfile(varPath, 'ConvexPCthetaPD_2020'),'theta')

%% Routines

function d = dpsi(x,gamma)
%    d = ( ( 1-x .^ (1/(gamma+1)) ) .^ (1+gamma) )...
%        .* ( x.^(-gamma/(gamma+1)) );
    d = (1+gamma).*( ( 1-x ) .^ (gamma) );
    d(~isfinite(d)) = 0;
    d(isnan(d)) = 0;
end

function [f,F,x] = BG(bp,cp,bn,cn)

N=2^14;
eta=0.15;
lambda = 2*pi/(N*eta);
beta = lambda*N/2;
u=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
x = -beta+lambda*[0:N-1];

phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
f = max(real(fft((1/pi)*exp(1i*u*beta).*phi*eta)),0);
f = f/(sum(f)*lambda);
F = min(cumsum(f)*lambda,1);
end

function f = BGDensity(x,fx,y)
f = interp1(x,fx,y,'linear','extrap');
%f = f.*(y<max(x)).*(y>min(x));
end

function F = BGcdf(x,Fx,y)
F = interp1(x,Fx,y,'linear','extrap');
%F = F.*(y<max(x)).*(y>min(x))+(y>max(x));
end

function b = beta(gamma,gammal,gammau,alpha)
    b = ( exp(alpha(1)*(gamma-gammal/10))-1-alpha(1)*(gamma-gammal/10) ) ./...
            (gammau*1.1-gamma).^(1+alpha(2));
    
    b(~isfinite(b)) = 0;
    b(isnan(b)) = 0;
end


function t = H(gamma,gammal,gammau,alpha,P,theta,r,T,f,F,x)
%     fun = @(y) y .* ...%dpsi( BGcdf(x,F,y) , gamma ).*...
%                     BGDensity(x,f,y);
%                 
%     t = integral(fun,min(x),max(x));
    lambda = x(2)-x(1);
    t = sum((exp(x)-1).*dpsi(F,gamma).*f)*lambda;
    t = (1-theta)*P*(exp(r*T)-1)+P*theta*t+beta(gamma,gammal,gammau,alpha);
end