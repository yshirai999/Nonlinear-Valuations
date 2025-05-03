%% Portfolio Choice

clear
clc
close all

%% Testing Parameters.
% Rate and Term
r = 0.0001;
dt = 252;

% BG parameters
bp = 0.0076;
cp = 1.5592*dt;
bn = 0.0181;
cn = 0.6308*dt;

[fBG,FBG,x] = BG(bp,cp,bn,cn); %Fourier inversion

% Rebate Parameters
gammau = 0.75;
gammal = 0.05;

chi = [0.1,0.2,0.3];
chi2 = 2;
M = length(chi);
gammaplot = linspace(gammal*1.001,gammau*0.999,50);

figure
hold on
for j=1:M
    alpha = [chi(j),chi2];
    b = beta(gammaplot,gammal,gammau,alpha);
    plot(gammaplot,b)
    %leg1{j} = strcat('$\chi = $', num2str(chi(j)));
end
%legend = legend(leg1,'interpreter','latex','location','best');
xlabel('$\theta$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')

%% Maximization

N = 100;
p_l = 10;
p_u = 1000;
p = linspace(p_l,p_u,N)';
p = [0;p]; N=N+1;

b = zeros(N,M);
gamma = zeros(N,M);
L = zeros(N,M);

strplot = {'-.','--','-'};

figure
hold on
for j=1:M
    alpha = [chi(j),chi2];
    for i = 1:N
        fprintf('j = %d,i = %d\n',j,i)
    
        % Maximization
    
        options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
        f = @(gamma)H(gamma,gammal,gammau,alpha,p(end),p(i)/p(end),r,dt,fBG,FBG,x);
    
        fgammaplot = zeros(size(gammaplot));
    
        for gammaind = 1:length(gammaplot)
            fgammaplot(gammaind) = f(gammaplot(gammaind));
        end
    
        [~,gammaind] = min(fgammaplot);
        gamma(i,j) = gammaplot(gammaind);
    
        %[gamma(i,j),L(i,j)] = fminsearch(f,gamma(i,j),options);
        %gamma(i,j) = max(min(gamma(i,j),gammaplot(end)),gammaplot(1));
            
        L(i,j) = f(gamma(i,j));
        b(i,j) = beta(gamma(i,j),gammau,gammal,alpha);
    end
    
    plot(p/p(end),L(:,j),strplot{j})
    [~,theta] = max(L(:,j));
    theta = theta/N;
    leg{j} = strcat('$\chi = $', num2str(chi(j)));
end

legend = legend(leg,'interpreter','latex','location','best');
xlabel('$\theta$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\Spectral Martingale Measures');
str=strcat('ConvexPortfolioChoicePD');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

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
F = cumsum(f)*lambda;
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
