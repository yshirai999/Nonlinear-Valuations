%% Portfolio Choice
clear
clc
close all

%% Testing Parameters.
% Rate and Term
r = 0.0001;
dt = 252;

% BG parameters
LL = 1;
bp = 6.779746e-03;
cp = 2.877843e+01*dt;
bn = 5.210701e-02;
cn = 6.110235e+00*dt;

bn = 6.779746e-03; %short position
cn = 2.877843e+01*dt;
bp = 5.210701e-02;
cp = 6.110235e+00*dt;

bp = 0.0075;
cp = 1.5592*dt;
bn = 0.0181;
cn = 0.6308*dt;

[fBG,FBG,x] = BG(bp,cp,bn,cn); %Fourier inversion

% Coefficient of Absolute Risk Aversion
eps = [1,2,3];
M = length(eps);

%% Maximization

N = 100;
p_l = 10;
p_u = 1000;
p = linspace(p_l,p_u,N)';
p = [0;p]; N=N+1;

L = zeros(N,M);

strplot = {'-.','--','-'};

figure
hold on
for j=1:M
    for i = 1:N
        fprintf('j=%d, i = %d\n',j,i)
        a = ( (1-bp)^(-cp) )*( (1+bn)^(-cn) );
        if a < 1
            [bn,bp]=deal(bp,bn);
        end
    
        % Maximization
    
        options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
        f = @(theta)H(eps(j),p(end),theta,r,dt,fBG,x);
        L(i,j) = f(p(i)/p(end));
    end
    
    plot(p(:)/p(end),L(:,j),strplot{j})
    [~,theta] = max(L(:,j));
    theta = theta/N;
    leg{j} = strcat('$\epsilon = ', num2str(eps(j)),'$');
end

legend = legend(leg,'interpreter','latex','location','best');
xlabel('$\theta$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\Spectral Martingale Measures');
str=strcat('ConvexPortfolioChoiceCARABG');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

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
    t = sum( ( ( (1-theta)*(exp(r*T))+theta*(exp(x)) ).^(-eps) ).*f ) * lambda;
    t = log(P)-(1/eps)*log(t);
end
