%% Portfolio Choice
clear
clc
close all

% This file plots the different rebate functions for the entropic and spectral case.
%% Testing Parameters.
% Rate and Term
r = 0.0005;
T = 10;

% BG parameters
LL = 250;
bp = 7.5964e-03;
cp = 1.5592e+00;
bn = 1.8120e-02;
cn = 6.3083e-01;

mup = bp.*cp;
sigmap = bp.*sqrt(cp);
mun = bn.*cn;
sigman = bn.*sqrt(cn);

vis_c = 0;
vis_p = 0;

% Distortion Parameters
cu = 5;
cl = 0.05;
gamma = 0.25;

% Entropic Rebate Parameters
eps = 2;

% Spectral Rebate Parameters

chi = 4;
chi2 = 0.5;
alpha = [chi,chi2];

%% Maximization

c = linspace(cl*1.5,cu*0.9,1000);
% b = beta(cplot,cu,cl,alpha);
% plot(cplot,b)

N = 1;
p = 100;
theta = 0.5;

b = @(c)H(c,gamma,bp,cp,bn,cn,p,T,eps,theta);
    
ebc = zeros(size(c));
sbc = zeros(size(c));
for i = 1:length(c)
    ebc(i) = b(c(i));
    sbc = beta(c,cu,cl,alpha);
end
    
vizPath = getPath('Visualization');
figure
hold on
plot(c,ebc)
xlabel('$c$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')
str=strcat('EntropicRebate');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');
hold off

figure
hold on
plot(c,sbc)
xlabel('$c$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')
str=strcat('SpectralRebate');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'pdf');
hold off


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

function t = H(c,gamma,bp,cp,bn,cn,P,eps,T,theta)
    funp = @(y) dpsin( cp*expint(y/bp), c ).*...
                    kappa(y,bp,cp,bn,cn);
                
    funn = @(y) dpsip( cn*expint(y/bn), c, gamma).*...
                    kappa(-y,bp,cp,bn,cn);
        
    Ip = integral(funp,0,Inf);
    In = integral(funn,0,Inf);
    t = Ip - In;
    a = log( ( (1-bp).^(-cp) ) .* ( (1+bn).^(-cn) ) );
    
    t = theta*P*T*(t+a)+(theta*P*T/eps)*(cp*log(1+eps*theta*bp)+cn*log(1-eps*theta*bn));
    
end

function b = beta(c,cu,cl,alpha)
    b = alpha(1)*exp( 1./(c-cl).^alpha(2) ).*exp(-1./(cu-c).^alpha(2)).*(c<cu);
    
%     b = ( exp(alpha(1)*(c-cu))-1-alpha(1)*(c-cu) ) ./...
%             (c-cl).^alpha(2);
    
    b(~isfinite(b)) = 0;
    b(isnan(b)) = 0;
end
