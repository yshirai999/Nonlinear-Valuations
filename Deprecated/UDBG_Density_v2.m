%% Distorted Bilateral Gamma Levy and Probability Densities
clear
clc
close all

% DBG Parameters
c = 0.5;
gamma = 0.25;
bp = 4.466081e-03;
cp = 1.072381e+02;
bn = 5.256185e-02;
cn = 4.187231e+00;

% Montecarlo Integration Parameters
Nsim = 10000; rng('default'); rng(1);
U = rand(Nsim,1);

% Levy Density

xp = [0.0001:0.0001:0.01]; xn = [-0.01:0.0001:-0.0001]; x = [xn';xp'];

k = [-cn*exp(bn*xn)./xn, cp*exp(-bp*xp)./xp]';

psi_U = [cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn,...
        (1/(1+gamma)) .* ( (1 - exp(-c*cp*expint(xp/bp))).^(-gamma/(1+gamma)) )...
                      .*  cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp]';

psi_L = [- (1/(1+gamma)) .* ( (1 - exp(-c*cn*expint(-xn/bn))).^(-gamma/(1+gamma)) )...
                         .*  cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn,...
                    - cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp]';

figure
hold on
grid on
box on
plot(x,k,'-');
plot(x,k+psi_U,'--');
plot(x,k+psi_L,'-.')
legend('Non-Distorted','Upper','Lower','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('DBG_Levy');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

psi_Un = cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn;
psi_Up = (1/(1+gamma)) .* ( (1 - exp(-c*cp*expint(xp/bp))).^(-gamma/(1+gamma)) )...
                       .*  cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp;

psi_Ln = - (1/(1+gamma)) .* ( (1 - exp(-c*cn*expint(-xn/bn))).^(-gamma/(1+gamma)) )...
                          .*  cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn;
psi_Lp = - cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp;

figure
hold on
grid on
box on
plot(xp,psi_Up,'b');
plot(xn,psi_Un,'b');
legend('$\psi^U$','','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('Psi_U');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

figure
hold on
grid on
box on
plot(xp,psi_Lp,'b');
plot(xn,psi_Ln,'b');
legend('$\psi^L$','','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('Psi_L');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

%% Probability Density

B = 1000000;
N=2^15;
eta=B/N;
lambda = 2*pi/(eta*N);
beta = lambda*N/2;
u=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
x = -beta+lambda*[0:N-1];

yp = -log(U)*bp;
yn = -log(U)*bn;

phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
% phi = w .* exp( sum(...
%             (exp(1i*yp*u)-1)...
%                 .*  bp*cp ./ yp...
%             +(exp(-1i*yn*u)-1) .*  bn*cn ./ yn ) /Nsim);
f = max(real(fft((1/pi)*exp(1i*u*beta).*phi*eta)),0);

phi_u = phi .* exp( -1i*u .* GhatU(bp,cp,bn,cn,c,gamma,1/c,1,u) );
f_u = max(real(fft((1/pi)*exp(1i*u*beta).*phi_u*eta)),0);

phi_l = phi .* exp( -1i*u .* GhatL(bp,cp,bn,cn,c,gamma,1/c,1,u) );
f_l = max(real(fft((1/pi)*exp(1i*u*beta).*phi_l*eta)),0);

figure
hold on
grid on
box on
[~,Nmax] = min((x-0.1).^2);
[~,Nmin] = min((x+0.1).^2);
plot(x(Nmin:Nmax),f(Nmin:Nmax),'-');
plot(x(Nmin:Nmax),f_u(Nmin:Nmax),'--');
plot(x(Nmin:Nmax),f_l(Nmin:Nmax),'-.')
legend('Non-Distorted','Upper','Lower','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('DBG_densities');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

% %% Drift
% 
% omega = ((1-bp).^(-cp)).*((1+bn).^(-cn));
% omegaU = omega * exp( sum(...
%             (exp(yp)-1) * (1/(1+gamma))...
%                 .* ( (1 - exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma)) )...
%                 .*  cp .* exp(-c*cp*expint(yp/bp)-yp/bp) ./ yp...
%             +(exp(-yn)-1) .*  cn .* exp(-c*cn*expint(yn/bn)-yn/bn) ./ yn) /Nsim);
% omegaL = omega * exp( sum(...
%             (exp(yp)-1) .*  cp .* exp(-c*cp*expint(yp/bp)-yp/bp) ./ yp...
%             +(exp(-yn)-1) * (1/(1+gamma))...
%                 .* ( (1 - exp(-c*cn*expint(yn/bn))).^(-gamma/(1+gamma)) )...
%                 .*  cn .* exp(-c*cn*expint(yn/bn)-yn/bn) ./ yn) /Nsim);
% 
% fprintf('[omega^U,omega,omega^L] = [%d,%d,%d]',omegaU,omega,omegaL)

%% Routines

function Ghat = GhatU(bp,cp,bn,cn,c,gamma,a,b,u)

N = 2^16;
eta = 0.0001;
lambda = 2*pi/(eta*N);
beta = lambda*N/2;
y=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
xi = -beta+lambda*[0:N-1];

Gn = [b/c, (b/c) .* (1 - exp(-c*cn*expint(y(2:end)/bn)))];
Gp = [a, a .* (1 - exp(-c*cp*expint(y(2:end)/bp))).^(1/(1+gamma))];
Ghat = fft((1/(2*pi))*exp(1i*y*beta).*Gp.*w*eta)-fft((1/(2*pi))*exp(1i*y*beta).*Gn.*w*eta);
Ghat = interp1(xi,Ghat,u);
figure
plot(u,abs(1i*u.*Ghat))


% G = @(y) exp(1i*y*u).*(a*(1-exp(-c*cp*expint(y/bp))).^(1/(1+gamma)).*(y>0) +...
%          (b/c)*(1-exp(-c*cn*expint(y/bn))).*(y<0));
% 
% GHat = integral(real(G),-Inf,Inf,'ArrayValued',true)+integral(imag(G),-Inf,Inf,'ArrayValued',true);
% plot(u,GHat)

end

function Ghat = GhatL(bp,cp,bn,cn,c,gamma,a,b,u)

N = 2^16;
eta = 0.0001;
lambda = 2*pi/(eta*N);
beta = lambda*N/2;
y=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
xi = -beta+lambda*[0:N-1];

Gn = [a, a .* (1 - exp(-c*cn*expint(y(2:end)/bn))).^(1/(1+gamma))];
Gp = [b/c, (b/c) .* (1 - exp(-c*cp*expint(y(2:end)/bp)))];
Ghat = fft((1/(2*pi))*exp(1i*y*beta).*Gp.*w*eta)-fft((1/(2*pi))*exp(1i*y*beta).*Gn.*w*eta);
Ghat = interp1(xi,Ghat,u);
figure
plot(u,abs(1i*u.*Ghat))

end