%% Distorted Bilateral Gamma Levy and Probability Densities
clear
clc
close all

% DBG Parameters
% param0 = [3.294343e-02,1.580335e+00,5.256898e-02,5.982879e+00,...
%           4.121823e-01,1.685782e+08,9.379992e+00,5.716245e-02];
param0 = [0.005,1.544,0.010,0.664,...
          1,0.25,1,1];
% param0 = [6.779746e-03,2.877843e+01,5.210701e-02,6.110235e+00,...
%           1.172443e+02,7.720422e-11,8.529200e-03,1]; %Nsim = 1000
% param0 = [6.779746e-03,2.877843e+01,5.210701e-02,6.110235e+00,...
%           1.171026e+02,1.655907e+00,8.539520e-03,1]; %Nsim = 5000
      
param0 = [6.779746e-03,2.877843e+01,5.210701e-02,6.110235e+00,...
          0.01,0.25,100,0.1]; 

param0 = [0.0075,1.5592,0.0181,0.6308,...
          0.01,0.25,100,0.1];      
      
bp = param0(1);
cp = param0(2);
bn = param0(3);
cn = param0(4);
c = param0(5);
gamma = param0(6);
a = param0(7);
b = param0(8);

T = 0.0833;
T = 1;
%% Levy Density

xp = [0.001:0.001:0.3]; xn = [-0.3:0.001:-0.001]; x = [xn';xp'];

k = [-cn*exp(xn/bn)./xn, cp*exp(-xp/bp)./xp]';

psi_U = [b*cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn,...
        (a*c/(1+gamma)) .* ( (1 - exp(-c*cp*expint(xp/bp))).^(-gamma/(1+gamma)) )...
                        .*  cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp]';

psi_L = [- (a*c/(1+gamma)) .* ( (1 - exp(-c*cn*expint(-xn/bn))).^(-gamma/(1+gamma)) )...
                         .*  cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn,...
                    - b*cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp]';

vizPath = getPath('Visualization');

figure
hold on
grid on
box on
plot(x,log(k),'-',LineWidth=1);
plot(x,log(k+psi_U),'--',LineWidth=1);
plot(x,log(k+psi_L),':',LineWidth=1)
%legend('\bf{Q}','$\overline{\mathbb{Q}}$','$\underline{\mathbb{Q}}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('DBG_Levy');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

% psi_Un = cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn;
% psi_Up = (1/(1+gamma)) .* ( (1 - exp(-c*cp*expint(xp/bp))).^(-gamma/(1+gamma)) )...
%                        .*  cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp;
% 
% psi_Ln = - (1/(1+gamma)) .* ( (1 - exp(-c*cn*expint(-xn/bn))).^(-gamma/(1+gamma)) )...
%                           .*  cn .* exp(-c*cn*expint(-xn/bn)+xn/bn) ./ xn;
% psi_Lp = - cp .* exp(-c*cp*expint(xp/bp)-xp/bp) ./ xp;
% 
% figure
% hold on
% grid on
% box on
% plot(xp,psi_Up,'b');
% plot(xn,psi_Un,'b');
% legend('$\psi^U$','','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% str=strcat('Psi_U');
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off
% 
% figure
% hold on
% grid on
% box on
% plot(xp,psi_Lp,'b');
% plot(xn,psi_Ln,'b');
% legend('$\psi^L$','','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% str=strcat('Psi_L');
% fname=str;
% saveas(gcf, fullfile(vizPath, fname), 'epsc');
% hold off

%% Probability Density

% Montecarlo Integration Parameters
Nsim = 1000;
q = qrandstream('halton',1,'Skip',1e3,'Leap',1e2);
U = qrand(q,Nsim);
yp = -log(U)*bp;
yn = -log(U)*bn;

N=2^16;
eta=0.15;
lambda = 2*pi/(N*eta);
beta = lambda*N/2;
u=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
x = -beta+lambda*[0:N-1];

phi = ((1-1i*u*bp).^(-T*cp)).*((1+1i*u*bn).^(-T*cn)).*w;
% phi = w .* exp( sum(...
%             (exp(1i*yp*u)-1)...
%                 .*  bp*cp ./ yp...
%             +(exp(-1i*yn*u)-1) .*  bn*cn ./ yn ) /Nsim);
f = max(real(fft((1/pi)*exp(1i*u*beta).*phi*eta)),0);

phi_u = phi .* exp( T * sum(...
            (exp(1i*yp*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma)) )...
                .*  bp*cp .* exp(-c*cp*expint(yp/bp)) ./ yp...
            -(exp(-1i*yn*u)-1) .*  b*bn*cn .* exp(-c*cn*expint(yn/bn)) ./ yn) /Nsim);
f_u = max(real(fft((1/pi)*exp(1i*u*beta).*phi_u*eta)),0);

phi_l = phi .* exp( T * sum(...
            -(exp(1i*yp*u)-1) .*  b*bp*cp .* exp(-c*cp*expint(yp/bp)) ./ yp...
            +(exp(-1i*yn*u)-1) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cn*expint(yn/bn))).^(-gamma/(1+gamma)) )...
                .*  bn*cn .* exp(-c*cn*expint(yn/bn)) ./ yn) /Nsim);
f_l = max(real(fft((1/pi)*exp(1i*u*beta).*phi_l*eta)),0);
%%
figure
hold on
grid on
box on
[~,Nmax] = min((x-0.1).^2);
[~,Nmin] = min((x+0.1).^2);
plot(x(Nmin:Nmax),f(Nmin:Nmax),'-',LineWidth=1);
plot(x(Nmin:Nmax),f_u(Nmin:Nmax),'--',LineWidth=1);
plot(x(Nmin:Nmax),f_l(Nmin:Nmax),':',LineWidth=1)
%legend('$\mathbb{Q}$','$\overline{\mathbb{Q}}$','$\underline{\mathbb{Q}}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('DBG_densities');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

figure
hold on
grid on
box on
[~,Nmax] = min((x-0.25).^2);
[~,Nmin] = min((x+0.5).^2);
plot(x(Nmin:Nmax),f_u(Nmin:Nmax)-f(Nmin:Nmax),'-',LineWidth=1);
plot(x(Nmin:Nmax),f_l(Nmin:Nmax)-f(Nmin:Nmax),'-',LineWidth=1);
legend('$f_X^U-f_X$','$f_X^L-f_X$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('DBG_densities_Delta');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
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