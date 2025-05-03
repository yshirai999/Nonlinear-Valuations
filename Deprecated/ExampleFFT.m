%% Trial

clear
clc
close all

% DBG Parameters
param0 = [2.530257e-03,1.605663e+02,5.258761e-02,5.922168e+00,...
          1.241251e+00,3.655374e-02,1.575594e-09,1.068307e-01];

bp = param0(1);
cp = param0(2);
bn = param0(3);
cn = param0(4);
c = param0(5);
gamma = param0(6);
a = param0(7);
b = param0(8);

%% FFT
N = 2^12;
eta = 0.15;
lambda = 2*pi/(eta*N);
beta = lambda*N/2;
y=[0:N-1]*eta;
w = ones(1,N); w(1)=1/2;
xi = -beta+lambda*[0:N-1];

Gn = [b/c, (b/c) .* (1 - exp(-c*cn*expint(2*pi*y(2:end)/bn)))];
Gp = [a, a .* (1 - exp(-c*cp*expint(2*pi*y(2:end)/bp))).^(1/(1+gamma))];
GhatFFT = N*(ifft((1/(2*pi))*exp(2*pi*1i*y*beta).*Gp.*w*eta)+ifft((1/(2*pi))*exp(2*pi*1i*y*beta).*Gn.*w*eta));

%% Montecarlo Integration I
% Parameters
Nsim = 100;
q = qrandstream('halton',1,'Skip',1e3,'Leap',1e2);
U = qrand(q,Nsim);
yp = -log(U)*bp;
yn = -log(U)*bn;

% Summation
GhatMI = sum((exp(1i*yp*xi)) * (a*c/(1+gamma))...
                .* ( (1 - exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma)) )...
                .*  bp*cp .* exp(-c*cp*expint(yp/bp)) ./ yp...
            - exp(-1i*yn*xi) .*  b*bn*cn .* exp(-c*cn*expint(yn/bn)) ./ yn) /Nsim;

%% Montecarlo Integration II
% Parameters
Nsim = 100;
q = qrandstream('halton',1,'Skip',1e3,'Leap',1e2);
U = qrand(q,Nsim);
yp = -log(U)*bp;
yn = -log(U)*bn;

% Summation
GhatMII = sum( exp(1i*yp*xi)...
                .* ( (1 - exp(-c*cp*expint(yp/bp))).^(a*c/(1+gamma)) )...
                .* bp .* exp(yp*bp)...
            + exp(-1i*yn*xi) .* b .* ( 1 - exp(-c*cn*expint(yn/bn) ) )...
                .* exp(yn*bn) ) /Nsim;

%% Visualization

figure
subplot(3,3,1)
plot(xi,real(-2*pi*1i*xi.*GhatFFT))
subplot(3,3,2)
plot(xi,imag(-2*pi*1i*xi.*GhatFFT))
subplot(3,3,3)
plot(xi,abs(-2*pi*1i*xi.*GhatFFT))

subplot(3,3,4)
plot(xi,real(GhatMI))
subplot(3,3,5)
plot(xi,imag(GhatMI))
subplot(3,3,6)
plot(xi,abs(GhatMI))

subplot(3,3,7)
plot(xi,real(-1i*xi.*GhatMII))
subplot(3,3,8)
plot(xi,imag(1i*xi.*GhatMII))
subplot(3,3,9)
plot(xi,abs(-1i*xi.*GhatMII))
