%% BG parameters
clear
clc
close all

%% Load Data
Y = load('Y');
Y = Y.Y;

rng('default'); rng(1);
unif = rand(1,100);

N = 250;

options = optimset('MaxFunEval',5000,'MaxIter',1000,'Display','off','TolFun',1e-10,'TolX',1e-10);
optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',1000,'Display','iter','FunctionTolerance',1e-10);

%% Generalized Method of Moments
SY = 20; startdate = strcat('20',num2str(SY),'0103');
SD = str2double(startdate);
imin = find(Y(:,1)==SD); %20080103
imax = length(Y);
[n,~] = size(Y(imin:end,:));
theta = zeros(n,6);
mu_mod = zeros(n,1);
mu_modY = zeros(n,1);
mu_mkt = zeros(n,1);
RCi = zeros(n,2);
ci = 0.5;
gammai = 0.25;

for i = imin:imin
    bp = Y(i,5); cp = Y(i,6); bn = Y(i,7); cn = Y(i,8);
    theta0 = [ci,gammai,bp,cp,bn,cn];
    %theta0 = [0.5,0.25,0.0042, 1.7842, 0.0059, 1.0625];
    ind = [5:5:N];
    NN = length(ind);
    YY = zeros(NN+1,3);
    YY(1,1) = Y(i-N,1);
    YY(1,2) = min(Y(i-N-4:i-N,2));
    YY(1,3) = max(Y(i-N-4:i-N,3));
    YY(1,4) = mean([YY(1,2),YY(1,3)]);
    for j=2:NN+1
        YY(j,1) = Y(i-N+ind(j-1),1);
        YY(j,2) = min(Y(i-N+ind(j-1)-4:i-N+ind(j-1),2));
        YY(j,3) = max(Y(i-N+ind(j-1)-4:i-N+ind(j-1),3));
        YY(j,4) = mean([YY(j,2),YY(j,3)]);
    end
    u = ( log(YY(2:NN+1,3)') - log(YY(1:NN,3)') )/5;
    l = ( log(YY(2:NN+1,2)') - log(YY(1:NN,2)') )/5;
%    u = ( log(Y(i-N+1:i,3)') - log(Y(i-N:i-1,3)') );
%    l = ( log(Y(i-N+1:i,2)') - log(Y(i-N:i-1,2)') );
%    m = ( log(YY(2:NN+1,4)') - log(YY(1:NN,4)') );
    m = ( log(Y(i-N+1:i,4)') - log(Y(i-N:i-1,4)') );
    
    obj = @(theta) objM(m,theta,N,mean(m));
    
    [theta(i-imin+1,3:6),err] = fminunc(obj,theta0(3:end),options);
    %[theta(i-imin+1,:)] = patternsearch(obj,theta0,[],[],[],[],[],[],[],optionspsrch);
    
    bpi = abs(theta(i-imin+1,3));
    cpi = abs(theta(i-imin+1,4));
    bni = abs(theta(i-imin+1,5));
    cni = abs(theta(i-imin+1,6));
    
%     obj = @(theta) objU(u,theta,unif,bp,cp,bn,cn)+objL(l,theta,unif,bp,cp,bn,cn);
%     [theta(i-imin+1,1:2)] = fminsearch(obj,theta0(1:2),options);
%     
%     ci = abs(theta(i-imin+1,1));
%     gammai = abs(theta(i-imin+1,2));
    
%     mu_mod(i) = log(((1-bpi)^(-cpi))*(1+bni)^(-cni)); mu_modY(i) = log(((1-bp)^(-cp))*(1+bn)^(-cn)); mu_mkt(i) = mean(m);
%     RCi(i,:) = [RCU(theta(i-imin+1,:),unif), RCL(theta(i-imin+1,:),unif)];
%     fprintf('Date = %d: mod: mu_u = %d, mu = %d, mu_l = %d\n', Y(i,1), mu_modY(i)-RCi(i,1), mu_modY(i), mu_modY(i)+RCi(i,2) );
%     fprintf('Date = %d: mkt: mu_u = %d, mu = %d, mu_l = %d\n', Y(i,1), mean(u), mean(m), mean(l) );
    fprintf('Date = %d: (bp,cp,bn,cn) = (%d,%d,%d,%d)\n', Y(i,1), bpi, cpi, bni, cni)
    fprintf('Date = %d: (bp,cp,bn,cn) = (%d,%d,%d,%d)\n', Y(i,1), bp, cp, bn, cn)
    fprintf('Error = %d\n', err)
%     fprintf('Date = %d: %d, %d, %d\n', Y(i,1), mu_mod(i), mu_modY(i), mu_mkt(i))
end

[s,Pihat,Pi] = plotM(m,[bp,cp,bn,cn],N);
[si,Pihati,Pii] = plotM(m,[bpi,cpi,bni,cni],N);

figure
hold on
plot(s,Pihat)
plot(s,Pihati)
plot(s,Pi)
legend('$\hat{\pi}$', '$\hat{\pi}_i$', '$\pi$','interpreter','latex')
hold off

%% Visualization
% mu_mod = mu_mod(imin:end);
% mu_modY = mu_modY(imin:end);
% mu_mkt = mu_mkt(imin:end);
% RCi = RCi(imin:end);
d = datetime(Y(imin:1:imax,1),'ConvertFrom','yyyymmdd');
% figure
% plot(d,theta(:,3))
% legend('$b_p$','interpreter','latex')
% hold off
% figure
% plot(d,theta(:,4))
% legend('$b_p$','interpreter','latex')
% hold off
% figure
% plot(d,theta(:,5))
% legend('$b_p$','interpreter','latex')
% hold off
% figure
% plot(d,theta(:,6))
% legend('$b_p$','interpreter','latex')
% hold off

%% Save
% SPYMD = abs(theta(:,1:2));
SPYBG = abs(theta(:,3:6));
%save(strcat('SPYBG',num2str(SY),'DM'),'SPYBG');
%save(strcat('SPYMD',num2str(SY),'LS'),'SPYMD');

%% Routines

function I = RCU(theta,unif)
    N = length(unif);
    c = abs(theta(1)); gamma = abs(theta(2)); bp = abs(theta(3)); cp = abs(theta(4)); bn = abs(theta(5)); cn = abs(theta(6)); 
    yp = sort(-log(unif)*bp);
    yn = sort(-log(unif)*bn);
    I = sum ( (cp*bp)*(1/(1+gamma))*(exp(yp)-1)...
            .* (1-exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma))...
            .* (exp(-c*cp*expint(yp/bp))) ./ yp...
             - (cn*bn)*(exp(-yn)-1)...
            .* (exp(-c*cn*expint(yn/bn))) ./ yn )/N;
end

function I = RCL(theta,unif)
    N = length(unif);
    c = abs(theta(1)); gamma = abs(theta(2)); bp = abs(theta(3)); cp = abs(theta(4)); bn = abs(theta(5)); cn = abs(theta(6)); 
    yp = sort(-log(unif)*bp);
    yn = sort(-log(unif)*bn);
    I = sum ( (cp*bp)*(exp(yp)-1)...
            .* (exp(-c*cp*expint(yp/bp))) ./ yp...
             - (cn*bn)*(1/(1+gamma))*(exp(-yn)-1)...
            .* (1-exp(-c*cn*expint(yn/bn))).^(-gamma/(1+gamma))...
            .* (exp(-c*cn*expint(yn/bn))) ./ yn )/N;
end

function T = objU(u,theta,unif,bp,cp,bn,cn)
    mu = log ( ((1-bp)^(-cp))*((1+bn)^(-cn)) );
    RC = RCU([theta,bp,cp,bn,cn],unif);
    %T = ( norm ( ( u - mu + RC ) * ( u.^([0;1;2;3;4;5;6]) )' ) )^2;
    T =  ( u - mu + RC ) * ( u - mu + RC )';
end

function T = objL(l,theta,unif,bp,cp,bn,cn)
    mu = log ( ((1-bp)^(-cp))*((1+bn)^(-cn)) );
    RC = RCL([theta,bp,cp,bn,cn],unif);
    %T = norm ( ( ( l - mu + RC ) * ( l.^([0;1;2;3;4;5;6]) )' ) )^2;
    T = ( l - mu - RC ) * ( l - mu - RC )';
end

function T = objM(m,theta,N,mu)
    s = m;
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,99,99)/100;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(1); cp = theta(2); bn = theta(3); cn = theta(4);
    B = 10000;
    N=2^13;
    eta=B/N;
    lambda = 2*pi/B;
    b = lambda*N/2;
    u=[0:N-1]*eta;
    w = ones(1,N); w(1)=1/2;
    x = -b+lambda*[0:N-1];
    phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
    f = max(real(fft((1/pi)*exp(1i*u*b).*phi*eta)),0);
    f = [0,interp1(x,f,s)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)];
    Pihat = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
    Pihat(s>0) = 1-Pihat(s>0);
    
    mutheta = log ( ( (1-bp)^(-cp) ) * ( (1+bn)^(-cn) ) );
    
    T = sum( ((Pi-Pihat).^2)./(Pi.*(1-Pi)) );% + ((mu-mutheta)^2)*1e+02;
end

function [s,Pihat,Pi] = plotM(m,theta,N)
    s = m;
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,999,999)/1000;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(1); cp = theta(2); bn = theta(3); cn = theta(4);
    B = 10000;
    N=2^10;
    eta=B/N;
    lambda = 2*pi/B;
    bb = lambda*N/2;
    u=[0:N-1]*eta;
    w = ones(1,N); w(1)=1/2;
    x = -bb+lambda*[0:N-1];
    phi = ((1-1i*u*bp).^(-cp)).*((1+1i*u*bn).^(-cn)).*w;
    f = max(real(fft((1/pi)*exp(1i*u*bb).*phi*eta)),0); plot(x,f)
    f = [0,interp1(x,f,s)]; ds = [s(2)-s(1),s(2:end)-s(1:end-1)]; %plot(s,f(2:end));
    Pihat = cumsum( (f(2:end)+f(1:end-1)) .* ds ) / 2;
    Pihat(s>0) = 1-Pihat(s>0);
    
end
