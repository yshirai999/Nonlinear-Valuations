%% Digital Moments
clear
clc
close all

%% Load Data
Y = load('Y');
Y = Y.Y;

SY = 08; startdate = strcat('20',num2str(SY,'%02.f'),'0110');
SD = str2double(startdate);
% SPYBG = load(strcat('SPYBG',num2str(SY),'DM'));
% SPYBG = SPYBG.SPYBG;

N = 250;

options = optimset('MaxFunEval',5000,'MaxIter',1000,'Display','off','TolFun',1e-10,'TolX',1e-10);
optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',1000,'Display','off','FunctionTolerance',1e-10);

%% Generalized Method of Moments
imin = find(Y(:,1)==SD);
imax = length(Y);
n = length(imin:1:imax);
theta = zeros(n,8);
mu_mod = zeros(n,3);
mu_mkt = zeros(n,3);
RC = zeros(n,2);
%theta(:,5:8) = SPYBG;
theta(:,5:8) = Y(imin:imax,5:8);
d = Y(imin:5:imax,1);

ci = 1;
gammai = 0.2;
bi = 1;
ai = 1;

for i = imin:5:imax
    theta0 = [ci,gammai,bi,ai];
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
    m = ( log(Y(i-N+1:i,4)') - log(Y(i-N:i-1,4)') );
    
%     bpi = abs(theta(i-imin+1,5));
%     cpi = abs(theta(i-imin+1,6)); 
%     bni = abs(theta(i-imin+1,7));
%     cni = abs(theta(i-imin+1,8)); 
    bp = Y(i,5);
    cp = Y(i,6);
    bn = Y(i,7);
    cn = Y(i,8);
    obj = @(theta) objU(u,[theta,bp,cp,bn,cn],NN)+objL(l,[theta,bp,cp,bn,cn],NN);
    
    [theta(i-imin+1,1:4)] = abs(fminsearch(obj,theta0(1:4),options));
%     [theta(i-imin+1,1:4)] = patternsearch(obj,theta0(1:4),[],[],[],[],[],[],[],optionspsrch);
    
%     ci = abs(theta(i-imin+1,1));
%     gammai = abs(theta(i-imin+1,2));
%     bi = abs(theta(i-imin+1,3));
%     ai = abs(theta(i-imin+1,4));
    
%     mui = ((1-bpi)^(-cpi))*(1+bni)^(-cni);
    mu = ((1-bp)^(-cp))*(1+bn)^(-cn); 
    RC(i-imin+1,:) = [RCU(theta(i-imin+1,:)), RCL(theta(i-imin+1,:))];
    mu_mod(i-imin+1,:) = [log(mu)-RC(i-imin+1,1),log(mu),log(mu)+RC(i-imin+1,2)];
    mu_mkt(i-imin+1,:) = [mean(u),mean(m),mean(l)];
    fprintf('Date = %d: (bp,cp,bn,cn) = (%d,%d,%d,%d), c = %d, gamma = %d, b = %d, a = %d\n', Y(i,1), theta(i-imin+1,5:8), theta(i-imin+1,1:4))
    fprintf('Date = %d: mu_u = %d, mu = %d, mu_l = %d\n', Y(i,1), mu_mod(i-imin+1,1), mu_mod(i-imin+1,2), mu_mod(i-imin+1,3))
    fprintf('Date = %d: mu_u = %d, mu = %d, mu_l = %d\n\n', Y(i,1), mu_mkt(i-imin+1,1), mu_mkt(i-imin+1,2), mu_mkt(i-imin+1,3))
end

theta = abs(theta(1:5:n,1:4));
RC = RC(1:5:n,:);
mu_mod = mu_mod(1:5:n,:);
mu_mkt = mu_mkt(1:5:n,:);
d = datetime(d,'ConvertFrom','yyyymmdd');

%% Visualization
figure
hold on
box on
grid on
plot(d,RC(:,1))
plot(d,RC(:,2))
plot(d,mu_mod(:,2))
legend('$RC^U$','$RC^L$','$\mu$','interpreter','latex')
hold off


%fprintf('corr_u = %d, corr_l = %d\n', corr(mu_mod(:,1),mu_mkt(:,1)), corr(mu_mod(:,3),mu_mkt(:,3)));
fprintf('corr_u = %d, corr_l = %d\n', corr(RC(:,1),mu_mkt(:,2)), corr(RC(:,2),mu_mkt(:,3)));

%% Save
SPYMD = abs(theta(:,1:4));
save(strcat('SPYMD',num2str(SY),'DMv1'),'SPYMD');
save(strcat('RC',num2str(SY),'DMv1'),'RC');
save(strcat('mu_mod',num2str(SY)),'mu_mod');
save(strcat('mu_mkt',num2str(SY)),'mu_mkt');


%% Routines

function I = RCU(theta,unif)
    N = length(unif);
    c = abs(theta(1)); gamma = abs(theta(2));
    b = abs(theta(3));
    a = abs(theta(4));
    bp = abs(theta(5)); cp = abs(theta(6)); bn = abs(theta(7)); cn = abs(theta(8)); 
    yp = sort(-log(unif)*bp);
    yn = sort(-log(unif)*bn);
    I = sum ( (a*c)*(cp*bp)*(1/(1+gamma))*(exp(yp)-1)...
        .* (1-exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma))...
            .* (exp(-c*cp*expint(yp/bp))) ./ yp...
             - (b)*(cn*bn)*(exp(-yn)-1)...
            .* (exp(-c*cn*expint(yn/bn))) ./ yn )/N;
end

function I = RCL(theta,unif)
    N = length(unif);
    c = abs(theta(1)); gamma =  abs(theta(2));
    b = abs(theta(3)); 
    a = abs(theta(4));
    bp = abs(theta(5)); cp = abs(theta(6)); bn = abs(theta(7)); cn = abs(theta(8)); 
    yp = sort(-log(unif)*bp);
    yn = sort(-log(unif)*bn);
    I = sum ( (b)*(cp*bp)*(exp(yp)-1)...
            .* (exp(-c*cp*expint(yp/bp))) ./ yp...
             - (a*c)*(cn*bn)*(1/(1+gamma))*(exp(-yn)-1)...
             .* (1-exp(-c*cn*expint(yn/bn))).^(-gamma/(1+gamma))...
            .* (exp(-c*cn*expint(yn/bn))) ./ yn )/N;
end

function T = objU(u,theta,unif,N)
    s = u + RCU(theta,unif);
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,999,999)/1000;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(5); cp = theta(6); bn = theta(7); cn = theta(8);
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
    
    T = sum( ((Pi-Pihat).^2)./(Pi.*(1-Pi)) );
end

function T = objL(l,theta,unif,N)
    s = l - RCL(theta,unif);
    s = sort(s); pi_N = linspace(1,N,N)/N;
    Pi = linspace(1,999,999)/1000;
    s = interp1(pi_N,s,Pi,'linear','extrap');
    
    Pi(s>0) = 1-Pi(s>0);

    bp = theta(5); cp = theta(6); bn = theta(7); cn = theta(8);
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
    
    T = sum( ((Pi-Pihat).^2)./(Pi.*(1-Pi)) );
end