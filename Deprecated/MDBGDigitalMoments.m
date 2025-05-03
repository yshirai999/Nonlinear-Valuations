%% MDBG parameters
clear
clc
close all

%% Load Data
Y = load('Y');
Y = Y.Y;

SY = 20; startdate = strcat('20',num2str(SY,'%02.f'),'0110');
SD = str2double(startdate);

N = 250;

options = optimset('MaxFunEval',5000,'MaxIter',1000,'Display','off','TolFun',1e-10,'TolX',1e-10);
optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',1000,'Display','iter','FunctionTolerance',1e-10);

%% Generalized Method of Moments
imin = find(Y(:,1)==SD);
imax = length(Y);
n = length(imin:1:imax);
theta = zeros(n,8);
mu_mod = zeros(n,3);
mu_mkt = zeros(n,3);
RC = zeros(n,2);

d = Y(imin:5:imax,1);

ci = 1;
gammai = 0.2;
bi = 1;
ai = 1;
bp = Y(imin,5); cp = Y(imin,6); bn = Y(imin,7); cn = Y(imin,8);
for i = imin:5:imax
    
    theta0 = [ci,gammai,bi,ai,bp,cp,bn,cn];
    ind = [5:5:N];
    NN = length(ind);
    YY = zeros(NN+1,4);
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
    m = ( log(YY(2:NN+1,4)') - log(YY(1:NN,4)') )/5;%( log(Y(i-N+1:i,4)') - log(Y(i-N:i-1,4)') );
    
    obj = @(theta) objU(u,theta,NN)+objL(l,theta,NN);
    
    [theta(i-imin+1,:)] = abs(fminunc(obj,theta0,options));
%     [theta(i-imin+1,:)] = patternsearch(obj,theta0,[],[],[],[],[],[],[],optionspsrch);
    
    bp = abs(theta(i-imin+1,5));
    cp = abs(theta(i-imin+1,6));
    bn = abs(theta(i-imin+1,7));
    cn = abs(theta(i-imin+1,8));
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

theta = abs(theta(1:5:n,:));
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
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('SPY_RCvsMktReturn_DigitalMoments_MDBG');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off


%fprintf('corr_u = %d, corr_l = %d\n', corr(mu_mod(:,1),mu_mkt(:,1)), corr(mu_mod(:,3),mu_mkt(:,3)));
fprintf('corr_u = %d, corr_l = %d\n', corr(RC(:,1),mu_mkt(:,2)), corr(RC(:,2),mu_mkt(:,3)));

%% Save
SPYMDBG = abs(theta);
save(strcat('SPYMDBG',num2str(SY),'DM'),'SPYMDBG');
save(strcat('RCMDBG',num2str(SY),'MD'),'RC');
save(strcat('mu_modMDBG',num2str(SY)),'mu_mod');
save(strcat('mu_mkt',num2str(SY)),'mu_mkt');

%% Routines

function I = RCU(theta)
    c = abs(theta(1)); gamma = abs(theta(2));
    b = abs(theta(3));
    a = abs(theta(4));
    bp = abs(theta(5)); cp = abs(theta(6)); bn = abs(theta(7)); cn = abs(theta(8)); 
    
    L = @(w) cp*expint(log(1+w)/bp);
    K = @(w) cn*expint(-log(1-w)/bn);
    Gp = @(y) a*(1-exp(-c*y)).^(1/(1+gamma));
    Gm = @(y) (b/c)*(1-exp(-c*y));
    funp = @(w) Gp(L(w));
    funm = @(w) Gm(K(w));

    Ip = integral(funm,0,1);
    Im = integral(funp,0,Inf);
    I = Ip+Im;
end

function I = RCL(theta)
    c = abs(theta(1)); gamma =  abs(theta(2));
    b = abs(theta(3)); 
    a = abs(theta(4));
    bp = abs(theta(5)); cp = abs(theta(6)); bn = abs(theta(7)); cn = abs(theta(8)); 
    
    L = @(w) cp*expint(log(1+w)/bp);
    K = @(w) cn*expint(-log(1-w)/bn);
    Gp = @(y) a*(1-exp(-c*y)).^(1/(1+gamma));
    Gm = @(y) (b/c)*(1-exp(-c*y));
    funp = @(w) Gm(L(w));
    funm = @(w) Gp(K(w));

    Ip = integral(funm,0,1);
    Im = integral(funp,0,Inf);
    I = Ip+Im;
end

function T = objU(u,theta,N)
    s = u + RCU(theta);
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

function T = objL(l,theta,N)
    s = l - RCL(theta);
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