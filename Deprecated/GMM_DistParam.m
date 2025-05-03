%% Distortion Parameter Estimation

clear
clc

%% Load Data and Set Inputs

Y = load('Y');
Y = Y.Y;

N = 100;
rng('default'); rng(1);
U = rand(1,N);

tau = 250;

options = optimset('MaxFunEval',5000,'MaxIter',1000,'Display','iter','TolFun',1e-10,'TolX',1e-10);

%% Generalized Method of Moments

imin = find(Y(:,1)==20080103);
[n,~] = size(Y(imin:end,:));
theta = ones(n,2).*[0.01,0.2];
for i = imin:imin
    bp = Y(i,5); cp = Y(i,6); bn = Y(i,7); cn = Y(i,8);
    mu = log ( ( (1-bp)^(-cp) ) * ( (1+bn)^(-cn) ) );
    dY = ( log(Y(i,4)) - log(Y(i-tau,4)) ) / tau;
    dU = ( log(Y(i,3)) - log(Y(i-tau,3)) ) / tau;
    dL = ( log(Y(i,2)) - log(Y(i-tau,2)) ) / tau;
    obj = @(theta) (dU-mu+RCU(bp,cp,bn,cn,abs(theta(1)),abs(theta(2)),U,N))^2+(dL-mu-RCL(bp,cp,bn,cn,abs(theta(1)),abs(theta(2)),U,N))^2;
    [theta(i-imin+1,:)] = fminunc(obj,theta(i-imin+1,:),options);
    theta(i-imin+1,:) = [abs(theta(i-imin+1,1)),abs(theta(i-imin+1,2))];
    ci = theta(i-imin+1,1);
    gammai = theta(i-imin+1,2);
    RCi = [RCU(bp,cp,bn,cn,abs(theta(i-imin+1,1)),abs(theta(i-imin+1,2)),U,N), RCL(bp,cp,bn,cn,abs(theta(i-imin+1,1)),abs(theta(i-imin+1,2)),U,N)];
end

%% Routines

function I = RCU(bp,cp,bn,cn,c,gamma,U,N)
    yp = sort(-log(U)*bp);
    yn = sort(-log(U)*bn);
    I = 0;
    for i = 1:N
        I = I + ( (cp*bp)*(exp(yp(i))-1)*(1/(1+gamma))*(1-exp(-c*cp*expint(yp(i)/bp)))*(exp(-c*cp*expint(yp(i)/bp)))...
            +(cn*bn)*(1-exp(-yn(i)))*(exp(-c*cn*expint(yn(i)/bn))) )/N;
    end
end

function I = RCL(bp,cp,bn,cn,c,gamma,U,N)
    yp = sort(-log(U)*bp);
    yn = sort(-log(U)*bn);
    I = 0;
    for i = 1:N
        I = I + ( (cp*bp)*(exp(yp(i))-1)*(exp(-c*cp*expint(yp(i)/bp)))...
            +(cn*bn)*(1-exp(-yn(i)))*(1/(1+gamma))*(1-exp(-c*cn*expint(yn(i)/bn)))*(exp(-c*cn*expint(yn(i)/bn))) )/N;
    end
end