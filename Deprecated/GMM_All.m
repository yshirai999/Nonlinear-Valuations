%% Distortion Parameter Estimation

clear
clc

%% Load Data and Set Inputs

Y = load('Y');
Y = Y.Y;

N = 100;
rng('default'); rng(1);
unif = rand(1,N);

tau = 250;

options = optimset('MaxFunEval',5000,'MaxIter',1000,'Display','iter','TolFun',1e-10,'TolX',1e-10);
optionspsrch = optimoptions('patternsearch','MaxFunctionEvaluation',5000,'MaxIterations',1000,'Display','iter','FunctionTolerance',1e-10);
%% Generalized Method of Moments

imin = find(Y(:,1)==20080103);
[n,~] = size(Y(imin:end,:));
theta = zeros(n,6);
for i = imin:imin
    theta0 = [0.5,0.2,Y(i,5:8)];
    %Y = [ones(tau,1),Y(i-tau:i-1,4),Y(i-tau:i-1,4).^2] / tau;
    U = [ones(tau,1),Y(i-tau:i-1,3),Y(i-tau:i-1,3).^2,Y(i-tau:i-1,3).^3,...
                     Y(i-tau:i-1,3).^4,Y(i-tau:i-1,3).^5,Y(i-tau:i-1,3).^6] / tau;
    L = [ones(tau,1),Y(i-tau:i-1,2),Y(i-tau:i-1,2).^2,Y(i-tau:i-1,2).^3,...
                     Y(i-tau:i-1,2).^4,Y(i-tau:i-1,2).^5,Y(i-tau:i-1,2).^6] / tau;
    
    %dY = ( log(Y(i-tau+1:i,4)') - log(Y(i-tau:i-1,4)') );
    dU = ( log(Y(i-tau+1:i,3)') - log(Y(i-tau:i-1,3)') );
    dL = ( log(Y(i-tau+1:i,2)') - log(Y(i-tau:i-1,2)') );
    
    obj = @(theta) sum( ( (dU-RCU(theta,unif,N))*U ).^2+( (dL-RCL(theta,unif,N))*L ).^2);
    
    [theta(i-imin+1,:)] = fminsearch(obj,theta0,options);
    %[theta(i-imin+1,:)] = patternsearch(obj,theta0,[],[],[],[],[],[],[],optionspsrch);
    ci = abs(theta(i-imin+1,1));
    gammai = abs(theta(i-imin+1,2));
    RCi = [RCU(theta(i-imin+1,:),unif,N), RCL(theta(i-imin+1,:),unif,N)];
end

%% Routines

function I = RCU(theta,unif,N)
    c = abs(theta(1)); gamma = abs(theta(2)); bp = abs(theta(3)); cp = abs(theta(4)); bn = abs(theta(5)); cn = abs(theta(6)); 
    yp = sort(-log(unif)*bp);
    yn = sort(-log(unif)*bn);
    mu = ( (1-bp)^(-cp) ) * ( (1+bn)^(-cn) ) - 1;
    I = 0;
    for i = 1:N
        I = I + ( (cp*bp)*(exp(yp(i))-1)*(1/(1+gamma))*(1-exp(-c*cp*expint(yp(i)/bp)))*(exp(-c*cp*expint(yp(i)/bp)))...
            +(cn*bn)*(1-exp(-yn(i)))*(exp(-c*cn*expint(yn(i)/bn))) )/N;
    end
    I = mu - I;
end

function I = RCL(theta,unif,N)
    c = abs(theta(1)); gamma = abs(theta(2)); bp = abs(theta(3)); cp = abs(theta(4)); bn = abs(theta(5)); cn = abs(theta(6)); 
    yp = sort(-log(unif)*bp);
    yn = sort(-log(unif)*bn);
    I = 0;
    mu = ( (1-bp)^(-cp) ) * ( (1+bn)^(-cn) ) - 1;
    for i = 1:N
        I = I + ( (cp*bp)*(exp(yp(i))-1)*(exp(-c*cp*expint(yp(i)/bp)))...
            +(cn*bn)*(1-exp(-yn(i)))*(1/(1+gamma))*(1-exp(-c*cn*expint(yn(i)/bn)))*(exp(-c*cn*expint(yn(i)/bn))) )/N;
    end
    I = mu + I;
end