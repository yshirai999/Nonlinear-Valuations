clear
clc

bp = 0.005; 
cp = 1.544;
bn = 0.010;
cn = 0.664;
c = 1;
gamma = 0.25;
b = 1;
a = 1;

N = 1000000;
u = rand(N,1);

y = log(-u);
% Iv2 = sum ( a*(1-exp(-c*cp*expint(log(1+y)/bp))).^(1/(1+gamma)) .* exp(y)...
%              + (b/c)*(1-exp(-c*cn*expint(-log(1-u)/bn))) )/N;

L = @(w) cp*expint(log(1+w)/bp);
K = @(w) cn*expint(-log(1-w)/bn);
Gp = @(y) a*(1-exp(-c*y)).^(1/(1+gamma));
Gm = @(y) (b/c)*(1-exp(-c*y));
funp = @(w) Gp(L(w));
funm = @(w) Gm(K(w));

Ipv2 = integral(funm,0,1,'RelTol',0,'AbsTol',1e-12);
Imv2 = integral(funp,0,Inf,'RelTol',0,'AbsTol',1e-12);
Iv2 = Ipv2+Imv2;

yp = -log(u)*bp;
yn = -log(u)*bn;
Ipv1 = sum ( (a*c)*(bp*cp)*(1/(1+gamma))*(exp(yp)-1)...
            .* (1-exp(-c*cp*expint(yp/bp))).^(-gamma/(1+gamma))...
            .* (exp(-c*cp*expint(yp/bp))) ./ yp )/N;%...
Imv1 = sum ( - (b)*(bn*cn)*(exp(-yn)-1)...
             .* (exp(-c*cn*expint(yn/bn))) ./ yn )/N;
Iv1 = Ipv1+Imv1;

% dGp = @(x) (a*c) * (1/(1+gamma)) * ( (1-exp(-c*x)).^(-gamma/(1+gamma)) ) .* exp(-c*x); %x>0
% dGm = @(x) b*exp(-c*x); %x>0
% nup = @(y) cp*expint(y/bp); %y>0
% num = @(y) cn*expint(-y/bn); %y<0
% dnup = @(y) -cp*exp(-y/bp)./y; %y>0
% dnum = @(y) cn*exp(-abs(y)/bn)./abs(y); %y<0
% funp = @(y) -(exp(y)-1).*dGp(nup(y)).*dnup(y);
% funm = @(y) -(exp(-y)-1).*dGm(num(y)).*dnum(y);
% Iv1 = integral(funp,0,Inf)+integral(funm,-Inf,0);

% y = [0.1:0.1:10];
% Psi = dnup(1)
% plot(y,Psi)

fprintf('I: v1 = %d, v2 = %d\n', Iv1, Iv2)