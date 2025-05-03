%% Example
clear 
clc
close all


% s = 10;
% x = [-10:0.01:5];
% rng('default')
% theta = rand(1,1); theta = [theta,1-theta];
% %theta = [0.5,0.5];
% 
% f1 = (s+1)/theta(1)-(theta(2)/theta(1))*exp(x);
% f2 = (s+1)/theta(2)-(theta(1)/theta(2))*exp(x);
% ind1 = (f1>0);
% ind2 = (f2>0);
% 
% figure 
% hold on
% plot(x(ind1),log(f1(ind1)))
% plot(log(f2(ind2)),x(ind2))
% legend('f1','f2')
% hold off

%% Sigman

% bp = 6.779746e-03;
% cp = 2.877843e+01;
% bn = 5.210701e-02;
% cn = 6.110235e+00;
% 
% bp = 7.5964e-03;
% cp = 1.5592e+00;
% bn = 1.8120e-02;
% cn = 6.3083e-01;
% 
% N = 100;
% mup = bp*cp;
% sigmap = bp*sqrt(cp);
% mun = bn*cn;
% sigman = linspace(0.005,0.25,N)';
% 
% bp = zeros(N,1);
% cp = zeros(N,1);
% bn = zeros(N,1);
% cn = zeros(N,1);
% 
% for i = 1 : N
%     bp(i) = ( (sigmap.^2)./mup  )';
%     cp(i) = ( (mup.^2)./(sigmap.^2) )';
%     bn(i) = ( (sigman(i).^2)./mun )';
%     cn(i) = ( (mun.^2)./(sigman(i).^2) )';
% end
% 
% a = log( ((1-bp).^(-cp)).*((1+bn).^(-cn)) );
% 
% figure
% plot(sigman,a)
% xlabel('$\sigma_{n}$', 'Interpreter', 'latex')

%% Sigmap

% bp = 6.779746e-03;
% cp = 2.877843e+01;
% bn = 5.210701e-02;
% cn = 6.110235e+00;
% 
% N = 100;
% mup = bp*cp;
% sigmap = linspace(0.005,0.25,N)'; %3.6370e-02
% mun = bn*cn;
% sigman = bn*sqrt(cn);
% 
% bp = zeros(N,1);
% cp = zeros(N,1);
% bn = zeros(N,1);
% cn = zeros(N,1);
% 
% for i = 1 : N
%     bp(i) = ( (sigmap(i).^2)./mup  )';
%     cp(i) = ( (mup.^2)./(sigmap(i).^2) )';
%     bn(i) = ( (sigman.^2)./mun )';
%     cn(i) = ( (mun.^2)./(sigman.^2) )';
% end
% 
% a = log( ((1-bp).^(-cp)).*((1+bn).^(-cn)) );
% 
% figure
% plot(sigmap,a)
% xlabel('$\sigma_{p}$', 'Interpreter', 'latex')

%%
% c = 0.01:0.01:10;
% p = abs(100*rand(1,1));
% Gamman = (1./c).*(1-exp(-c*p));
% figure
% box on
% grid on
% plot(c,Gamman)

%% 
% close all
% clear
% clc
% 
% cu = 0.5;
% cl = 0.005;
% gamma = 0.25;
% alpha = 1;
% 
% c = 0.05:0.001:0.15;
% 
% b = alpha*exp( 1./(c-cl) ).*exp(-1./(cu-c)).*(c<cu);
% b(~isfinite(b)) = 0;
% b(isnan(b)) = 0;
% 
% plot( c, b)

%% 
% close all
% clear
% clc
% 
% cu = 0.99;
% cl = 0.1;
% gamma = 0.25;
% a = 1;
% b = 1;
% 
% c = 0.1:0.1:50;
% 
% x = 0.01:0.01:0.1;
% 
% figure
% hold on
% box on
% grid on
% 
% for i=1:length(x)
%     %reb = exp( 1./(c-cl) ).*exp(-1./(cu-c)).*(c<cu);
%     psip = b.*exp(-c*x(i));%-reb;
%     psin = -((a.*c.*exp(-c*x(i)))/(1+gamma)).*(1-exp(-c*x(i))).^(-gamma/(1+gamma));
%     plot(c,psin)    
% end
% hold off

%%
close all
clear
clc

gamma = 0.2;
c = 0.0001:0.0001:0.1;
f1 = (1./c).*(1-exp(-c));
f2 = (1./c).*(1-exp(-c)).^(1/(1+gamma));

figure
hold on
plot(c,f2)
%legend('$f_1$', '$f_2$', 'interpreter','latex')