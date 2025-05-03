% close all
% 
% x = 0.001:0.001:10;
% U = expint(x/4)-expint(x);
% U_approx = exp(-x/4).*log(1+4./x)-0.5*exp(-x).*log(1+2./x);
% 
% figure
% hold on
% box on
% grid on
% plot(x,U)
% plot(x,U_approx)
% legend('x','xapprox')
% 
% b = 0.9;
% c = 1;
% bp = 0.0075;
% bp1 = bp*exp(b/c)-0.001;
% 
% cp = 1.5592;
% psi = -b*exp(-c*cp*expint(x/bp));
% psi1 = exp(-x*(1/bp1-1/bp))-1;
% 
% figure
% hold on
% box on
% grid on
% plot(x,psi)
% plot(x,psi1)
% legend('$\psi$','$\psi_1$','Interpreter','latex')
% 
% figure
% hold on
% box on
% grid on
% plot(x,log(psi1-psi))
% legend('$\Delta$','Interpreter','latex')

%% Distortions on 31 December 2020 Calibrated to Option Prices
close all
Delta = 0.000001;
y = Delta:Delta:0.05;
BG = [0.0039,614.5672,0.0979,3.7175];
MD = [0.0019,1.0328,0.0004,0.0065];
nup = BG(2)*expint(y/BG(1));
num = BG(4)*expint(y/BG(3));
% Gammap = MD(3)*(1-exp(-MD(1)*BG(2)*expint(y/BG(1)))).^(1/(1+MD(2)));
% Gammam = (MD(4)/MD(1))*exp(-MD(1)*BG(4)*expint(y/BG(3)));
% 
% figure
% hold on
% box on
% grid on
% plot(y,(Gammap-Gammam))
% legend('$\Delta$','Interpreter','latex')

%% Gamma_-
Gammap = (MD(4)/MD(1))*exp(-MD(1)*nup);
Gammam = (MD(4)/MD(1))*exp(-MD(1)*num);

fprintf('int_0^{infty} Gamma_-(y)dy = %d\n', sum(Gammap)*Delta)
fprintf('int{-infty}^0 Gamma_-(y)dy = %d\n', sum(Gammam)*Delta)

VizPath=NonlinearPricing.Functions.getPath('Visualization');

figure
hold on
box on
grid on
plot(y,Gammap)
plot(y,Gammam)
legend('$\Lambda_-(\nu([y,\infty)))$','$\Lambda_-(\nu((-\infty,-y]))$','Interpreter','latex')
title("$\Lambda_-$",'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('Gammam_Options');
fname=str;
saveas(gcf, fullfile(VizPath, fname), 'epsc');
hold off

figure
hold on
box on
grid on
plot(y,Gammap-Gammam)
legend('$\Delta$','Interpreter','latex')
title("$\Delta\Lambda_-$",'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('Gammam_Options_Delta');
fname=str;
saveas(gcf, fullfile(VizPath, fname), 'epsc');
hold off

%% Gamma_+
Gammap = MD(3)*(1-exp(-MD(1)*nup)).^(1/(1+MD(2)));
Gammam = MD(3)*(1-exp(-MD(1)*num)).^(1/(1+MD(2)));

fprintf('int_0^{infty} Gamma_+(y)dy = %d\n', sum(Gammap)*Delta)
fprintf('int{-infty}^0 Gamma_+(y)dy = %d\n', sum(Gammam)*Delta)

figure
hold on
box on
grid on
plot(y,Gammap)
plot(y,Gammam)
legend('$\Lambda_+(\nu([y,\infty)))$','$\Lambda_+(\nu((-\infty,-y]))$','Interpreter','latex')
title("$\Lambda_+$",'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('Gammap_Options');
fname=str;
saveas(gcf, fullfile(VizPath, fname), 'epsc');
hold off

figure
hold on
box on
grid on
plot(y,Gammap-Gammam)
legend('$\Delta$','Interpreter','latex')
title("$\Delta\Lambda_+$",'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
str=strcat('Gammap_Options_Delta');
fname=str;
saveas(gcf, fullfile(VizPath, fname), 'epsc');
hold off

%%
close all