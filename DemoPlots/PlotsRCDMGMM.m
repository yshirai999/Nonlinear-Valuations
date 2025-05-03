%% Plots
clear
clc
close all

%% Load Data
SY = 8;
d = load(strcat('d',num2str(SY)));
RCDM = load(strcat('RC',num2str(SY),'DM'));
RCGMM = load(strcat('RC',num2str(SY),'GMM'));
mu_mod = load(strcat('mu_mod',num2str(SY)));
mu_mkt = load(strcat('mu_mkt',num2str(SY)));


%% Visualization
mu_mod = mu_mod.mu_mod;
mu_mkt = mu_mkt.mu_mkt;
d = d.d;

RC = RCDM.RC;
figure
hold on
box on
grid on
plot(d,RC(:,1))
plot(d,RC(:,2))
plot(d,mu_mkt(:,2))
legend('$RC^U$','$RC^L$','$\mu^Y$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('SPY_RCvsMktReturn_DM');
fname=str;
%saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

% figure
% hold on
% box on
% grid on
% plot(d,RC(:,1))
% plot(d,RC(:,2))
% plot(d,mu_mod(:,2))
% legend('$RC^U$','$RC^L$','$\mu$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
% str=strcat('SPY_RCvsModReturn_DM');
% fname=str;
% saveas(gcf, fullfile(fpath, fname), 'epsc');
% hold off
% 
% figure
% hold on
% box on
% grid on
% plot(d,RC(:,1))
% plot(d,RC(:,2))
% plot(d,mu_mkt(:,2))
% plot(d,mu_mod(:,2))
% legend('$RC^U$','$RC^L$','$\mu$ (market)','$\mu$ (model)','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
% str=strcat('SPY_RCvsMktModReturn_DM');
% fname=str;
% saveas(gcf, fullfile(fpath, fname), 'epsc');
% hold off

fprintf('DM, 2008-2020: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(:,1),mu_mkt(:,2)), corr(RC(:,2),mu_mkt(:,2)));
%fprintf('DM, 2008-2020: corr_u_mod = %d, corr_l_mod = %d\n', corr(RC(:,1),mu_mod(:,2)), corr(RC(:,2),mu_mod(:,2)));
fprintf('DM, Aug-Oct 2008: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(31:41,1),mu_mkt(31:41,2)), corr(RC(31:41,2),mu_mkt(31:41,2)));
fprintf('DM, 2008: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(1:50,1),mu_mkt(1:50,2)), corr(RC(1:50,2),mu_mkt(1:50,2)));
fprintf('DM, 2008-2009: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(1:101,1),mu_mkt(1:101,2)), corr(RC(1:101,2),mu_mkt(1:101,2)));
fprintf('DM, 2008-2010: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(1:151,1),mu_mkt(1:151,2)), corr(RC(1:151,2),mu_mkt(1:151,2)));
fprintf('DM, Feb-Mar 2020: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(607:614,1),mu_mkt(607:614,2)), corr(RC(607:614,2),mu_mkt(607:614,2)));

RC = RCGMM.RC;
figure
hold on
box on
grid on
plot(d,RC(:,1))
plot(d,RC(:,2))
plot(d,mu_mkt(:,2))
legend('$RC^U$','$RC^L$','$\mu^Y$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('SPY_RCvsMktReturn_GMM');
fname=str;
%saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

% figure
% hold on
% box on
% grid on
% plot(d,RC(:,1))
% plot(d,RC(:,2))
% plot(d,mu_mod(:,2))
% legend('$RC^U$','$RC^L$','$\mu$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
% str=strcat('SPY_RCvsModReturn_GMM');
% fname=str;
% saveas(gcf, fullfile(fpath, fname), 'epsc');
% hold off
% 
% figure
% hold on
% box on
% grid on
% plot(d,RC(:,1))
% plot(d,RC(:,2))
% plot(d,mu_mkt(:,2))
% plot(d,mu_mod(:,2))
% legend('$RC^U$','$RC^L$','$\mu$ (market)','$\mu$ (model)','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
% str=strcat('SPY_RCvsMktModReturn_GMM');
% fname=str;
% saveas(gcf, fullfile(fpath, fname), 'epsc');
% hold off

fprintf('GMM, 2008-2020: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(:,1),mu_mkt(:,2)), corr(RC(:,2),mu_mkt(:,2)));
%fprintf('GMM, 2008-2020: corr_u_mod = %d, corr_l_mod = %d\n', corr(RC(:,1),mu_mod(:,2)), corr(RC(:,2),mu_mod(:,2)));
fprintf('GMM, Aug-Oct 2008: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(31:41,1),mu_mkt(31:41,2)), corr(RC(31:41,2),mu_mkt(31:41,2)));
fprintf('GMM, 2008-2009: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(1:101,1),mu_mkt(1:101,2)), corr(RC(1:101,2),mu_mkt(1:101,2)));
fprintf('GMM, 2008-2010: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(1:151,1),mu_mkt(1:151,2)), corr(RC(1:151,2),mu_mkt(1:151,2)));
fprintf('GMM, Feb-Mar 2020: corr_u_mkt = %d, corr_l_mkt = %d\n', corr(RC(607:614,1),mu_mkt(607:614,2)), corr(RC(607:614,2),mu_mkt(607:614,2)));