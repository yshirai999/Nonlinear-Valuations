%% Correlation

clear
clc
close all

%% Load data

varPath = getPath('VarArchive');

YPD = load(fullfile(varPath, 'ConvexPCthetaPD_2020'));
YCARA = load(fullfile(varPath, 'ConvexPCthetaCARA_2020'));
YCRRA = load(fullfile(varPath, 'ConvexPCthetaCRRA_2020'));
YMD = load(fullfile(varPath, 'ConvexPCthetaMD_2020'));

YPD = YPD.theta;
YMD = YMD.theta;
YCRRA = YCRRA.theta;
YCARA = YCARA.theta;
rho_CRRACARA = corr(YCRRA,YCARA);
rho_CRRAMD = corr(YCRRA,YMD);
rho_CARAMD = corr(YCARA,YMD);
rho_CARAPD = corr(YCARA,YPD);
rho_CRRAPD = corr(YCRRA,YPD);
rho_MDPD = corr(YPD,YMD);

fprintf('CRRACARA = %d, CRRAMD = %d, CARAMD = %d\n', rho_CRRACARA,rho_CRRAMD,rho_CARAMD)
fprintf('CARAPD = %d, CRRAPD = %d, MDPD = %d\n', rho_CARAPD,rho_CRRAPD,rho_MDPD)
