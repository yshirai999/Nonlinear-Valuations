%% Monetary Base Plot
clear
close all
clc

%% Load and Visualize
XX = load('MonetaryBase.mat');
XX = XX.XX;
dates = XX.BOGMBASE;
M = XX.MonetaryBaseTotalMillionsOfDollarsMonthlyNotSeasonallyAdjusted;

figure
hold on
box on
grid on
plot(dates,M)
title('Total Monetary Base (Source: FRED)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\Spectral Martingale Measures');
str=strcat('MonetaryBase');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

Y = load('Y');
Y = Y.Y;
SY = 10; startdate = strcat('20',num2str(SY,'%02.f'),'0105');
SD = str2double(startdate);
N = 250;
enforcegamma = 1;
enforceb = 1;
Delta = 5;
RC = load(strcat('RC',num2str(SY),'DM',num2str(Delta),num2str(enforcegamma),num2str(enforceb),num2str(N))); RC = RC.RC;
imin = find(Y(:,1)==SD);
imax = length(Y);
d = Y(imin:Delta:imax,1);
d = datetime(d,'ConvertFrom','yyyymmdd');
RC = interp1(d,RC,dates);
RCDelta = RC(:,1)-RC(:,2);
figure
hold on
box on
grid on
plot(dates,RCDelta)
title('Distance between upper and lower DM drivers','Interpreter','latex')


M(1)=M(3);
M(2)=M(3);
MNorm = M/M(1);
RCDelta(1) = RCDelta(4);
RCDelta(2) = RCDelta(4);
RCDelta(3) = RCDelta(4);
RCDeltaNorm = RCDelta/RCDelta(1);


close all
figure
hold on
box on
grid on
title('Monetary Base and Difference between DM drivers (normalized)', 'Interpreter','latex')
plot(dates(3:end),RCDeltaNorm(3:end),'-','linewidth',1)
plot(dates(3:end),MNorm(3:end),':','linewidth',1)
x = [1 1 180 180];
y = [-1.5 3.5 3.5 -1.5];
filler1 = fill(x,y,'r','FaceAlpha',.3);
x = [350 350 600 600];
y = [-1.5 3.5 3.5 -1.5];
filler1 = fill(x,y,'r','FaceAlpha',.3);
x = [1000 1000 1800 1800];
y = [-1.5 3.5 3.5 -1.5];
filler2 = fill(x,y,'r','FaceAlpha',.3);
x = [3650 3650 4000 4000];
y = [-1.5 3.5 3.5 -1.5];
filler3 = fill(x,y,'r','FaceAlpha',.3);
legend({'$\Delta$ Drivers', 'M','QE'},'Interpreter','latex','location','north')
set(gca,'TickLabelInterpreter','latex')
set(gca, 'LooseInset', get(gca,'TightInset'))
fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\Spectral Martingale Measures');
str=strcat('MonetaryBase');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

%%
TrailingPer = 20;
rho = zeros(134-TrailingPer,1);
for x = TrailingPer+1:134
    rho(x-TrailingPer) = corr(RCDeltaNorm(x-TrailingPer:x),MNorm(x-TrailingPer:x));
end
% figure
% hold on
% box on
% grid on
% plot(dates(TrailingPer+1:134),rho)
%% Save
save("MonetaryBase","XX")
