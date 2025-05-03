%% Digital Moments
clear
clc
close all

%% Load Data
ticker = {'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY', 'SPY'};


SY = 09; startdate = strcat('20',num2str(SY,'%02.f'),'0106');
SD = str2double(startdate);

N = 250;

%% Estimation
Delta = 5;

for jj = 1:length(ticker)
    Y = load(strcat('Y',ticker{jj}));
    Y = Y.Y;
    imin = find(Y(:,1)==SD);
    imax = length(Y);
    d = Y(imin:Delta:imax,1);
    n = length(imin:Delta:imax);
    mu_mkt = zeros(n,3);
    
    for i = imin+Delta:Delta:imax
        ind = [Delta:Delta:N];
        NN = length(ind);
        YY = zeros(NN+1,3);
        YY(1,1) = Y(i-N,1);
        YY(1,2) = min(Y(i-N-Delta+1:i-N,2));
        YY(1,3) = max(Y(i-N-Delta+1:i-N,3));
        YY(1,4) = mean([YY(1,2),YY(1,3)]);
        for j=2:NN+1
            YY(j,1) = Y(i-N+ind(j-1),1);
            YY(j,2) = min(Y(i-N+ind(j-1)-Delta+1:i-N+ind(j-1),2));
            YY(j,3) = max(Y(i-N+ind(j-1)-Delta+1:i-N+ind(j-1),3));
            YY(j,4) = mean([YY(j,2),YY(j,3)]);
        end
        u = ( log(YY(2:NN+1,3)') - log(YY(1:NN,3)') )*252/Delta;
        l = ( log(YY(2:NN+1,2)') - log(YY(1:NN,2)') )*252/Delta;
        %m = ( log(YY(2:NN+1,4)') - log(YY(1:NN,4)') )*252/Delta;
        m = ( log(Y(i-N+1:i,4)') - log(Y(i-N:i-1,4)') )*252;
        
        %fprintf('NN = %d, u = %d', NN, length(u))

        mu_mkt((i-imin)/Delta+1,:) = [mean(u),mean(m),mean(l)];
        %mu_mkt((i-imin)/Delta+1,:) = [ExpEmpLoss(u,NN),ExpEmpLoss(m,N),ExpEmpLoss(l,NN)];
    end
    %fprintf(strcat(ticker{jj},' & %4.2f & %4.2f \\\\ \n'),mean(mu_mkt(:,[1,3]))*100);
    fprintf(strcat(ticker{jj},' & %4.4f & %4.4f & %4.4f \\\\ \n'),mean(mu_mkt));
end

%% SPY
Y = load(strcat('Y','SPY'));
Y = Y.Y;
d = datetime(Y(:,1),'ConvertFrom','yyyymmdd'); 
if SY == 8 
    d = d(2:end);
end
figure
hold on
box on
grid on
plot(d(505:end),Y(505:end,2))
title('SPY','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\Spectral Martingale Measures');
str=strcat('SPY');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

%% Routines
function mu = ExpEmpLoss(r,N)
%fprintf('%d,%d',N,length(r))
s = r;% + RCU(theta);
s = sort(s); pi_N = linspace(1,N,N)/N;
Pi = linspace(1,999,999)/1000;
s = interp1(pi_N,s,Pi,'linear','extrap');

Pi(s>0) = 1-Pi(s>0);

Pim = Pi(s<0);
Pip = Pi(s>0);
sm = s(s<0);
sp = s(s>0);

mu = -Pim(2:end)*(sm(2:end)-sm(1:end-1))'+Pip(2:end)*(sp(2:end)-sp(1:end-1))';

end