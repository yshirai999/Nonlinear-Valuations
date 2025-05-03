%% DataSet BG
close all
clear
clc

%% Parameters
q = 5;
L = 2^12;

rng('default'); rng(1);

%% Load Data
d = 6;
J = 4;
I = 46;
n = 3273;

T = string(zeros(I*J,1));
X = zeros(n,d);
M = zeros(n,J);

% Path to the Data folder
[parentFolder, ~, ~] = fileparts(pwd); % Get the parent folder of the current folder
dataFolder = fullfile(parentFolder, 'Data'); % Construct the path to the Data folder

for j=1:J
    str = fullfile(dataFolder, strcat('BGP', num2str(j), '.mat'));
    dat = load(str);
    str = strcat('d',num2str(j),'r');
    XX = dat.(str);
    for i = 1:I
        [n,~] = size(XX(i).parms);
        T(I*(j-1)+i) = string(XX(i).ticker);
        X(1:n,:,I*(j-1)+i) = XX(i).parms; 
        M(1:n,:,I*(j-1)+i) = [X(1:n,3,I*(j-1)+i).*X(1:n,2,I*(j-1)+i),... %c_pb_p
                                sqrt(X(1:n,3,I*(j-1)+i)).*(X(1:n,2,I*(j-1)+i)),... %sqrt(c_p)b_p
                                X(1:n,5,I*(j-1)+i).*X(1:n,4,I*(j-1)+i),... %c_nb_n
                                sqrt(X(1:n,5,I*(j-1)+i)).*(X(1:n,4,I*(j-1)+i))]; %sqrt(c_n)b_n
    end
end

%% Dataset

dates = X(:,1,1);

%save('dates','dates')
%save('T','T')

[n,d,I] = size(M);
M2D = transpose(reshape(permute(M,[2,1,3]),d,n*I));
MOut2D = rmoutliers(M2D,'percentiles',[0.5*q,100-0.5*q]);
%MOut2D = rmoutliers(M2D,'median'); %this eliminates too many

%% SPY data
SPY = fullfile(dataFolder, 'SPY.xlsm'); % Construct the full path to SPY.xlsm
SPY = xlsread(SPY);
SPY = SPY(:,2:end);
ind = (SPY(:,1)>20080102);
SPY = SPY(ind,:);

SPYBG = X(:,2:5,148);

Y = [SPY,zeros(size(SPY))];
for i = 1:n
    ii= find(SPY(:,1)==dates(i));
    Y(ii,5:8)=SPYBG(i,:);
end

save('Y','Y')

%% ETF data

ticker = {'XLB', 'XLE', 'XLF', 'XLI', 'XLK', 'XLP', 'XLU', 'XLV', 'XLY'};

for i=1:9
    dataFile = fullfile(dataFolder, strcat(ticker{i}, '_20082020')); % Construct the full path
    ETFi = load(dataFile); % Load the file from the Data folder
    ETFi = ETFi.P;
    ETFi = ETFi(:,2:end);
    ind = (ETFi(:,1)>20080102);
    ETFi = ETFi(ind,:);

    ETFiBG = X(:,2:5,173+i);

    Y = [ETFi,zeros(size(ETFi))];
    for j = 1:n
        ii= find(ETFi(:,1)==dates(j));
        Y(ii,5:8)=ETFiBG(j,:);
    end
    
    save(strcat('Y',ticker{i}),'Y')
end

