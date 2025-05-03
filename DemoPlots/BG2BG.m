%% BG2BG Distortion
clear
clc
close all
% Creates a plot of the distortion function for the BG2BG model
% with the parameters from the calibration to Option data

%% Parameters
param0 = [3.876142e-03,6.145676e+02,9.791253e-02,3.717541e+00,3.924071e-03,9.727626e-02,3.827786e-03,9.854453e-02];
bp = param0(1);
cp = param0(2);
bn = param0(3);
cn = param0(4);
bpU = param0(5);
bnU= param0(6);
bpL = param0(7);
bnL = param0(8);

x = 0.1:0.1:30;
Up = zeros(size(x)); Un = zeros(size(x));
Ep = 1;
En = 1;
for i = 1:length(x)
    Ep = expintinv(x(i)*cp,Ep);
    En = expintinv(x(i)*cn,En);
    Up(i) = (1/cp)*expint( Ep*bp/bpU ) - x(i);
    Un(i) = -(1/cn)*expint( En*bn/bnU ) + x(i);
    fprintf('step i=%d of %d: U = %d\n', i,length(x),Up(i))
end

%% Visualization

% Dynamically construct the path to the visualization folder
[parentFolder, ~, ~] = fileparts(pwd); % Get the parent folder of the current folder
visFolder = fullfile(parentFolder, 'Visualization'); % Path to the Visualization folder

% Ensure the Visualization folder exists
if ~exist(visFolder, 'dir')
    mkdir(visFolder);
end

figure
hold on
grid on
box on
plot(x,Up,'-');
% plot(x,Un,'--');
% plot(x,-log(bp/bpU)*ones(size(x))/cp)
% plot(x,-log(bn/bnU)*ones(size(x))/cn)
%legend('$\Upsilon^+$','$\Upsilon^-$','interpreter','latex')
%set(gca,'TickLabelInterpreter','latex','yscale','log')
fname = 'BG2BG';
saveas(gcf, fullfile(visFolder, fname), 'epsc');

hold off

% x = 0.1:0.1:10;
% Up = -(1/cpU)*expint( expintinv(x*cp)*(bp/bpU) ) + x;
% figure
% hold on
% plot(x,Up);
% hold off

%% Routines

function x = expintinv(y,x0)
options = optimset('MaxFunEval',5000,'MaxIter',5000,'Display','off','TolFun',1e-9,'TolX',1e-9);
x = real(fminsearch( @(x)(abs(real(expint(x)-y))), x0, options )); 
end