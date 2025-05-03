%% Straddle
clear
clc
close all

%% Parameters
% Market parameters
T = 1; % option maturity in months
S0 = 1;
K = S0; % strike price

% BG parameters
bp = 0.0067; %risk neutral 2 January 2020
cp = 28.7780;
bn = 0.0521;
cn = 6.1102;

% bp = 0.0075; %physical 31 December 2020
% cp = 1.5592;
% bn = 0.0181;
% cn = 0.6308;

omega = ((1-bp)^(-T*cp))*((1+bn)^(-T*cn));

% Distortion parameter
w = 0.01;
aw = 100;
bw = 0.1;
gamma = 0.25;

% Numerical scheme parameters
N = 400; % space discretization
M = 50; % time discretization
Smax = 5*S0;
Smin = 0.01;
Xmax = log ( Smax );
Xmin = log ( Smin );
dx = ( Xmax - Xmin ) / N;
dt = T / M;

% g_p = expint ( [ 1 : 1 : N ] * dx / b_p );
% g_n = expint ( [ 1 : 1 : N ] * dx / b_n );
% gg_n = expint ( [ 1 : 1 : N ] * dx * ( 1 / b_n + 1 ) );
% sigma2 = c_p * b_p ^ 2 * ( - dx * exp ( - dx / b_p ) / b_p + 1 -...
%     exp ( - dx / b_p ) ) +  c_n * b_n ^ 2 *...
%     ( - dx * exp ( - dx / b_n ) / b_n + 1 - exp ( - dx / b_n ) );
% omega = c_n * g_n ( 1 ) - c_n * expint ( dx * ( 1 / b_n + 1 ) )...
%     + c_p * g_p ( 1 ) - c_p * expint ( dx * ( 1 / b_p - 1 ) );

% MC Simulation
Nsim = 100;
U = rand(1,Nsim);
Yp = - bp * log ( U );
Yn = - bn * log ( U );
MC = 0;

% Initial condition
X = ( Xmin : dx : Xmax )';
u = zeros ( M + 1, N + 1 );
u ( 1, 1 : end ) = max ( exp( X ) - K, 0 ) + max ( K - exp ( X ), 0 );
uj = u ( 1, 1 : end )';

%% NonLinear PIDE

% Stiffness matrix
a = - ones(1,N+1) / dt;
b = ones(1,N+1) * ( 1 + omega * dt / dx );
c = - ones(1,N+1) / dt;

% Neumann boundary conditions
b(2) = b(2) + 2 * a(2) / (1+0.5*dx);
b(N) = b(N) + 2 * c(N) / (1-0.5*dx);
c(2) = c(2) - a(2) * (1-0.5*dx) / (1+0.5*dx);
a(N) = a(N) - c(N) * (1+0.5*dx) / (1-0.5*dx);

R = zeros(1,N+1);
RC = zeros(1,N+1);
for j = 1:M
    for i = 2:N %317:317
        if MC == 0
            R ( i ) = RfunDist ( i, dx, X, uj, bp, cp, bn, cn, w, gamma, aw, bw );
        else
            R ( i ) = RfunDistMC ( i, X, uj, bp, cp, bn, cn, w, gamma, aw, bw, Yp, Yn, Nsim );        
        end
        if j == M
        RC ( i ) = R ( i ) - Rfun ( i, dx, X, uj, bp, cp, bn, cn );
        end
    end
        uj = uj + dt * R';
        TriDiagSolver ( a, b, c, uj, N ); % solve A*U_{j+1}=U_j+dx*R, U_{j}=u_{j,1},...,u_{j,N};
        uj(1) = 2 * uj ( 2 ) / ( 1 + 0.5 * dx )...
            - ( 1 - 0.5 * dx ) / ( 1 + 0.5 * dx ) * uj ( 3 );
        uj ( N + 1 ) = - ( 1 - 0.5 * dx ) * uj ( N - 1 ) / ( 1 + dx / 2 )...
             + 2 * uj ( N ) / ( 1 + 0.5 * dx );
        u(j+1,1:end)=uj';
        j
end

%% Visualization
% Comparison
S = exp(X);

figure
hold on
box on
grid on
plot ( S, abs(uj) );
% plot ( S, abs(u ( end - 12, : )) );
% plot ( S, abs(u ( 25, : )) );
% plot ( S, abs(u ( 12, : )) );
% plot ( S, abs(u ( 1, : )) );
%legend('$t = 0$', '$t = 0.25T$', '$t = 0.5T$', '$t = 0.75T$', '$t = T$', 'interpreter', 'latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
fname = strcat('Straddle',num2str(N),num2str(M));
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

figure
box on; clf; hold on; grid;
ma = max ( max ( abs( u ) ) );
mi = min ( min ( abs( u ) ) );
contourf ( S, [0 : dt : T]', abs( u ), linspace ( mi, ma, 20 ) );
%colormap(gray)
colorbar
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
fname = strcat('StraddleContourf',num2str(N),num2str(M));
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

% Risk Charge
figure
hold on
box on
grid on
plot(S,RC)
hold off
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
fname = strcat('Driver_M',num2str(N),num2str(M));
saveas(gcf, fullfile(fpath, fname), 'epsc');

%close all
