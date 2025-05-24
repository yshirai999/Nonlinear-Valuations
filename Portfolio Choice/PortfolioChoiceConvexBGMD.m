%% Portfolio Choice

clear
clc
close all

% Testing Parameters.
% Rate and Term
r = 0.0001;
dt = 252;

% BG parameters
LL = 1;
bp = 6.779746e-03;
cp = 2.877843e+01*dt;
bn = 5.210701e-02;
cn = 6.110235e+00*dt;

bn = 6.779746e-03; %short position
cn = 2.877843e+01*dt;
bp = 5.210701e-02;
cp = 6.110235e+00*dt;

bp = 0.0075;
cp = 1.5592*dt;
bn = 0.0181;
cn = 0.6308*dt;

% Rebate and Distortion Parameters
cu = 2;
cl = 0.25;
gamma = 0.01;

chi = [1,5,10];
chi2 = 2*[1,1,1];
M = length(chi);

figure
hold on
cplot = linspace(cl*1.001,cu*0.999,50);

for j=1:M
    alpha = [chi(j),chi2(j)];
    b = beta(cplot,cu,cl,alpha);
    plot(cplot,b)
    %leg{j} = strcat('$\chi = $', num2str(chi(j)));
end
%legend = legend(leg,'interpreter','latex','location','best');
xlabel('$\theta$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')

%% Maximization

vizPath = getPath('Visualization');

N = 100;
p_l = 10;
p_u = 1000;
p = linspace(p_l,p_u,N)';

b = zeros(N,M);
c = zeros(N,M);
L = zeros(N,M);

strplot = {'-.','--','-'};

figure
hold on
for j=1:M
    alpha = [chi(j),chi2(j)];
    for i = 1:N
        fprintf('j = %d,i=%d\n',j,i)
        a = log( ( (1-bp).^(-cp) ) .* ( (1+bn).^(-cn) ) );
    
        % Maximization
    
        options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
        f = @(c)H(c,gamma,bp,cp,bn,cn,cu,cl,alpha,p(i));
    
        fcplot = zeros(size(cplot));
    
        for cind = 1:length(cplot)
            fcplot(cind) = f(cplot(cind));
        end
    
        [~,cind] = min(fcplot);
        c(i,j) = cplot(cind);
    
        %[c(i,j),L(i,j)] = fminsearch(f,c(i,j),options);
        L(i,j) = f(c(i,j));
        b(i,j) = beta(c(i,j),cu,cl,alpha);
        L(i,j) = L(i,j)/LL;
    end
    
    plot(p/p(end),(p(end)-p)*dt*r+p*a+L(:,j),strplot{j})
    [~,theta] = max((p(end)-p)*dt*r+p*a+L(:,j));
    theta = theta/N;
    leg{j} = strcat('$\chi = $', num2str(chi(j)));
end

%plot(p,p*a)
%leg{j+1} = strcat('$\varpi a $');
legend = legend(leg,'interpreter','latex','location','best');
xlabel('$\theta$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')
str=strcat('ConvexPortfolioChoiceMD');
fname=str;
saveas(gcf, fullfile(vizPath, fname), 'epsc');
hold off

% for i = 1:N-1
%     fprintf('(%d,%d) ',p(i),p(i+1))
%     for j=1:M
%         fprintf('& %d ',p(i+1)*a+L(i+1,j)-p(i)*a-L(i,j) )
%     end
%     fprintf('\n')
% end

%[p*a,L]

%% Routines

function d = dpsin(x,c)
    d = 0.1*exp( -c * x );
end
 
function d = dpsip(x,c,gamma)
    d = (1/(1+gamma))...
        .* ( ( 1-exp( - c * x ) ) .^ (-gamma/(1+gamma)) )...
        .* exp( -c * x );
    d(~isfinite(d)) = 0;
end

function k = kappa(x,bp,cp,bn,cn)
    k = (exp(x)-1).*...
        ( cp * ( exp(-abs(x)/bp) ./ abs(x) ) .* (x>0)...
        + cn * ( exp(-abs(x)/bn) ./ abs(x) ) .* (x<0));
    k(isnan(k))=0;
end

function t = H(c,gamma,bp,cp,bn,cn,cu,cl,alpha,P)
    funp = @(y) dpsin( cp*expint(y/bp), c ).*...
                    kappa(y,bp,cp,bn,cn);
                
    funn = @(y) dpsip( cn*expint(y/bn), c, gamma).*...
                    kappa(-y,bp,cp,bn,cn);
        
    Ip = integral(funp,0,Inf);
    In = integral(funn,0,Inf);
    t = Ip - In;
    t = -P*t+beta(c,cu,cl,alpha);
end

function b = beta(c,cu,cl,alpha)
    b = alpha(1)*exp( 1./(c-cl).^alpha(2) ).*exp(-1./(cu-c).^alpha(2)).*(c<cu);
    
    b = ( exp(alpha(1)*(c-cu))-1-alpha(1)*(c-cu) ) ./ (c-cl).^alpha(2);
    
    b(~isfinite(b)) = 0;
    b(isnan(b)) = 0;
end




