%% Portfolio Choice

clear
clc
close all

%% Testing Parameters.
% BG parameters
BGParam = 'P';
dt = 1;
LL = 1;
bp = 7.5964e-03;
cp = 1.5592e+00*dt;
bn = 1.8120e-02;
cn = 6.3083e-01*dt;

% Distortion Parameters

if BGParam == 'N'
    cu = 1;
    cl = 0.05;
    gamma = 0.01;
    chi = [10,50,250,1250];%[10,50,100,500];
    chi2 = 2*[1,1,1,1];
elseif BGParam == 'P'
    cu = 1000;
    cl = 200;
    gamma = 0.01;
    chi = [10,50,250,1250];%[10,50,100,500];
    chi2 = 2*[1,1,1,1];
end
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

%% Maximization

cplot = linspace(cl*1.001,cu*0.999,25);
% b = beta(cplot,cu,cl,alpha);
% plot(cplot,b)

N = 100;
if BGParam == 'N'
    p_l = 1000;
    p_u = 100000;
elseif BGParam == 'P'
    p_l = 1000;
    p_u = 100000;
end
p = linspace(p_l,p_u,N)';

b = zeros(N,M);
c = zeros(N,M);
L = zeros(N,M);

strplot = {':','-.','--','-'};
figure
hold on
for j=1:M
    alpha = [chi(j) ,chi2(j)];
    for i = 1:N
        fprintf('j = %d, i = %d\n',j,i)
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
        %c(i,j)=200;
        L(i,j) = f(c(i,j));
        b(i,j) = beta(c(i,j),cu,cl,alpha);
        L(i,j) = L(i,j)/LL;
    end
    
    plot( p/p(1),(p*a+L(:,j))/p(1),strplot{j} )
    leg{j} = strcat('$\chi = $', num2str(chi(j)));
end

%plot(p,p*a)
%leg{j+1} = strcat('$\varpi a $');
legend = legend(leg,'interpreter','latex','location','best');
xlabel('$\varpi$', 'Interpreter', 'latex')
grid on
box on
set(gca,'TickLabelInterpreter','latex')
visPath = getPath('Visualization'); % Get the path to the Visualization folder
str=strcat('DRoS',BGParam);
fname=str;
saveas(gcf, fullfile(visPath, fname), 'pdf');
hold off

for i = 1:N-1
    fprintf('(%d,%d) ',p(i),p(i+1))
    for j=1:M
        fprintf('& %d ',p(i+1)*a+L(i+1,j)-p(i)*a-L(i,j) )
    end
    fprintf('\n')
end

%[p*a,L]

%% Routines

function d = dpsin(x,c)
    d = exp( -c * x );
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
    
    b = ( exp(alpha(1)*(c-cu))-1-alpha(1)*(c-cu) ) ./...
            (c-cl/2).^alpha(2);
    
    b(~isfinite(b)) = 0;
    b(isnan(b)) = 0;
end


