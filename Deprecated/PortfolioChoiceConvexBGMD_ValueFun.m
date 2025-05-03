%% Portfolio Choice

clear
clc
close all

%% Testing Parameters.
% MBG parameters
D = 2;

BGParam = 'NP'; %NN, NP, PP

bp = [3.6476e-03, 4.3758e-03];
cp = [1.9298e+00, 2.5562e+00];
bn = [5.2109e-03, 7.1023e-03];
cn = [1.0983e+00, 1.3074e+00];

% bp = 7.5964e-03;
% cp = 1.5592e+00;
% bn = 1.8120e-02;
% cn = 6.3083e-01;

mup = bp.*cp;
sigmap = bp.*sqrt(cp);
mun = bn.*cn;
sigman = bn.*sqrt(cn);

%%
if strcmp(BGParam,'NN')
    sigmap(1) = [sigmap(1)/2];
elseif strcmp(BGParam,'NP')
    sigmap(1) = sigmap(1)*2;
elseif strcmp(BGParam,'PP')
    sigmap = [sigmap(1)*1.5,sigmap(2)*2];
end

bp = ( (sigmap.^2)./mup )';
cp = ( (mup.^2)./(sigmap.^2) )';
bn = ( (sigman.^2)./mun )';
cn = ( (mun.^2)./(sigman.^2) )';

vis_c = 0;
vis_p = 0;

% Distortion Parameters
cu = 0.5;
cl = 0.05;
gamma = 0.25;
alpha = 1;

%% Maximization

thetaplot = linspace(0.01,0.99,50);
cplot = linspace(0.1,0.4,5);

N = 10;
p_l = 1;
p_u = 1000;
p = linspace(p_l,p_u,N)';

theta = zeros(N,1);
b = zeros(N,1);

Fopt = zeros(N,1);
F = zeros(N,1);
c = zeros(length(thetaplot),N);
%thetai = 0.5;


try
    Fopt = load(strcat('Fopt_ConvexBG_',num2str(N),num2str(p_l*10),num2str(p_u)));
    Fopt = Fopt.Fopt;
    F = load(strcat('F_ConvexBG_',num2str(N),num2str(p_l*10),num2str(p_u)));
    F = F.F;
    c = load(strcat('c_ConvexBG_',num2str(N),num2str(p_l*10),num2str(p_u)));
    c = c.c;
catch
    for i = 1:N
        tic
        a = log( ( (1-bp).^(-cp) ) .* ( (1+bn).^(-cn) ) );
    
        % Maximization
    
        options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
        f = @(theta,c)H([theta;1-theta],c,gamma,bp,cp,bn,cn,a,cu,cl,alpha,p(i));
    
        ff = zeros(size(thetaplot));
        fiiplot = zeros(size(cplot));
    
        for ii = 1:length(thetaplot)
            fii = @(c) f(thetaplot(ii),c);
            for cind = 1:length(cplot)
                fiiplot(cind) = fii(cplot(cind));
            end
            [~,cind] = min(fiiplot);
            c(ii,i) = cplot(cind);
            [c(ii,i),ff(ii)] = fminsearch(fii,c(ii,i),options);
        end
        [~,ii] = max(ff);
        thetai = thetaplot(ii);
    
        fi = @(theta)-f(theta,c(ii,i));
        [thetai,Fopt(i),exitflagBG] = fminsearch(fi,thetai,options); 
        %thetai = thetai-sign(thetai)*floor(abs(thetai));
        thetai = thetai-floor(thetai); %Long ONLY
        theta(i) = thetai;
    
        cputime = toc;
        
        F(i) = f(c(25,i),thetaplot(25));
    
        b(i) = beta(c(ii,i),cu,cl,alpha);
    
        fprintf('i = %d, Fopt = %d, beta = %d, theta = %d, Expected Time to Complete < %d\n', i, Fopt(i), b(i), thetai, (N-i)*cputime)
    end
end

plot(p,F)

%% 

save(strcat('Fopt_ConvexBG_',num2str(N),num2str(p_l*10),num2str(p_u)),'Fopt')
save(strcat('F_ConvexBG_',num2str(N),num2str(p_l*10),num2str(p_u)),'F')
save(strcat('c_ConvexBG_',num2str(N),num2str(p_l*10),num2str(p_u)),'c')

T = 100; %time horizon in days
dt = 1;
tt = dt:dt:T;

NN = 10;
p_l = 1;
p_u = 19;
pp = linspace(p_l,p_u,NN)';

PIopt = zeros(NN,length(tt));
PIopt(:,length(tt)) = pp;
PI = zeros(NN,length(tt));
PI(:,length(tt)) = pp;
for t = length(tt)-1:-1:1
    for i = 1:NN
        [~,ind] = min( (PIopt(i,t+1)-(p+dt*Fopt)).^2 );
        PIopt(i,t) = p(ind);
        [~,ind] = min( (PI(i,t+1)-(p+dt*F)).^2 );
        PI(i,t) = p(ind);
    end
end

%plot(pp,PI(:,1));

%% Visualization

if vis_p == 1
    figure
    hold on
    box on
    grid on
    plot(pp,PIopt(:,1))
    %ylabel('$\theta_*(\varpi)$', 'Interpreter', 'latex')
    title(strcat('Value Function for $c_u = ',num2str(cu),'$'),'interpreter','latex')
    xlabel('$\varpi$', 'Interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
    fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
    str=strcat('OptimalControlConvexBG_ValueFun');
    fname=str;
    saveas(gcf, fullfile(fpath, fname), 'epsc');
    hold off
end 

%% Routines

function n = nup(p,theta,bp,cp,bn,cn)
    n = zeros(size(p));
    ind = (p>0);
    for k=1:length(theta)
        if theta(k)>0
            n(ind) = n(ind)+cp(k)*expint(log(p(ind)/theta(k)+1)/bp(k));
            n(isnan(n)) = 0;
        elseif theta(k)<0
            indn = (theta(k)<-p);
            ind = logical(indn.*ind);
            n(ind) = n(ind)+cn(k)*expint(-log(p(ind)/theta(k)+1)/bn(k));
        end
    end
end

function n = nun(p,theta,bp,cp,bn,cn)
    n = zeros(size(p));
    ind = (p<0);
    for k=1:length(theta)
        if theta(k)<0
            n(ind) = n(ind)+cp(k)*expint(log(p(ind)/theta(k)+1)/bp(k));
            n(isnan(n)) = 0;
        elseif theta(k)>0
            indn = (theta(k)>-p);
            indk = logical(indn.*ind);
            n(indk) = n(indk) + cn(k)*expint(-log(p(indk)/theta(k)+1)/bn(k));
        end
    end
end

function d = dpsin(x,c)
    d = 0.1*exp( -c * x );
end
 
function d = dpsip(x,c,gamma)
    d = (1/(1+gamma))...
        .* ( ( 1-exp( - c * x ) ) .^ (-gamma/(1+gamma)) )...
        .* exp( -c * x );
    d(~isfinite(d)) = 0;
end

function k = kappa(x,theta,bp,cp,bn,cn)
    k = theta.*(exp(x)-1).*...
        ( cp * ( exp(-abs(x)/bp) ./ abs(x) ) .* (x>0)...
        + cn * ( exp(-abs(x)/bn) ./ abs(x) ) .* (x<0));
    k(isnan(k))=0;
end

function t = H(theta,c,gamma,bp,cp,bn,cn,a,cu,cl,alpha,P)
    %theta = theta-sign(theta).*floor(abs(theta));
    theta = theta-floor(theta);
    t = zeros(2,1);
    for j = 1:2
        funp = @(yj) dpsin( nup(theta(j)*(exp(yj)-1),theta,bp,cp,bn,cn) , c ).*...
                    kappa(yj,theta(j),bp(j),cp(j),bn(j),cn(j));
                
        funn = @(yj) dpsip( nun(theta(j)*(exp(yj)-1),theta,bp,cp,bn,cn) , c , gamma).*...
                    kappa(yj,theta(j),bp(j),cp(j),bn(j),cn(j));
        
        Ip = integral(funp,0,Inf);
        In = integral(funn,-Inf,0);
        t(j) = Ip - In;
    end
    
    t = min( P*(theta'*a-sum(t))+beta(c,cu,cl,alpha),P*(theta'*a) );
    %t = P*(theta'*a-sum(t))+beta(c,cu,cl,alpha);
end

function b = beta(c,cu,cl,alpha)
    b = alpha*exp( 1./(c-cl) ).*exp(-1./(cu-c)).*(c<cu);
    b(~isfinite(b)) = 0;
    b(isnan(b)) = 0;
end


