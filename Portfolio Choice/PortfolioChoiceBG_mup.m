%% Portfolio Choice

clear
clc
close all

%% Testing Parameters.
% MBG parameters
D = 2;

Ej = 'neg';

if strcmp(Ej,'pos')
    %Dec 31 2020
    bp = 7.5964e-03;
    cp = 1.5592e+00;
    bn = 1.8120e-02;
    cn = 6.3083e-01;
    
    mup = bp*cp;
    sigmap = bp*sqrt(cp);
    mun = bn*cn;    
    sigman = bn*sqrt(cn);

    N = 100;
    E = ones(N,1);
    Emat = [E,E];
    mup = [linspace(0.001,0.025,N)',mup*E];
    sigmap = [sigmap, sigmap].*Emat;
    mun = [mun,mun].*Emat;
    sigman = [sigman, sigman].*Emat;
    
else
    % Jan 2 2020
    bp = 6.779746e-03;
    cp = 2.877843e+01;
    bn = 5.210701e-02;
    cn = 6.110235e+00;
    
    mup = bp*cp;
    sigmap = bp*sqrt(cp);
    mun = bn*cn;    
    sigman = bn*sqrt(cn);

    N = 100;
    E = ones(N,1);
    Emat = [E,E];
    mup = [linspace(0.1,0.3,N)',mup*E];
    sigmap = [sigmap, sigmap].*Emat;
    mun = [mun,mun].*Emat;
    sigman = [sigman, sigman].*Emat;

end

vis = 1;

bpvec = zeros(N,2);
cpvec = zeros(N,2);
bnvec = zeros(N,2);
cnvec = zeros(N,2);

for i = 1 : N
    bpvec(i,:) = ( (sigmap(i,:).^2)./mup(i,:) );
    cpvec(i,:) = ( (mup(i,:).^2)./(sigmap(i,:).^2) );
    bnvec(i,:) = ( (sigman(i,:).^2)./mun(i,:) );
    cnvec(i,:) = ( (mun(i,:).^2)./(sigman(i,:).^2) );
end
    
% Distortion Parameters
c = 0.01;
gamma = 0.25;
aw = 100;
bw = 0.1;

%% Maximization

theta = zeros(N,2);
Sharpe = zeros(N,2);

thetai = 0.5;
for i = 1:N
    % MBG Parameters
    
    bp = bpvec(i,:)';
    cp = cpvec(i,:)';
    bn = bnvec(i,:)';
    cn = cnvec(i,:)';
    
    a = log( ((1./(1-bp)).^cp).*((1./(1+bn)).^cn) );
    
    % Maximization
    
    options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
    f = @(theta)-H([theta;1-theta],c,gamma,aw,bw,bp,cp,bn,cn,a);

     thetaplot = 0.01:0.01:0.99;
     ff = zeros(size(thetaplot));
     for ii = 1:length(thetaplot)
         ff(ii) = f(thetaplot(ii));
     end
%     figure
%     plot(thetaplot,-ff);
    
     [~,ii] = min(ff);
     thetai = thetaplot(ii);
    
    [thetai,Sopt,exitflagBG] = fminsearch(f,thetai,options); thetai = thetai-floor(thetai);
    theta(i,:) = [thetai,1-thetai];
    
    Sharpe(i,:) = ( (mup(i,:)-mun(i,:))./sqrt((sigmap(i,:).^2+sigman(i,:).^2)) );
    fprintf('i = %d, theta_i = (%d,%d), Sharpe = (%d,%d), bn = (%d,%d), cn = (%d,%d)\n',...
                i, theta(i,1), theta(i,2), Sharpe(i,1), Sharpe(i,2), bn(1), bn(2), cn(1), cn(2))
    %C = [exp(-Sopt), exp(f(1)), exp(f(0.01))]
end    

%% Visualization

if vis == 1    
    figure
    hold on
    box on
    grid on
    plot(mup(:,1),theta(:,1))
    ylabel('$\theta_*(\mu_{p1})$', 'Interpreter', 'latex')
    xlabel('$\mu_{p1}$', 'Interpreter', 'latex')
    set(gca,'TickLabelInterpreter','latex')
    fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
    str=strcat('OptimalControlBGmup_mup',Ej);
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

function d = dpsin(x,c,bw)
    d = bw*exp( -c * x );
end
 
function d = dpsip(x,c,gamma,aw)
    d = (aw*c/(1+gamma))...
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

function t = H(theta,c,gamma,aw,bw,bp,cp,bn,cn,a)
    theta = theta-floor(theta);
    t = zeros(2,1);
    for j = 1:2
        funp = @(yj) dpsin( nup(theta(j)*(exp(yj)-1),theta,bp,cp,bn,cn) , c, bw ).*...
                    kappa(yj,theta(j),bp(j),cp(j),bn(j),cn(j));
                
        funn = @(yj) dpsip( nun(theta(j)*(exp(yj)-1),theta,bp,cp,bn,cn) , c , gamma, aw).*...
                    kappa(yj,theta(j),bp(j),cp(j),bn(j),cn(j));
        
        Ip = integral(funp,0,Inf);
        In = integral(funn,-Inf,0);
        t(j) = Ip - In;
    end
    t = theta'*a-sum(t);
end