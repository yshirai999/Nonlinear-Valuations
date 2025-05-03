%% Portfolio Choice

clear
clc
close all

%% Testing Parameters.
% MBG parameters
D = 2;

N = 10;
E = ones(N,1);
Emat = [E,E];
mup = [4295.2238, 4295.2238].*Emat;
sigmap = [800.6735, 800.6735].*Emat;
mun = [117.2783, 117.2783].*Emat;
sigman = [[linspace(47.4449/10,47.4449,N/2)';linspace(47.4449,10*47.4449,N/2)'], 47.4449*E];

bpvec = zeros(N,2);
cpvec = zeros(N,2);

bnvec = zeros(N,2);
cnvec = zeros(N,2);

for i = 1 : N
    bpvec(i,:) = ( mup(i,:)./(sigmap(i,:).^2) )';
    cpvec(i,:) = ( (mup(i,:).^2)./(sigmap(i,:).^2) )';
    bnvec(i,:) = ( mun(i,:)./(sigman(i,:).^2) )';
    cnvec(i,:) = ( (mun(i,:).^2)./(sigman(i,:).^2) )';
end

Cvec = [1, 0, 0, 1].*[Emat,Emat];
    
% Distortion Parameters
c = 0.5;
gamma = 0.25;
 
% QuasiRandom Sequence
Nsim = 300;
rng default  % For reproducibility
p = haltonset(4,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
U = net(p,Nsim);
Up = U(:,[2,1]); Um = U(:,3:4);

% Y = -log(Up')+log(Um');
% Upsilon = [-log(Up(:,1)');-log(Up(:,2)')+log(Um(:,2)')];

Y = 2*Up'-1;
Upsilon = 2*Up'-1;

e = ones(D,1);

%% Maximization

theta = zeros(N,2);
Sharpe = zeros(N,2);

for i = 1:N
    % MBG Parameters
    
    bp = bpvec(i,:)';
    cp = cpvec(i,:)';
    bn = bnvec(i,:)';
    cn = cnvec(i,:)';
    v = 1/min([cp;cn])+0.01;
    C = [Cvec(i,1:2);
         Cvec(i,3:4)];
     
    b = (bp-bn)/v;
    sigma = sqrt(2*bp.*bn/v);
    Sigma = diag(sigma)*C*diag(sigma);
    Sigmainv = inv(Sigma);
    detSigma = det(Sigma);
    
    a = ((1./(1-b+Sigma*e)).^(1/v)).*((1./(1-bp)).^(cp-1/v)).*((1./(1+bn)).^(cn-1/v));
    
    % Maximization
    
    options = optimset('MaxFunEval',5000,'MaxIter',100,'Display','off','TolFun',1e-9,'TolX',1e-9);
    f = @(theta)-H([abs(theta);1-abs(theta)],Y,c,gamma,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn,a);
    
    [thetai,Sopt,exitflagBG] = fminsearch(f,0.5,options); thetai = thetai-floor(thetai); theta(i,:) = [abs(thetai),1-thetai];
    
    % theta = [0.45;0.55];
    % y = [0.9070;
    %    -1.1208];
    % p = theta'*exp(y)-1;
    % test = Int(theta,p,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn);
    
    % x = linspace(0.001,0.005,100);
    % y = linspace(0.001,0.005,100);
    % x = x(x~=0);
    % y = y(y~=0);
    % [X,Y] = meshgrid(x,y);
    % Z = X+Y;
    % for i=1:length(x)
    %     for j = 1:length(y)
    %         Z(i,j) = VGkappa([x(i);y(j)],D,b,v,Sigmainv,detSigma);
    %     end
    % end
    % contour(X,Y,Z)
    
    Sharpe(i,:) = ( (mup(i,:)-mun(i,:))./sqrt((sigmap(i,:).^2+sigman(i,:).^2)) );
    fprintf('i = %d, theta_i = (%d,%d), Sharpe = (%d,%d), bn = (%d,%d), cn = (%d,%d)\n',...
                i, theta(i,1), theta(i,2), Sharpe(i,1), Sharpe(i,2), bn(1), bn(2), cn(1), cn(2))
    %C = [exp(-Sopt), exp(f(1)), exp(f(0.01))]
end    

%% Visualization
figure
hold on
grid on
plot(bnvec(:,1),theta(:,1))
ylabel('$\theta_*(\sigma_n)$','Interpreter','latex')
xlabel('$$b_n$$', 'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('OptimalControlsigman_bn');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

figure
hold on
grid on
plot(Sharpe(:,1),theta(:,1))
ylabel('$\theta_*(\sigma_n)$','Interpreter','latex')
xlabel('$\frac{\mu}{\sigma}$', 'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('OptimalControlsigman_Sharpe');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

figure
hold on
grid on
plot(mun(:,1),theta(:,1))
ylabel('$\theta_*(\sigma_n)$','Interpreter','latex')
xlabel('$\mu_n$', 'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('OptimalControlsigman_mun');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off

figure
hold on
grid on
plot(sigman(:,1),theta(:,1))
ylabel('$\theta_*(\sigma_n)$','Interpreter','latex')
xlabel('$$\sigma_n$$', 'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
fpath=('C:\Users\Yoshihiro Shirai\Desktop\PhD\Research\CDXO Nonlinear Valuation');
str=strcat('OptimalControlsigman_sigman');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');
hold off


%% Routines

function mtilde = VGkappa(y,D,b,v,Sigmainv,detSigma) %y is D x N
    mtilde = ( exp(b'*Sigmainv*y)' ./ ( v*(2*pi)^(D/2-1)*sqrt(abs(detSigma))*sqrt(diag((y'*Sigmainv)*y)) ) ).*...
            exp(-sqrt( (b'*Sigmainv*b+2/v)*(diag(y'*Sigmainv*y))));
end

function kappa = BGkappa(y,D,bp,cp,bn,cn,v)
    kappa = zeros(D,1);
    w = zeros(D,1);
    ind = 1:1:D;
    for j=1:D
        indj = ind(ind~=j);
        w(j) = prod(y(indj)==0);
        kappa(j) = ((cn(j)-1/v)/abs(y(j)))*exp(-abs(y(j))/bn(j))*(y(j)<0)+...
                 ((cp(j)-1/v)/y(j))*exp(-y(j)/bp(j))*(y(j)<0)*w(j);
    end
    kappa = sum(kappa);
end

function kappa = MBGkappa(y,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)
    kappa = VGkappa(y,D,b,v,Sigmainv,detSigma)+BGkappa(y,D,bp,cp,bn,cn,v);
end

function upsilon = F(theta,y)
    upsilon = zeros(size(y));
    upsilon(1) = theta'*exp(y)-1;
    for j = 2:length(y)
        upsilon(j) = y(j)-y(j-1);
    end
end

function y = Finv(theta,upsilon)
%     [D,~] = size(upsilon);
%     j = 1;
%     while theta(j)==0
%         j=j+1;
%     end
%     y(j,:) = log((upsilon(1,:)+1)/theta(j));
%     jj = j;
%     while jj>1
%         y(jj-1,:)=y(jj,:)-upsilon(jj-1,:);
%         jj = jj-1;
%     end
%     jj = j;
%     while jj < D
%         y(jj+1,:)=upsilon(jj+1,:)+y(jj,:);
%         jj = jj+1;
%     end
    y = zeros(size(upsilon));
    if theta(1)>0    
        y(1,:) = log( (upsilon(1,:)+1)/(theta(1)+theta(2)*exp(upsilon(2,:))) );
        y(2,:) = upsilon(2,:)+y(1,:);
    elseif theta(1)==0
        y(2,:) = log( (upsilon(1,:)+1)/(theta(2)+theta(1)*exp(upsilon(1,:))) );
        y(1,:) = y(2,:)-upsilon(2,:);
    end
end

function I = Int(theta,p,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)
    d = length(p);
    I = p;
    for j = 1:d
        [~,Nsim] = size(Upsilon);
        pvec = p(j)*ones(1,Nsim);
        y = Finv(theta,[pvec;Upsilon(2:end,:)]);
%         I(j) = (p(j)\abs(p(j)))*(1/(p(j)+1))*...
%             sum(MBGkappa(y,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn).*...
%             exp(sum(abs(Upsilon(2:end,:))))/2)/Nsim;
        I(j) = -sign(p(j))*(p(j)+1)*sum(MBGkappa(y,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn))/Nsim;
    end
end


function psi1 = psihat(theta,p,c,gamma,Y,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)
    [~,Nsim] = size(Y);
    I = @(s)Int(theta,s,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn);
    k = I(p);
    if p > 0
%         nu = sum( MBGkappa(Y,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)'.*(theta'*exp(Y)>p+1).*...
%                     exp(sum(abs(Y)))/2)/Nsim;
        nu = sum( MBGkappa(Y,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)'.*(theta'*exp(Y)>p+1) )/Nsim;
        psi1 = -exp(-c*nu)*k;
    elseif p < 0
%         nu = max(sum( MBGkappa(Y,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)'.*(theta'*exp(Y)<p+1).*...
%                     exp(sum(abs(Y)))/2)/Nsim, 1e-15);
        nu = max(sum( MBGkappa(Y,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)'.*(theta'*exp(Y)<p+1) ));
        psi1 = ( (1/(1+gamma))*(1-exp(-c*nu)).^(-gamma/(1+gamma)) )*...
                    exp(-c*nu)*k;
    end
end

function psi2 = psi(theta,y,c,gamma,Y,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)
    psi2 = (theta'*exp(y)) *...
            psihat(theta,theta'*exp(y)-1,c,gamma,Y,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn) /...
            Int(theta,theta'*exp(y)-1,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn);
end

function t = H(theta,Y,c,gamma,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn,a)
    theta = theta-floor(theta);
    [~,Nsim] = size(Y);
    t = zeros(Nsim,1);
    for j = 1:Nsim
        t(j) = psi(theta,Y(:,j),c,gamma,Y,Upsilon,D,b,v,Sigmainv,detSigma,bp,cp,bn,cn)/Nsim;
    end
    %[t,[1:1:Nsim]']
    t = t(~isnan(t));
    t = t(~isinf(t));
    t = theta'*a-sum(t(~isnan(t)));
end