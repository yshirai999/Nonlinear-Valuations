function R = RfunDist ( i, dx, x, uj, bp, cp, bn, cn, w, gamma, aw, bw )
    f = @(y) zji(i,y,x,uj).*(1+psiU( i, y, x, uj, bp, cp, bn, cn, w, gamma, aw, bw))...
                .*kappa(y, bp,cp,bn,cn);
    Rp = integral(f,dx,Inf);
    Rm = integral(f,-Inf,-dx);
    R = Rp + Rm;
end


function z = zji ( i, y, x, uj )
    z = interp1(x,uj,x(i)+y,'spline','extrap') - uj(i);
end

function k = kappa ( y, bp, cp, bn ,cn )
    k = ( cp * exp ( - abs(y) / bp) ./ abs(y) ) .* ( y>0 )...
            + ( cn * exp ( - abs(y) / bn) ./ abs(y) ) .* ( y<0 );
end

function psi = psiU( i, y, x, uj, bp, cp, bn, cn, w, gamma, aw, bw)
    
    N = length(x);
    yk = x - x(i);
    
    s = zji ( i, y, x, uj );
    sk = zji ( i, yk, x, uj );
    
    [~,j0] = min ( sk );
    
    if j0 == N
        indp = N*ones(size(s));
    else
        [~,indp] = min( abs ( - sk(j0+1:end) ) ); %NB: s is a row vector and yk is column.        
        indp = indp + j0; %NB: size(indp) = size(s);
    end
    if j0 == 1
        indn = ones(size(s));
    else
        [~,indn] = min( abs( s - sk(1:j0-1) ) ); %NB: size(indn) = size(s);
    end
    
    xi1 = transpose(nu(yk(indp),bp,cp,bn,cn));
    xi2 = transpose(nu(yk(indn),bp,cp,bn,cn));
    xi1(~isfinite(xi1)) = 0;
    xi2(~isfinite(xi2)) = 0;
    psi = ( (dGammap(xi1,w,gamma,aw)) ).*( s>0 ) ...
            - ( (dGamman(xi2,w,bw)) ).*(s<0);
    psi(isnan(psi)) = 0;
%     figure
%     plot(y,psi);
end

function dG = dGammap (xi,w,gamma,aw)
    dG = (aw*w/(1+gamma))*((1-exp(-w*xi)).*(-gamma/(1+gamma))).*exp(-w*xi);
end

function dG = dGamman (xi,w,bw)
    dG = bw*exp(-w*xi);
end

function nueval = nu(x,bp,cp,bn,cn)
    nuevalp = ( cp*expint(x/bp) ).*(x>0); 
    nuevaln = ( cn*expint(abs(x)/bn) ).*(x<0);
    nuevalp(~isfinite(nuevalp)) = 0;
    nuevaln(~isfinite(nuevaln)) = 0;
    nueval = nuevalp+nuevaln;
end    

