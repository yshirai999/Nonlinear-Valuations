function R = Rfun ( i, dx, x, uj, bp, cp, bn, cn )
    f = @(y) zji ( i, y, x, uj ) .* kappa ( y, bp, cp, bn, cn );
    Rp = integral(f,dx,Inf);
    Rm = integral(f,-Inf,-dx);
    R = Rp + Rm;
end

function z = zji ( i, y, x, uj )
%     K = length(y);
%     for k = 1 : K
%         z = sum( (uj-uj(i)).*([x(2:end);Inf]>x(i)+y(k)).*([-Inf;x(1:end-1)]<x(i)+y(k)) );
%     end
    
%     y = min(max(x)-x(i),y);
%     y = max(min(x)-x(i),y);
    z = interp1(x,uj,x(i)+y,'spline','extrap') - uj(i);
    %z(isnan(z)) = 0;
    %plot(y,z)
end

function k = kappa ( y, bp, cp, bn ,cn )
    k = ( cp * exp ( - abs(y) / bp) ./ abs(y) ) .* ( y>0 )...
            + ( cn * exp ( - abs(y) / bn) ./ abs(y) ) .* ( y<0 );
end
