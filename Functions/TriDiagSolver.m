function TriDiagSolver( a, b, c, u, N )
% Solves Ax = u when A is tridiagonal. The solution x is stored in the
% initial vector u. The function leaves the vectors a, b and c unchanged.
% The vectors a, b c and u are in R^{N+1}. A is an N-1 x N-1 matrix. The
% vector x is in R^{N-1}, and its entries are saved in u(2),...,U(N).

aa = a;
bb = b;
cc = c;

for i = N:-1:3
    xMult = aa(i-1) / bb(i);
    bb(i-1) = bb(i-1) - xMult * cc(i-1);
    u(i-1) = u(i-1) -xMult * u(i);
end

u(2) = u(2) / bb(1);
for i = 3:N
    u(i) = ( u(i) - cc(i) * u(i-1) ) / bb(i);
end

end

