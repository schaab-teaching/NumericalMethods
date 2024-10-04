function g = KF(A, G, p)

g0 = zeros(G.J, 2);
g0(G.a == p.amin, :) = 1/2/G.da;

g = g0;

AT = A';

for n = 1:p.maxitKF
    
    B = 1/p.DeltaKF * speye(G.J*2) - AT;
    g_new = B \ (1/p.DeltaKF * g(:));

    % Compare Vnew to V0 and update:
    g_change = g_new - g(:);
    g = reshape(g_new, [G.J, 2]);

    dist = max(max(abs(g_change)));
    if dist < p.tolKF; break; end

    if mod(n, 1) == 0, fprintf('KFE: %.i    Remaining Gap: %.2d\n', n, dist); end

end
    
end