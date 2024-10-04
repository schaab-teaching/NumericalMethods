function S = stationary(r, G, p)

% pseudo-code
% given r, run VFI to get c / s as well as generator A
% compute g
% evaluate market clearing S


%% VFI

G.income = r * G.a + p.zz;
c0 = G.income;

V0 = p.u(c0) / p.rho;

% State-constrained boundary conditions:
left_bound  = p.u1(G.income(G.a == p.amin, :));
right_bound = p.u1(G.income(G.a == p.amax, :));
for j = 1:2
    BC{1}.left.type = '0'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j) * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
end

% Value function iteration:
[V, A, c, s] = HJB(V0, c0, G, p);


%% KF EQUATION
g = KF(A, G, p);



end
