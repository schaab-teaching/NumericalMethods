%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PREAMBLE
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 
clc

addpath(genpath('../../lib'));


%% DEFINE PARAMETERS

p = define_parameters();


%% INITIALIZE GRID

G = setup_grid(p.l, 0, p.amin, p.amax, 'NamedDims', {1}, 'Names', {'a'});
figure; scatter(G.a, 0*G.a, 'filled');


%% VFI

G.income = p.r * G.a + p.zz;
c0 = G.income;

V0 = p.u(c0) / p.rho;

figure; scatter(G.a, V0);


% State-constrained boundary conditions:
left_bound  = p.u1(G.income(G.a == p.amin, :));
right_bound = p.u1(G.income(G.a == p.amax, :));
for j = 1:2
    BC{1}.left.type = '0'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j) * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
end

V = V0;

for iter = 1:p.maxit
    
    % Value function derivatives:
    [VaF, VaB] = deal(zeros(G.J, 2));
    for j = 1:2
        VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
        VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
    end
    Va0 = c0.^(-p.gamma);

    % VaF(G.a == p.amax, :) = p.u1(G.income(G.a == p.amax, :));
    % VaB(G.a == p.amin, :) = p.u1(G.income(G.a == p.amin, :));
    
    % Upwinding:
    cF = VaF.^(-1/p.gamma);
    cB = VaB.^(-1/p.gamma);

    sF = G.income - cF;
    sB = G.income - cB;

    IF = sF > 0;
    IB = (sB < 0) .* ~IF;
    I0 = ~IF .* ~IB;

    s = sF .* IF + sB .* IB;
    c = cF .* IF + cB .* IB + c0 .* I0;
    u = p.u(c);

    % Generator:
    Aa{1} = FD_operator(G, s(:, 1), zeros(G.J, 1), 1, '1');
    Aa{2} = FD_operator(G, s(:, 2), zeros(G.J, 1), 1, '2');
    
    Az = [-speye(G.J)*p.lambda,  speye(G.J)*p.lambda; ...
           speye(G.J)*p.lambda, -speye(G.J)*p.lambda];

    A = blkdiag(Aa{1}, Aa{2}) + Az;

    % Solve the linear system of equations to get Vnew:
    % B * Vnew = u + 1/Delta * Vold
    % A x = b => x = A\b
    B = (p.rho + 1/p.Delta) * speye(G.J*2) - A;
    V_new = B \ (u(:) + 1/p.Delta * V(:));

    % Compare Vnew to V0 and update:
    V_change = V_new - V(:);
    V = reshape(V_new, [G.J, 2]);
    
    dist = max(max(abs(V_change)));
    if dist < p.tol; break; end

    if mod(iter, 1) == 0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end

end

