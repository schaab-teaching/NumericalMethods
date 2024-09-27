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
left_bound  = p.u1(G.income(G.grid(:, 1) == 0, :));
right_bound = p.u1(G.income(G.grid(:, 1) == 1, :));
for j = 1:2
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j) * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
end

for iter = 1:p.maxit

    V = V0;
    
    % Value function derivatives:
    [VaF, VaB] = deal(zeros(G.J, 2));
    for j = 1:2
        VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
        VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
    end
    
    % Upwinding:

    % FOC:

    % Generator:

    % Solve the linear system of equations to get Vnew:

    % Compare Vnew to V0 and update:

end