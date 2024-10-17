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
% figure; scatter(G.a, 0*G.a, 'filled');


%% STATIONARY EQUILIBRIUM

% guess r 
% write a function that computes S(r)

% Newton over S(r) = 0

r0 = p.r; 

S0 = stationary(r0, G, p);

% newton-solve over stationary(r0)
f = @(x) stationary(x, G, p);

options = optimset('Display', 'iter');
r = fsolve(f, r0, options);

% Rerun stationary with correct interest rate:
[~, ss] = stationary(r, G, p);


%% FIGURES
figure; plot(G.a, ss.g);
figure; plot(G.a, ss.s);