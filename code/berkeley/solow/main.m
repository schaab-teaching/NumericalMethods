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

param = define_parameters();


%% GRID INITIALIZATION

G = setup_grid(12, 0, param.tmin, param.tmax, 'NamedDims', {1}, 'Names', {'t'});


%% SOLOW MODEL

% Given K0
% K+1 = K + dt*(I - delta K)

K = zeros(G.J, 1);
K(1) = param.K0;

for i = 1:G.J-1
    
    K(i+1) = K(i) + G.dt * (param.s * param.A * K(i)^param.alpha ...
                  - param.delta * K(i));

end


figure; plot(G.t, K);