function p = define_parameters()


%% GRID PARAMETERS
p.amin = -1;
p.amax = 10; 

p.l = 7;


%% TUNING PARAMETERS
p.maxit = 100;
p.Delta = 1000;
p.tol = 1e-8;

p.maxitKF = 100;
p.DeltaKF = 1000;
p.tolKF = 1e-8;


%% ECONOMIC PARAMETERS

% Preferences:
p.rho = 0.02;
p.gamma = 2;

p.u = @(x) 1/(1 - p.gamma) * x.^(1-p.gamma);
p.u1 = @(x) x.^(-p.gamma);

% Exogenous constant interest rate:
p.r = 0.01;

% Earnings process:
p.zz = [0.8, 1.2];
p.lambda = 1/3;


end