%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HW3Q3_cd.m
%   Problem 3 (c) & (d): Neoclassical Growth with Variable Labor
%   Nonlinear vs. Linear Approximation using Dynare
%
%   We assume a central planner solves:
%     max_{c_t, l_t, k_{t+1}}  sum_{t=0 to infty} beta^t [ ln(c_t) + gamma(1-l_t) ]
%   subject to
%     k_{t+1} + c_t = k_t^alpha * l_t^{1-alpha} + (1-delta)*k_t.
%   We set alpha=0.36, beta=0.95, delta=0.025, and choose l_ss=1/3 => solve for
%   (k_ss, c_ss, gamma) that satisfies the steady state conditions from part (b).
%
%   Part (c): Nonlinear time paths from 10% of steady-state capital
%   Part (d): Linear approximation, compare paths with (c)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% 1. Parameters
alpha  = 0.36;
beta   = 0.95;
delta  = 0.025;
ell_ss = 1/3;   % We impose steady-state labor = 1/3

% 2. Solve for the steady state (k_ss, c_ss, gamma)
%    from the three equations in part (b):
%       (i)   alpha * k_ss^(alpha-1) * ell_ss^(1-alpha) + 1 - delta = 1/beta
%       (ii)  c_ss = k_ss^alpha * ell_ss^(1-alpha) - delta * k_ss
%       (iii) c_ss = (1-alpha)*k_ss^alpha / [ gamma * ell_ss^alpha ]

% (i) => alpha * k_ss^(alpha - 1) * ell_ss^(1 - alpha) = 1/beta - (1-delta).
lhs = (1/beta) - (1 - delta);
% => k_ss^(alpha-1) = lhs / [ alpha * ell_ss^(1-alpha) ]
k_ss = ( lhs / ( alpha * ell_ss^(1-alpha) ) )^(1/(alpha-1));  % watch sign of (alpha-1) < 0

% (ii):
c_ss = k_ss^alpha * ell_ss^(1 - alpha) - delta * k_ss;

% (iii) => gamma = ( (1-alpha)*k_ss^alpha ) / ( c_ss * ell_ss^alpha )
gamma = ( (1 - alpha)*k_ss^alpha ) / ( c_ss * ell_ss^alpha );

fprintf('STEADY STATE RESULTS (Part b):\n');
fprintf('k_ss    = %8.4f\n', k_ss);
fprintf('c_ss    = %8.4f\n', c_ss);
fprintf('ell_ss  = %8.4f\n', ell_ss);
fprintf('gamma   = %8.4f\n\n', gamma);

%% 3. Write a Dynare .mod file in memory for Nonlinear + Linear Simulations
% We'll name it "HW3Q3cd_temp.mod". It will define the system:
%   1) Euler:  (1/c) = beta*(1/c(+1)) * [ alpha*k(+1)^(alpha-1)*l(+1)^(1-alpha) + 1-delta ]
%   2) Labor supply:  c - ( (1-alpha)*k^alpha / (gamma*l^alpha ) ) = 0
%   3) Resource: k(+1) + c - [ k^alpha*l^(1-alpha) + (1-delta)*k ] = 0

tempModFile = 'HW3Q3cd_temp.mod';
fid = fopen(tempModFile, 'w');

fprintf(fid, 'var k c l;\n');
fprintf(fid, 'predetermined_variables k;\n');
fprintf(fid, 'parameters alpha beta delta gamma;\n\n');

fprintf(fid, 'alpha = %.6f;\n', alpha);
fprintf(fid, 'beta  = %.6f;\n', beta);
fprintf(fid, 'delta = %.6f;\n', delta);
fprintf(fid, 'gamma = %.6f;\n', gamma);

fprintf(fid, '\nmodel;\n');
fprintf(fid, '// 1) Euler\n');
fprintf(fid, '(1/c) = beta*(1/c(+1))*( alpha*k(+1)^(alpha-1)*l(+1)^(1-alpha) + 1 - delta );\n\n');
fprintf(fid, '// 2) Labor supply\n');
fprintf(fid, 'c - ( (1-alpha)*k^alpha/( gamma*l^alpha ) ) = 0;\n\n');
fprintf(fid, '// 3) Resource constraint\n');
fprintf(fid, 'k(+1) + c - ( k^alpha*l^(1-alpha) + (1-delta)*k ) = 0;\n');
fprintf(fid, 'end;\n\n');

% initval: start at 10% of k_ss, near c_ss, near l_ss
% (No reason we cannot start l near the SS as well. 
%  If you want, you could shift them a bit.)
fprintf(fid, 'initval;\n');
fprintf(fid, '  k = %.6f;\n', 0.1*k_ss);
fprintf(fid, '  c = %.6f;\n', 0.5*c_ss);
fprintf(fid, '  l = %.6f;\n', 0.5*ell_ss);
fprintf(fid, 'end;\n\n');

% endval: converge to SS
fprintf(fid, 'endval;\n');
fprintf(fid, '  k = %.6f;\n', k_ss);
fprintf(fid, '  c = %.6f;\n', c_ss);
fprintf(fid, '  l = %.6f;\n', ell_ss);
fprintf(fid, 'end;\n\n');

fprintf(fid, 'resid;\n\n');
% Perfect foresight setup for 50-60 periods
simulationPeriods = 200;
fprintf(fid, 'perfect_foresight_setup(periods=%d);\n', simulationPeriods);
fprintf(fid, 'perfect_foresight_solver;\n\n');

% Request saving the variables k, c, l in separate series
fprintf(fid, 'rplot k c l;\n');

% For convenience, we next do linear_approximation in the same mod file 
% by re-running the perfect_foresight_setup + solver, but with "linear_approximation"
fprintf(fid, '\n// Now do linear approximation\n');
fprintf(fid, 'perfect_foresight_setup(periods=%d);\n', simulationPeriods);
fprintf(fid, 'perfect_foresight_solver(linear_approximation);\n');
fprintf(fid, 'rplot k c l;\n');

fclose(fid);

%% 4. Run Dynare on that .mod file
clearvars -except alpha beta delta gamma k_ss c_ss ell_ss tempModFile simulationPeriods
dynare(tempModFile, 'noclearall');

% 5. Parse Solutions
allVars = M_.endo_names;   % Should be {'k','c','l'} in that order
kIndex  = find(strcmp(allVars,'k'));
cIndex  = find(strcmp(allVars,'c'));
lIndex  = find(strcmp(allVars,'l'));

% Nonlinear solution chunk:
nonlin_k = oo_.endo_simul(kIndex, 1:(simulationPeriods+1));
nonlin_c = oo_.endo_simul(cIndex, 1:(simulationPeriods+1));
nonlin_l = oo_.endo_simul(lIndex, 1:(simulationPeriods+1));

% Linear approximation chunk:
% Run Dynare for Linear Approximation
dynare PQuestion3D.mod noclearall;

% Extract Linear Approximation Results
lin_k = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
lin_c = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';
lin_l = oo_.endo_simul(strcmp(M_.endo_names, 'l'), :)';

% Ensure timeVec matches the length of results
timeVec = 0:(length(nonlin_k) - 1);

% Trim lin_k, lin_c, lin_l to match the length of timeVec
lin_k_trimmed = lin_k(1:length(timeVec));
lin_c_trimmed = lin_c(1:length(timeVec));
lin_l_trimmed = lin_l(1:length(timeVec));

% Overlay Nonlinear and Linear Results for Capital
figure('Name', 'Capital: Nonlinear vs Linear');
plot(timeVec, nonlin_k, 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear'); hold on;
plot(timeVec, lin_k_trimmed, 'r--', 'LineWidth', 2, 'DisplayName', 'Linear');
title('Capital: Nonlinear vs. Linear Approximation');
xlabel('Time');
ylabel('k_t');
legend('Location', 'Best');
grid on;

% Restrict timeVec and corresponding results 
timeRestricted = timeVec(2:end); % Exclude time 0
nonlin_c_restricted = nonlin_c(2:end); % Consumption: Nonlinear 
lin_c_restricted = lin_c_trimmed(2:length(timeRestricted)+1); % Consumption: Linear 

nonlin_l_restricted = nonlin_l(2:end); % Labor: Nonlinear
lin_l_restricted = lin_l_trimmed(2:length(timeRestricted)+1); % Labor: Linear 

% Overlay Nonlinear and Linear Results for Consumption 
figure('Name', 'Consumption: Nonlinear vs Linear (1 to 60)');
plot(timeRestricted, nonlin_c_restricted, 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear'); hold on;
plot(timeRestricted, lin_c_restricted, 'r--', 'LineWidth', 2, 'DisplayName', 'Linear');
title('Consumption: Nonlinear vs. Linear Approximation');
xlabel('Time');
ylabel('c_t');
legend('Location', 'Best');
grid on;

% Overlay Nonlinear and Linear Results for Labor (Periods 1 to 60)
figure('Name', 'Labor: Nonlinear vs Linear');
plot(timeRestricted, nonlin_l_restricted, 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear'); hold on;
plot(timeRestricted, lin_l_restricted, 'r--', 'LineWidth', 2, 'DisplayName', 'Linear');
title('Labor: Nonlinear vs. Linear Approximation');
xlabel('Time');
ylabel('l_t');
legend('Location', 'Best');
grid on;
