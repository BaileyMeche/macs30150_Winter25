%------- problem3.m

clear; clc; close all;

%---- Parameters for Parts (a) & (b)
alpha = 0.36;         % Production elasticity
beta  = 0.90;         % Discount factor (baseline for parts a & b)
sigma = 0.01;         % Std. dev. of shock ln(z)
N     = 20;           % Number of grid points
m     = 3;            % ±3 standard deviations for the grid

% To use tauchen.m, we need to match the AR(1) form:
%       x_{t+1} = (1-ρ)*μ + ρ x_t + ε_{t+1},
% with ρ = α and μ chosen so that the constant term ln(αβ) appears.
% Since ln(k_{t+1}) = α ln(k_t) + ln(αβ) + ε, we set:
mu_tauchen = log(alpha * beta) / (1 - alpha);

%% (a) Finite-State Approximation using Tauchen's Method
% [ incomplete ] 

%% (b) Compute the Ergodic (Stationary) Distribution
% Find π such that π'Π = π'. One common way is to iterate on an initial guess.
pi_stationary = ones(N,1) / N;  % Start with a uniform distribution
tol = 1e-8;
err = 1;
while err > tol
    pi_new = Pi' * pi_stationary;
    err = max(abs(pi_new - pi_stationary));
    pi_stationary = pi_new;
end

% --- Interpretation for Part (b) ---
% The vector pi_stationary is the ergodic (stationary) distribution for ln(k),
% which, when exponentiated, gives the distribution over the capital stock.
fprintf('\n--- Part (b): Ergodic Distribution of Capital ---\n');
fprintf('Ergodic distribution (π):\n');
disp(pi_stationary);

% Plot the ergodic distribution against the actual capital levels (levels = exp(ln(k))).
figure;
bar(exp(kGrid), pi_stationary, 'FaceColor', [0.2 0.6 0.5]);
xlabel('Capital Stock, k (levels)');
ylabel('Ergodic Probability');
title('Ergodic Distribution of Capital');
grid on;

%% (c) Compute the Expected Capital Stock for Different Beta Values
% For β = 0.95, 0.97, 0.99, compute the expected value of capital using π.
% Also compute the deterministic steady state: k* = (αβ)^(1/(1-α))
beta_values = [0.95, 0.97, 0.99];
expected_capital = zeros(length(beta_values),1);
deterministic_k = zeros(length(beta_values),1);

fprintf('\n--- Part (c): Expected Capital Stock for Different Beta Values ---\n');
for i = 1:length(beta_values)
    beta_temp = beta_values(i);
    
    % Update the unconditional mean for ln(k) to match ln(αβ_temp)
    mu_tauchen_temp = log(alpha * beta_temp) / (1 - alpha);
    
    % Compute the deterministic steady state capital (levels):
    % k* = (αβ_temp)^(1/(1-α))
    deterministic_k(i) = (alpha * beta_temp)^(1/(1 - alpha));
    
    % Get the grid and transition matrix for the new β value:
    [kGrid_temp, Pi_temp] = tauchen(N, mu_tauchen_temp, alpha, sigma, m);
    
    % Compute the stationary distribution for the new β
    pi_temp = ones(N,1) / N;
    err_temp = 1;
    while err_temp > tol
        pi_new_temp = Pi_temp' * pi_temp;
        err_temp = max(abs(pi_new_temp - pi_temp));
        pi_temp = pi_new_temp;
    end
    
    % The expected capital is computed in levels:
    expected_capital(i) = sum(exp(kGrid_temp) .* pi_temp);
end

% Display the expected capital and the deterministic steady state for each beta.
for i = 1:length(beta_values)
    fprintf('For β = %.2f: Expected Capital = %.4f, Deterministic k* = %.4f\n', ...
        beta_values(i), expected_capital(i), deterministic_k(i));
end

