%% Problem 1 (b) 
function ngm_collocation()
    % Parameters
    alpha = 0.4;    % Capital share
    beta = 0.95;    % Discount factor
    delta = 1;   % Depreciation rate 0.03
    sigma = 1;      % Risk aversion 2
    n = 7;          % Number of nodes

    % Find steady state capital
    k_ss = steady_state(alpha, beta, delta);
    
    % Define approximation range
    k_min = 0.1 * k_ss;
    k_max = 1.2 * k_ss;
    
    % Initial guess for coefficients (adjust length to n)
    a0 = transpose(linspace(0.1, 0.01, n));  % Adaptive initial guess
    
    % Solve for coefficients
    options = optimset('Display', 'iter', 'MaxFunEvals', 3000, 'MaxIter', 1000);  
    [a_sol, fval, exitflag] = fsolve(@(a) residual_function(a, n, k_min, k_max, alpha, beta, delta, sigma), a0, options);
    
    % Plot results
    plot_solution(a_sol, n, k_min, k_max, k_ss, alpha, delta);
end

function k_ss = steady_state(alpha, beta, delta)
    % Compute steady state capital
    k_ss = (alpha * beta /1)^(1/(1 - alpha));
end

function R = residual_function(a, n, k_min, k_max, alpha, beta, delta, sigma)
    % Get Chebyshev nodes
    z = cos((2*(1:n)' - 1)*pi/(2*n));
    k = ((z + 1) * (k_max - k_min)/2) + k_min;
    
    % Compute next period capital using policy function
    k_next = policy_function(k, a, k_min, k_max);
    k_next_next = policy_function(k_next, a, k_min, k_max);
    
    % Production and consumption
    y = k.^alpha;
    c = y + (1 - delta) * k - k_next;
    c_next = k_next.^alpha + (1 - delta) * k_next - k_next_next;
    
    % Euler equation residuals
    R = c.^(-sigma) - beta * c_next.^(-sigma) .* (alpha * k_next.^(alpha - 1) + 1 - delta);
end

function k_next = policy_function(k, a, k_min, k_max)
    % Map k to [-1, 1]
    x = 2 * (k - k_min) / (k_max - k_min) - 1;
    
    % Evaluate Chebyshev polynomials
    T = ones(length(k), length(a));  % Initialize T
    if length(a) > 1
        T(:, 2) = x;
    end
    for j = 3:length(a)
        T(:, j) = 2 * x .* T(:, j - 1) - T(:, j - 2);
    end
    
    %for i=0:length(x)-1
    %

    % Compute policy function
    k_next = sum(T .* a', 2);  % Correct matrix multiplication
end

function plot_solution(a, n, k_min, k_max, k_ss, alpha, delta)
    % Create grid for plotting
    k_grid = linspace(k_min, k_max, 200)';  % Make column vector
    
    % Compute policy function
    k_next = policy_function(k_grid, a, k_min, k_max);
    
    % Plot results
    figure;
    plot(k_grid, k_next, 'b-', 'LineWidth', 2);
    hold on;
    plot(k_grid, k_grid, 'r--');  % 45-degree line
    plot([k_ss k_ss], [k_min k_max], 'g:', 'LineWidth', 1.5);  % Steady state
    xlabel('Current Capital (k)');
    ylabel('Next Period Capital (k'')');
    title('Policy Function for Capital');
    legend('Policy Function', '45-degree line', 'Steady State');
    grid on;
end

%% Problem 2 (d-g)

% Parameters
beta = 0.95;
alpha = 0.36;
delta = 0.025;
sigma = 2;

% Find steady state k and c first
k_ss = ((1/beta - (1-delta))/(alpha))^(1/(alpha-1));
c_ss = k_ss^alpha - delta*k_ss;

% Calculate derivatives at steady state
f_prime = alpha * k_ss^(alpha-1);
f_double_prime = alpha * (alpha-1) * k_ss^(alpha-2);
u_double_prime = -sigma * c_ss^(-sigma-1);
u_prime = c_ss^(-sigma);

% Construct matrix A
A = [1/beta, -1;
     -(u_prime*f_double_prime/u_double_prime), 1 + (beta*u_prime*f_double_prime/u_double_prime)];

%%% (Part e)
% Find eigenvalues and eigenvectors
[V, D] = eig(A);

% Policy functions coefficients
V_inv = inv(V);
v21 = V_inv(2,1);
v22 = V_inv(2,2);
consumption_policy = v21/v22;
capital_policy = min(diag(D)); % lambda1

% Simulate transition dynamics (Part e)
T = 100;
k0 = 0.1 * k_ss;

k_path = zeros(T+1, 1);
c_path = zeros(T+1, 1);
k_path(1) = k0;

for t = 1:T
    k_dev = k_path(t) - k_ss;
    c_path(t) = c_ss + consumption_policy * k_dev;
    k_path(t+1) = k_ss + capital_policy * k_dev;
end
c_path(T+1) = c_ss + consumption_policy * (k_path(T+1) - k_ss);

% Figure for Part (e)
% Transition simulation using policy functions
figure;
subplot(2,1,1);
plot(0:T, k_path, 'r', 'DisplayName', 'Capital Path (Part e)');
title('Capital Path (Part e: Policy Function Simulation)');
ylabel('Capital');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T, c_path, 'r', 'DisplayName', 'Consumption Path (Part e)');
title('Consumption Path (Part e: Policy Function Simulation)');
ylabel('Consumption');
xlabel('Time');
legend show;
grid on;

% Run Dynare nonlinear simulation (Part f)
dynare PQuestion2F.mod noclearall;

% Extract Dynare results
dynare_k = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_c = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';

% Figure for Part (f)
% Compare MATLAB Part (e) with Dynare nonlinear results
figure;
subplot(2,1,1);
plot(0:T, k_path, 'r', 'DisplayName', 'MATLAB Capital Path (Part e)');
hold on;
plot(0:(length(dynare_k)-1), dynare_k, 'b--', 'DisplayName', 'Dynare Capital Path (Part f)');
title('Capital Path (Part f: Nonlinear Simulation Comparison)');
ylabel('Capital');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T, c_path, 'r', 'DisplayName', 'MATLAB Consumption Path (Part e)');
hold on;
plot(0:(length(dynare_c)-1), dynare_c, 'b--', 'DisplayName', 'Dynare Consumption Path (Part f)');
title('Consumption Path (Part f: Nonlinear Simulation Comparison)');
ylabel('Consumption');
xlabel('Time');
legend show;
grid on;

% Run Dynare linear simulation (Part g)
dynare PQuestion2G.mod noclearall;

% Extract Dynare results
dynare_k_linear = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_c_linear = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';

% Figure for Part (g)
% Compare MATLAB Part (e) with Dynare linear results
figure;
subplot(2,1,1);
plot(0:T, k_path, 'r', 'DisplayName', 'MATLAB Capital Path (Part e)');
hold on;
plot(0:(length(dynare_k_linear)-1), dynare_k_linear, 'b--', 'DisplayName', 'Dynare Capital Path (Part g)');
title('Capital Path (Part g: Linear Approximation Comparison)');
ylabel('Capital');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T, c_path, 'r', 'DisplayName', 'MATLAB Consumption Path (Part e)');
hold on;
plot(0:(length(dynare_c_linear)-1), dynare_c_linear, 'b--', 'DisplayName', 'Dynare Consumption Path (Part g)');
title('Consumption Path (Part g: Linear Approximation Comparison)');
ylabel('Consumption');
xlabel('Time');
legend show;
grid on;

fprintf('Policy functions:\n');
fprintf('c_t - c = %.4f(k_t - k)\n', consumption_policy);
fprintf('k_{t+1} - k = %.4f(k_t - k)\n', capital_policy);
