%% ----- QUESTION 1
close all; clear all; clc;

%% Part (a)
% Parameter setup
rho = 0.36;   % Capital share
disc = 0.98;  % Discount factor
dep = 0.03;   % Depreciation rate
riskAverse = 2;

% Steady-state capital computation
k_ss = ((rho * disc) / (1 - (1 - dep) * disc))^(1 / (1 - rho));

% Grid configuration
points = 500;
k_low = 0.001 * k_ss;
k_high = 1.5 * k_ss;
k_grid = linspace(k_low, k_high, points)';

% Wealth distribution
w_grid = k_grid.^rho + (1 - dep) * k_grid;

% Initialize storage for value function
v_func = zeros(points, 1);

% Compute initial values
for indx = 1:points
    v_func(indx) = util_fn(max(k_grid(indx)^rho, 1e-10), riskAverse);
end

% ****** Value Function Iteration ******
iter_max = 1000;
tolerance = 1e-10;
converged = false;

for step = 1:iter_max
    % Compute derivatives
    dV = zeros(points, 1);
    dV(1) = (v_func(2) - v_func(1)) / (w_grid(2) - w_grid(1));
    dV(end) = (v_func(end) - v_func(end-1)) / (w_grid(end) - w_grid(end-1));

    for indx = 2:points-1
        dV(indx) = (v_func(indx+1) - v_func(indx-1)) / (w_grid(indx+1) - w_grid(indx-1));
    end
    
    % Ensure positivity
    dV = max(dV, 1e-10);

    % Optimal consumption policy
    c_opt = max((disc * dV .* (rho * k_grid.^(rho-1) + 1 - dep)).^(-1 / riskAverse), 1e-10);
    
    % Implied wealth update
    w_new = c_opt + k_grid;

    % Value function update
    v_new = util_fn(c_opt, riskAverse) + disc * v_func;
    
    % Sorting to maintain monotonicity
    [w_new, indices] = sort(w_new);
    v_new = v_new(indices);
    
    % Remove duplicates
    [w_new, unique_indices] = unique(w_new);
    v_new = v_new(unique_indices);
    
    % Interpolation
    v_updated = interp1(w_new, v_new, w_grid, 'linear', 'extrap');

    % Convergence check
    if max(abs(v_updated - v_func)) < tolerance
        converged = true;
        break;
    end
    
    v_func = v_updated;
end

% ****** Policy Functions Computation ******
ks = zeros(points, 1);
for indx = 1:points
    objective = @(x) max(x^rho + (1 - dep) * x - w_new(indx), -1e10);
    ks(indx) = max(fsolve(objective, k_grid(indx), optimset('Display','off')), 1e-10);
end

% Ensure monotonicity
ks = max(ks, 1e-10);
[ks, sort_indices] = sort(ks);

% Interpolate onto regular grid
k_policy = interp1(ks, k_grid(sort_indices), k_grid, 'linear', 'extrap');
k_policy = max(k_policy, 1e-10);
c_policy = max(k_grid.^rho + (1 - dep) * k_grid - k_policy, 1e-10);
i_policy = k_policy - (1 - dep) * k_grid;

% ****** Graphing with Modified Style ******
figure('Position', [150, 150, 1300, 900]);

subplot(2,2,1);
plot(k_grid, v_func, 'r-', 'LineWidth', 1.8);
title('Value Function Approximation');
xlabel('Capital');
ylabel('Value Function');
grid minor;

subplot(2,2,2);
plot(k_grid, k_policy, 'g-', 'LineWidth', 1.8);
hold on;
plot(k_grid, k_grid, 'k-.', 'LineWidth', 1.2);
title('Capital Policy Mapping');
xlabel('Capital (Current)');
ylabel('Capital (Next Period)');
legend('Policy', '45-degree');
grid minor;

subplot(2,2,3);
plot(k_grid, c_policy, 'b-', 'LineWidth', 1.8);
title('Optimal Consumption Path');
xlabel('Capital');
ylabel('Consumption');
grid minor;

subplot(2,2,4);
plot(k_grid, i_policy, 'm-', 'LineWidth', 1.8);
title('Investment Policy');
xlabel('Capital');
ylabel('Investment');
grid minor;

% ****** Utility Function Definition ******
function utility = util_fn(consumption, gamma)
    consumption = max(consumption, 1e-10);
    if gamma == 1
        utility = log(consumption);
    else
        utility = (consumption.^(1 - gamma) - 1) / (1 - gamma);
    end
end

%% -------- Part (b) ----------
close all; clear all;clc;

capital_share = 0.36;  
discount_factor = 0.98;  
depreciation = 1;       
sigma = 1;  % Equivalent to gamma in GAUSS

% Compute steady-state capital
steady_k = ((capital_share * discount_factor) / (1 - (1 - depreciation) * discount_factor))^(1 / (1 - capital_share));

% ==== Discretization of State Space ====
grid_size = 500;
k_low = 0.001 * steady_k;
k_high = 1.5 * steady_k;
capital_grid = linspace(k_low, k_high, grid_size)';

% Wealth grid construction (GAUSS style)
wealth_grid = capital_grid.^capital_share + (1 - depreciation) * capital_grid;

% ==== Initialize Value Function ====
value_fn = zeros(grid_size, 1);
for idx = 1:grid_size
    value_fn(idx) = util_fn(capital_grid(idx)^capital_share, sigma);
end

% ==== Value Function Iteration ====
max_iterations = 1000;
threshold = 1e-10;
converged = false;

for iteration = 1:max_iterations
    % Compute derivatives (mimicking GAUSS dval)
    d_val = zeros(grid_size, 1);
    d_val(1) = (value_fn(2) - value_fn(1)) / (wealth_grid(2) - wealth_grid(1));
    d_val(end) = (value_fn(end) - value_fn(end-1)) / (wealth_grid(end) - wealth_grid(end-1));

    for idx = 2:grid_size-1
        d_val(idx) = (value_fn(idx+1) - value_fn(idx-1)) / (wealth_grid(idx+1) - wealth_grid(idx-1));
    end

    % Optimal consumption decision (GAUSS style)
    consumption_opt = (discount_factor * d_val .* (capital_share * capital_grid.^(capital_share-1) + 1 - depreciation)).^(-1 / sigma);

    % Compute updated wealth
    wealth_new = consumption_opt + capital_grid;

    % Update value function
    value_updated = util_fn(consumption_opt, sigma) + discount_factor * value_fn;

    % Interpolation to original grid
    value_fn_new = interp1(wealth_new, value_updated, wealth_grid, 'linear', 'extrap');
    
    % Check for convergence
    diff = max(abs(value_fn_new - value_fn));
    if diff < threshold
        converged = true;
        break;
    end
    value_fn = value_fn_new;
end

% ==== Compute Endogenous Capital Policy ====
endogenous_k = zeros(grid_size, 1);
for idx = 1:grid_size
    eqn_fn = @(x) x^capital_share + (1 - depreciation) * x - wealth_new(idx);
    endogenous_k(idx) = fsolve(eqn_fn, capital_grid(idx), optimset('Display', 'off'));
end

% ==== Policy Function Computation ====
k_policy = interp1(endogenous_k, capital_grid, capital_grid, 'linear', 'extrap');
consumption_policy = capital_grid.^capital_share + (1 - depreciation) * capital_grid - k_policy;
investment_policy = k_policy - (1 - depreciation) * capital_grid;

% ==== Analytical Solutions (σ=1, δ=1 Case) ====
k_policy_analytical = capital_share * discount_factor * capital_grid.^capital_share;
c_policy_analytical = capital_grid.^capital_share - k_policy_analytical;
i_policy_analytical = k_policy_analytical;

% ==== Graphical Representation (Modified Style) ====
figure('Position', [150, 150, 1400, 900]);

subplot(2,2,1);
plot(capital_grid, value_fn, 'r-', 'LineWidth', 2);
title('Value Function Approximation');
xlabel('Capital Stock');
ylabel('Value');
grid minor;

subplot(2,2,2);
plot(capital_grid, k_policy, 'g-', 'LineWidth', 2);
hold on;
plot(capital_grid, k_policy_analytical, 'b--', 'LineWidth', 2);
plot(capital_grid, capital_grid, 'k:', 'LineWidth', 1.2);
title('Capital Transition');
xlabel('Capital (t)');
ylabel('Capital (t+1)');
legend('Numerical Policy', 'Analytical Policy', '45-degree Line');
grid minor;

subplot(2,2,3);
plot(capital_grid, consumption_policy, 'm-', 'LineWidth', 2);
hold on;
plot(capital_grid, c_policy_analytical, 'c--', 'LineWidth', 2);
title('Consumption Policy');
xlabel('Capital Stock');
ylabel('Consumption');
legend('Numerical', 'Analytical');
grid minor;

subplot(2,2,4);
plot(capital_grid, investment_policy, 'b-', 'LineWidth', 2);
hold on;
plot(capital_grid, i_policy_analytical, 'r--', 'LineWidth', 2);
title('Investment Policy');
xlabel('Capital');
ylabel('Investment');
legend('Numerical', 'Analytical');
grid minor;
