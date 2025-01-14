% Parameters
sigma = 1;       % CRRA parameter
alpha = 0.36;    % Output elasticity of capital
beta = 0.98;     % Discount factor
delta = 1;       % Depreciation rate
T = 100;         % Time horizon
k0 = 0.1;        % Initial capital

% Production function
f = @(k) k.^alpha;
f_prime = @(k) alpha * k.^(alpha - 1);

% Utility function and its derivative (log utility for sigma = 1)
u_prime = @(c) 1 ./ c;

% Define the general G function as described in the image
function G = G_function(k_t, k_tp1, k_tp2, f, f_prime, u_prime, beta)
    c_t = f(k_t) - k_tp1;
    c_tp1 = f(k_tp1) - k_tp2;
    G = u_prime(c_t) - beta * u_prime(c_tp1) * f_prime(k_tp1);
end

% Nonlinear system of equations to solve
function F = ngm_system_b(k, T, f, f_prime, u_prime, beta)
    F = zeros(T, 1);
    for t = 1:T-1
        F(t) = G_function(k(t), k(t + 1), k(t + 2), f, f_prime, u_prime, beta);
    end
    % Boundary condition at T
    F(T) = k(T + 1); % k_T+1 = 0
end

% Initialize capital array (k_T+1 = 0)
k_guess = linspace(k0, 0, T + 1)';

% Solve the system using fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-8, 'TolX', 1e-8);
k_solution = fsolve(@(k) ngm_system_b(k, T, f, f_prime, u_prime, beta), k_guess, options);

% Compute consumption
c = zeros(T, 1);
for t = 1:T
    c(t) = f(k_solution(t)) - k_solution(t + 1);
end

% Analytical solution using general case
k_analytical_values = zeros(T + 1, 1);
k_analytical_values(1) = k0;
for t = 1:T
    if t < T
        k_analytical_values(t + 1) = beta * alpha * k_analytical_values(t)^(alpha);
    else
        k_analytical_values(t + 1) = 0; % Ensuring k_T+1 = 0 mathematically at the boundary
    end
end

% Plot the results
figure;
subplot(2, 1, 1);
plot(0:T, k_solution, 'b-', 'LineWidth', 1.5); hold on;
plot(0:T, k_analytical_values, 'r--', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Capital'); title('Capital Path: Numerical vs Analytical');
legend('Numerical', 'Analytical'); grid on;

subplot(2, 1, 2);
plot(1:T, c, 'r-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Consumption'); title('Consumption Path'); grid on;
