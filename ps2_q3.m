% Parameters
A = 1; s = 0.3; alpha = 0.36; delta = 0.08; np = 0.02; k0 = 0.01; T = 100; n = 6;

% Define the analytical solution
k_true = @(t) ((k0^(1-alpha) - s*A/(np + delta)) * exp(-(1-alpha)*(np + delta)*t) ...
              + s*A/(np + delta))^(1/(1-alpha));

% Rescale time domain to [0, T] and compute collocation points
k = 1:n;
x_k = cos((2*k - 1) * pi / (2*n)); % Chebyshev points on [-1, 1]
t_k = (x_k + 1) * (T / 2); % Map to [0, T]

% Define symbolic variables and the polynomial approximation
syms t a1 a2 a3 a4 a5 a6
a = [a1, a2, a3, a4, a5, a6];
k_hat = k0 + sum(a .* t.^(1:n)); % Polynomial approximation for k(t)

% Define the error function (squared differences)
error_func = @(a_vals) arrayfun(@(t_val) ...
    (k0 + sum(a_vals .* t_val.^(1:n)) - k_true(t_val))^2, linspace(0, T, 100));

% Sum of squared errors over the defined range
squared_error_func = @(a_vals) sum(error_func(a_vals));

% Initial guess for coefficients
initial_guess = [0.0562, 0.0544, -0.0284, 0.0197, -0.0100, 0.002];

% Solve using lsqnonlin
options = optimoptions('lsqnonlin', 'Display', 'iter', 'Diagnostics', 'on', ...
    'TolFun', 1e-8, 'TolX', 1e-8);
[solution, resnorm, residual, exitflag, output] = lsqnonlin(...
    @(a_vals) sqrt(error_func(a_vals)), initial_guess, [], [], options);

% Display results
disp('Solution coefficients (a1 to a6):');
disp(solution);
disp('Residual norm:');
disp(resnorm);

% Generate points for plotting
t_vals = linspace(0, T, 100); % Time domain for comparison
k_approx_vals = k0 + sum(solution' .* t_vals.^((1:n)'), 1); % Polynomial approximation
k_true_vals = arrayfun(k_true, t_vals); % Analytical solution

% Plot the results
figure;
plot(t_vals, k_approx_vals, 'b', 'LineWidth', 2); % Polynomial approximation
hold on;
plot(t_vals, k_true_vals, 'r--', 'LineWidth', 2); % Analytical solution
hold off;
legend('Polynomial Approximation (\hat{k}(t))', 'Analytical Solution (k_{true}(t))', ...
       'Location', 'Best');
xlabel('Time t');
ylabel('Capital k(t)');
title('Comparison of Polynomial Approximation and Analytical Solution');
grid on;

