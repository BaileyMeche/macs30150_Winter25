% Parameters
sigma = 2;       % CRRA parameter
alpha = 0.36;    % Output elasticity of capital
beta = 0.98;     % Discount factor
delta = 0.025;   % Depreciation rate
T = 100;         % Time horizon
k0 = 0.1;        % Initial capital

% Production function
f = @(k) k.^alpha + (1 - delta) * k;
f_prime = @(k) alpha * k.^(alpha - 1) + (1 - delta);

% Utility function and its derivative
u_prime = @(c) c.^(-sigma);

% Nonlinear system of equations to solve
function F = ngm_system(k, T, f, f_prime, u_prime, beta)
    F = zeros(T, 1);
    for t = 1:T-1
        % Euler equation
        c_t = f(k(t)) - k(t + 1);
        c_tp1 = f(k(t + 1)) - k(t + 2);
        F(t) = u_prime(c_t) - beta * u_prime(c_tp1) * f_prime(k(t + 1));
    end
    % Boundary condition at T
    F(T) = k(T + 1); % k_T+1 = 0
end

% Initialize capital array (k_T+1 = 0)
k_guess = linspace(k0, 0, T + 1)';

% Solve the system using fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-8, 'TolX', 1e-8);
k_solution = fsolve(@(k) ngm_system(k, T, f, f_prime, u_prime, beta), k_guess, options);

% Compute consumption
c = zeros(T, 1);
for t = 1:T
    c(t) = f(k_solution(t)) - k_solution(t + 1);
end

% Plot the results
figure;
subplot(2, 1, 1);
plot(0:T, k_solution, 'b-');
xlabel('Time'); ylabel('Capital'); title('Capital Path'); grid on;

subplot(2, 1, 2);
plot(1:T, c, 'r-');
xlabel('Time'); ylabel('Consumption'); title('Consumption Path'); grid on;
