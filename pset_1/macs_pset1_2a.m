% MATLAB Code for Solving Competitive Equilibrium

% Parameters
sigma_a = 2;   % CRRA coefficient for type a
sigma_b = 3;   % CRRA coefficient for type b
beta = 0.95;   % Discount factor
y0 = 1;        % Endowment in period 0
y1 = 2 * y0;   % Endowment in period 1 (2*y0)

% Define the system of equations as a nested function
function F = equilibrium_system(vars, y0, y1, sigma_a, sigma_b, beta)
    % Variables to solve for
    c0_a = vars(1); % Consumption in period 0 for type a
    c0_b = vars(2); % Consumption in period 0 for type b
    r = vars(3);    % Interest rate

    % Consumption in period 1 based on budget constraints
    c1_a = y1 + (y0 - c0_a) * (1 + r);
    c1_b = y1 + (y0 - c0_b) * (1 + r);

    % Euler equations
    eq1 = (c0_a^(-sigma_a)) - beta * (1 + r) * (c1_a^(-sigma_a));
    eq2 = (c0_b^(-sigma_b)) - beta * (1 + r) * (c1_b^(-sigma_b));

    % Market clearing condition
    eq3 = c0_a + c0_b - 2 * y0;

    % Output the system of equations
    F = [eq1; eq2; eq3];
end

% Initial guess for [c0_a, c0_b, r]
initial_guess = [0.5, 0.5, 0.1];

% Solve the system of equations using fsolve
options = optimoptions('fsolve', 'Display', 'iter');
[solution, fval, exitflag] = fsolve(@(vars) equilibrium_system(vars, y0, y1, sigma_a, sigma_b, beta), initial_guess, options);

% Extract results
c0_a = solution(1);
c0_b = solution(2);
r = solution(3);

% Display results
disp('Numerical Solution:');
fprintf('Consumption in period 0 (type a): %.4f\\n', c0_a);
fprintf('Consumption in period 0 (type b): %.4f\\n', c0_b);
fprintf('Equilibrium interest rate: %.4f\\n', r);

% Verify results
c1_a = y1 + (y0 - c0_a) * (1 + r);
c1_b = y1 + (y0 - c0_b) * (1 + r);
fprintf('Consumption in period 1 (type a): %.4f\\n', c1_a);
fprintf('Consumption in period 1 (type b): %.4f\\n', c1_b);
