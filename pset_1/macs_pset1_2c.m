% Parameters
beta = 0.95;                      % Discount factor
sigma_a = linspace(0.1, 2, 100);  % Grid for sigma_a
sigma_b = linspace(0.1, 2, 100);  % Grid for sigma_b
y0 = 1; y1 = 2;                   % Endowments (y0=1 implies y1=2*y0)
tol = 1e-6;                       % Tolerance for numerical solver

% Preallocate storage for results
r_values = zeros(length(sigma_a), length(sigma_b));
b_a = zeros(length(sigma_a), length(sigma_b));
b_b = zeros(length(sigma_a), length(sigma_b));

% Solve for r, b_a, and b_b
for i = 1:length(sigma_a)
    for j = 1:length(sigma_b)
        % Set sigma_a and sigma_b
        sa = sigma_a(i);
        sb = sigma_b(j);

        % Dynamically update initial guess for r
        r_initial_guess = 0.05; % Default initial guess
        if i > 1 && j > 1
            r_initial_guess = r_values(i-1, j); % Use prior result as guess
        end

        % Define equilibrium function
        eq_func = @(r) compute_equilibrium_r(r, beta, sa, sb, y0, y1);

        % Solve for r using fsolve
        try
            r_values(i, j) = fsolve(eq_func, r_initial_guess, optimset('Display', 'off', 'TolFun', tol));

            % Compute b_a and b_b
            b_a(i, j) = compute_b_a(r_values(i, j), beta, sa, y0, y1);
            b_b(i, j) = -b_a(i, j); % Market clearing condition

            % Check for feasibility and validity
            if ~isreal(b_a(i, j)) || ~isreal(b_b(i, j)) || b_a(i, j) < -y0 || b_a(i, j) > y0
                b_a(i, j) = NaN;
                b_b(i, j) = NaN;
                r_values(i, j) = NaN;
            end
        catch
            % Handle failed cases
            r_values(i, j) = NaN;
            b_a(i, j) = NaN;
            b_b(i, j) = NaN;
        end
    end
end

% Plot b_a and b_b
figure;
subplot(1, 2, 1);
surf(sigma_a, sigma_b, b_a, 'EdgeColor', 'none');
xlabel('\sigma_a'); ylabel('\sigma_b'); zlabel('b_a');
title('Bond Holdings for Type a');
colorbar;

subplot(1, 2, 2);
surf(sigma_a, sigma_b, b_b, 'EdgeColor', 'none');
xlabel('\sigma_a'); ylabel('\sigma_b'); zlabel('b_b');
title('Bond Holdings for Type b');
colorbar;

% Function to compute equilibrium conditions for r
function F = compute_equilibrium_r(r, beta, sa, sb, y0, y1)
    % Solve b_a from FOC
    b_a = compute_b_a(r, beta, sa, y0, y1);

    % Solve b_b using market clearing
    b_b = -b_a;

    % Type a consumption
    c0_a = y0 - b_a;
    c1_a = y1 + b_a * (1 + r);

    % Type b consumption
    c0_b = y0 - b_b;
    c1_b = y1 + b_b * (1 + r);

    % FOC for type a
    foc_a = c0_a^(-sa) - beta * (1 + r) * c1_a^(-sa);

    % FOC for type b
    foc_b = c0_b^(-sb) - beta * (1 + r) * c1_b^(-sb);

    % Return combined residuals
    F = foc_a + foc_b;
end

% Function to compute b_a for a given r
function b_a = compute_b_a(r, beta, sigma, y0, y1)
    % Solve FOC for type a
    eq_b = @(b) (y0 - b)^(-sigma) - beta * (1 + r) * (y1 + b * (1 + r))^(-sigma);
    b_guess = 0; % Initial guess for b_a
    
    % Solve for b_a using fsolve
    b_a = fsolve(eq_b, b_guess, optimset('Display', 'off'));

    % Ensure feasibility of bond holdings
    % Bond holdings must be within bounds [-y0, y0]
    if ~isreal(b_a) || b_a < -y0 || b_a > y0
        error('Bond holdings are infeasible or not real.');
    end

    % Comment: 
    % 1. This function computes the bond holdings for type a individuals given the interest rate r.
    % 2. It solves the first-order condition (FOC) derived from the optimization problem.
    % 3. Feasibility checks ensure b_a remains within the economic bounds of the problem.
end
