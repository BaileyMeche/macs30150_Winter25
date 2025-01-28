% Question 1b 
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
