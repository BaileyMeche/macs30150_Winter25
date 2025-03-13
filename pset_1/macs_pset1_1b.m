% Define the function and its derivative
f = @(x) (x - 10).*exp(-x.^2) + 5;    % Function
df = @(x) exp(-x.^2) - 2.*x.*(x - 10).*exp(-x.^2); % Derivative of f(x)

% Custom Newton-Raphson Method Function
function [x_star, success] = newton_raphson(f, df, x0, epsilon, delta, max_iter)
    % Initialization
    x_k = x0; % Initial guess
    success = false; % To track if a solution is found
    
    for k = 1:max_iter
        % Step 1: Compute iterate
        x_k_plus_1 = x_k - f(x_k) / df(x_k);
        
        % Step 2: Check stopping criterion
        if abs(x_k - x_k_plus_1) <= epsilon * (1 + abs(x_k_plus_1))
            % Step 3: Check if |f(x_k+1)| <= delta
            if abs(f(x_k_plus_1)) <= delta
                x_star = x_k_plus_1;
                success = true;
                return;
            else
                error('Newton-Raphson failed: |f(x_k+1)| > delta.');
            end
        end
        
        % Update x_k for the next iteration
        x_k = x_k_plus_1;
    end
    
    % If max iterations are reached without convergence
    error('Newton-Raphson failed: Maximum iterations reached.');
end

% Parameters
epsilon = 1e-6; % Tolerance for stopping criterion
delta = 1e-6;   % Tolerance for |f(x)|
max_iter = 100; % Maximum number of iterations

% Test with different starting values
starting_points = [-10, -5, -1 , 0,1, 5, 10]; % Different initial guesses
results = [];

fprintf('Newton-Raphson Method Results:\n');
for i = 1:length(starting_points)
    x0 = starting_points(i);
    try
        [root, success] = newton_raphson(f, df, x0, epsilon, delta, max_iter);
        fprintf('Starting point: %.2f -> Root found: x* = %.6f\n', x0, root);
        results = [results; x0, root];
    catch ME
        fprintf('Starting point: %.2f -> Failed: %s\n', x0, ME.message);
    end
end

% Display results
disp('Results (Starting Point -> Root):');
disp(results);
