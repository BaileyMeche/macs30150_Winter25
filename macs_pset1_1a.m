% Define the function
f = @(x) (x - 10).*exp(-x.^2) + 5;

% Custom Bisection Method Function (No Sign Check)
function x_star = custom_bisection(f, x_low, x_high, epsilon, delta)
    while true
        % Step 1: Compute midpoint
        x_mid = (x_low + x_high) / 2;
        
        % Step 2: Update bounds
        if f(x_low) * f(x_mid) < 0
            x_high = x_mid; % Root lies in [x_low, x_mid]
        else
            x_low = x_mid; % Root lies in [x_mid, x_high]
        end
        
        % Step 3: Check stopping criteria
        if abs(x_high - x_low) <= epsilon * (1 + abs(x_low) + abs(x_high)) || abs(f(x_mid)) <= delta
            x_star = x_mid; % Root found
            return;
        end
    end
end

% Function to find multiple zeroes
function zeroes = find_multiple_zeroes(f, x_start, x_end, epsilon, delta, step_size)
    zeroes = []; % Array to store roots
    search_range = x_start:step_size:x_end; % Divide the range into small intervals
    
    for i = 1:(length(search_range) - 1)
        x_low = search_range(i);
        x_high = search_range(i + 1);
        
        % Check if the interval contains a potential root
        if f(x_low) * f(x_high) <= 0 % Check for sign change
            % Use the bisection method to find the root
            try
                root = custom_bisection(f, x_low, x_high, epsilon, delta);
                
                % Check if the root is already found (to avoid duplicates)
                if isempty(zeroes) || all(abs(zeroes - root) > epsilon)
                    zeroes = [zeroes; root]; % Append root to the list
                end
            catch
                % If bisection fails, continue to the next interval
                continue;
            end
        end
    end
end

% Parameters for finding zeroes
x_start = -10; % Start of the range
x_end = 10; % End of the range
epsilon = 1e-6; % Stopping criterion for interval width
delta = 1e-6; % Stopping criterion for function value
step_size = 0.5; % Step size for dividing the range

% Call the function to find multiple zeroes
zeroes = find_multiple_zeroes(f, x_start, x_end, epsilon, delta, step_size);

% Display the results
fprintf('Zeroes of the function:\n');
disp(zeroes);

% Plot the function and zeroes
x_vals = linspace(x_start, x_end, 1000);
f_vals = f(x_vals);

figure;
plot(x_vals, f_vals, 'b-', 'LineWidth', 2);
hold on;
yline(0, 'r--', 'LineWidth', 1.5); % Zero line
plot(zeroes, zeros(size(zeroes)), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Mark zeroes
xlabel('x');
ylabel('f(x)');
title('Function and Zeroes');
grid on;
legend('f(x)', 'y=0', 'Zeroes');


