% Parameters
beta = 0.95; % Discount factor
y0 = 1;      % First-period endowment
y1 = 2;      % Second-period endowment
sigma_values = 0.1:0.1:2; % Values for sigma_a and sigma_b, evenly spaced from 0.1 to 2

% Preallocate results
num_sigma = length(sigma_values);
r_results = zeros(num_sigma, num_sigma);
c0_a_results = zeros(num_sigma, num_sigma);

% Define the system of equations as a function of sigma_a, sigma_b
for i = 1:num_sigma
    for j = 1:num_sigma
        sigma_a = sigma_values(i);
        sigma_b = sigma_values(j);

        % Define the system of equations
        equations = @(x) [
            (x(1))^(-sigma_a) - beta * (1 + x(2)) * (y1 + (y0 - x(1)) * (1 + x(2)))^(-sigma_a);
            (2 * y0 - x(1))^(-sigma_b) - beta * (1 + x(2)) * (y1 + (y0 - (2 * y0 - x(1))) * (1 + x(2)))^(-sigma_b)
        ];

        % Initial guesses for c_{i0}^a and r
        x0 = [0.5, 0.05];

        % Solve the system of equations
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8);
        try
            solution = fsolve(equations, x0, options);

            % Store results
            c0_a_results(i, j) = solution(1);
            r_results(i, j) = solution(2);
        catch
            % Handle cases where fsolve does not converge
            c0_a_results(i, j) = NaN;
            r_results(i, j) = NaN;
        end
    end
end

% Visualization of results
figure;
surf(sigma_values, sigma_values, r_results);
title('Equilibrium Interest Rate (r) as a Function of \sigma_a and \sigma_b');
xlabel('\sigma_a');
ylabel('\sigma_b');
zlabel('Equilibrium Interest Rate (r)');
colorbar;

figure;
surf(sigma_values, sigma_values, c0_a_results);
title('Equilibrium First-Period Consumption (c_{i0}^a)');
xlabel('\sigma_a');
ylabel('\sigma_b');
zlabel('Consumption c_{i0}^a');
colorbar;

% Save results to a CSV file for external analysis
results_table = array2table([reshape(repmat(sigma_values', 1, num_sigma), [], 1), ...
                             reshape(repmat(sigma_values, num_sigma, 1), [], 1), ...
                             reshape(r_results, [], 1), ...
                             reshape(c0_a_results, [], 1)], ...
                             'VariableNames', {'sigma_a', 'sigma_b', 'r', 'c0_a'});

writetable(results_table, 'equilibrium_results.csv');

% Output results as text
output_file = fopen('equilibrium_results.txt', 'w');
fprintf(output_file, 'Equilibrium Results:\n');
for i = 1:num_sigma
    for j = 1:num_sigma
        fprintf(output_file, 'sigma_a: %.2f, sigma_b: %.2f, r: %.4f, c0_a: %.4f\n', sigma_values(i), sigma_values(j), r_results(i, j), c0_a_results(i, j));
    end
end
fclose(output_file);

fprintf('Results saved to "equilibrium_results.csv" and "equilibrium_results.txt"\n');
