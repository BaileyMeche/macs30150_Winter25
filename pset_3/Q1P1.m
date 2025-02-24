
alpha = 0.4;
beta = 0.95;
delta = 1;
sigma = 2;
n = 7;

k_ss = ((alpha*beta)/(1-beta*(1-delta)))^(1/(1-alpha));
k_min = 0.1*k_ss;
k_max = 1.2*k_ss;

% A good guess for n roots
a0 = transpose(linspace(0.1, 0.01, n));

a_sol = fsolve(@(a) Resid(a, n, k_min, k_max, alpha, beta, delta, sigma), a0);

k_grid = linspace(k_min, k_max, 200)';    
k_next = policy(k_grid, a_sol, k_min, k_max);
plot(k_grid, k_next);
hold on;
plot(k_grid, k_grid);
plot([k_ss k_ss], [k_min k_max]);
xlabel('k');
ylabel('k''');
title('Capital Policy Function delta = 1');
legend('Policy Function', 'k=k''', 'kss');
grid on;

function R = Resid(a, n, k_min, k_max, alpha, beta, delta, sigma)
    % Find and rescale n roots
    z = cos((2*(1:n)' - 1)*pi/(2*n));
    k = ((z + 1) * (k_max - k_min)/2) + k_min;

    k_next = policy(k, a, k_min, k_max);
    k_next_next = policy(k_next, a, k_min, k_max);

    y = k.^alpha;
    c = y + (1 - delta) * k - k_next;
    c_next = k_next.^alpha + (1 - delta) * k_next - k_next_next;

    R = c.^(-sigma) - beta * c_next.^(-sigma) .* (alpha * k_next.^(alpha - 1) + 1 - delta);
end

function k_next = policy(k, a, k_min, k_max)
    x = 2 * (k - k_min) / (k_max - k_min) - 1;
    
    T = ones(length(k), length(a));
    T(:, 2) = x;

    for j = 3:length(a)
        T(:, j) = 2 * x .* T(:, j - 1) - T(:, j - 2);
    end
    
    k_next = sum(T .* a', 2);
end