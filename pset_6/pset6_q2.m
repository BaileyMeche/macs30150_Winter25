% pset6-q2.m

clear; close all; clc;

%% PARAMETERS
alpha    = 0.36;
beta     = 0.98;
delta    = 0.03;
gamma    = 2;          % CRRA risk aversion coefficient in u(c)= c^(1-gamma)/(1-gamma)
sigma_eps = 0.1;       % Std. dev. of the shock ε in ln(z)
rho      = 0.7;        % AR(1) coefficient
N_z      = 5;          % Number of states for the shock process
m        = 3;          % Bounds: grid covers [mu - m*sigma, mu + m*sigma]

% For the AR(1): ln(z_t+1) = rho*ln(z_t) + ε, the stationary std. is:
sigma_z = sigma_eps / sqrt(1 - rho^2);

%% (a) TAUCHEN APPROXIMATION AND ERGODIC DISTRIBUTION
disp('--------------------------------------------------');
disp('Problem 2(a): Tauchen approximation and ergodic distribution');

% Use Tauchen’s method to discretize ln(z)
% (Assuming tauchen.m is available; see the tauchen function below.)
[grid_lnz, P] = tauchen(N_z, 0, rho, sigma_z, m);
z_vals = exp(grid_lnz);  % Convert from ln(z) to z

% Compute the ergodic (stationary) distribution by finding the eigenvector
% associated with eigenvalue 1 of the transition matrix.
[V_eig, D] = eig(P');
[~, idx] = min(abs(diag(D) - 1));
ergodic_dist = V_eig(:, idx);
ergodic_dist = ergodic_dist / sum(ergodic_dist);  % Normalize

% Plot the ergodic distribution.
figure;
bar(z_vals, ergodic_dist, 'FaceColor',[0.2 0.6 0.8]);
xlabel('Productivity Shock, z');
ylabel('Ergodic Probability');
title('Ergodic Distribution of z (Tauchen Approximation)');
grid on;

%% (b) VALUE FUNCTION ITERATION (VFI) TO SOLVE THE MODEL
disp('--------------------------------------------------');
disp('Problem 2(b): Solving model via VFI');

% Compute the non-stochastic steady-state capital (k_ss) using the formula:
k_ss = (alpha * beta / (1 - beta*(1-delta)))^(1/(1-alpha));

% Set up a grid for capital with 1000 points in [0.01*k_ss, 4*k_ss]
N_k = 1000;
k_grid = linspace(0.01*k_ss, 4*k_ss, N_k)';

% Preallocate the value function and a policy function index matrix.
V         = zeros(N_k, N_z);    % Each column corresponds to a shock state.
policy_ind = zeros(N_k, N_z);    % Will store the index of optimal k' from the grid.

% Tolerance and maximum iterations for VFI.
tol      = 1e-6;
max_iter = 1000;
iter     = 0;
dist     = inf;

% Define the period utility function (CRRA). For c <= 0, assign a very low value.
u = @(c) (c > 0) .* (c.^(1-gamma))/(1-gamma) + (c <= 0)*(-1e10);

% VFI: For each state (k, z), choose k' to maximize:
%   u(c) + beta * E[V(k',z')]
% where c = z*k^alpha + (1-delta)*k - k'.
while dist > tol && iter < max_iter
    V_new = zeros(size(V));
    for iz = 1:N_z
        z_now = z_vals(iz);
        for ik = 1:N_k
            current_k = k_grid(ik);
            % For each candidate k', compute consumption from the resource constraint.
            cons = z_now * current_k^alpha + (1-delta)*current_k - k_grid;
            % Compute current period utility.
            util = u(cons);
            % Compute expected continuation value:
            EV = 0;
            for iz_next = 1:N_z
                EV = EV + P(iz, iz_next) * V(:, iz_next);
            end
            % Total return from choosing each possible k'.
            total_return = util + beta * EV;
            % Choose the k' that gives the highest value.
            [V_new(ik, iz), policy_ind(ik, iz)] = max(total_return);
        end
    end
    dist = max(max(abs(V_new - V)));
    V = V_new;
    iter = iter + 1;
    if mod(iter,50)==0
        fprintf('Iteration %d, distance = %e\n', iter, dist);
    end
end
fprintf('VFI converged in %d iterations with distance = %e\n', iter, dist);

% Extract the policy function for capital (k') using the index.
policy_k = k_grid(policy_ind);

% Plot the value function and policy function for each shock state.
figure;
subplot(1,2,1);
hold on;
for iz = 1:N_z
    plot(k_grid, V(:, iz), 'LineWidth', 2);
end
xlabel('Capital, k');
ylabel('Value Function, V(k,z)');
title('Value Function (VFI)');
legend(arrayfun(@(z) sprintf('z = %.2f', z), z_vals, 'UniformOutput', false));
grid on;

subplot(1,2,2);
hold on;
for iz = 1:N_z
    plot(k_grid, policy_k(:, iz), 'LineWidth', 2);
end
xlabel('Capital, k');
ylabel('Policy: k''(k,z)');
title('Policy Function (VFI)');
legend(arrayfun(@(z) sprintf('z = %.2f', z), z_vals, 'UniformOutput', false));
grid on;

%% (c) SIMULATION OF THE CAPITAL STOCK PATH
disp('--------------------------------------------------');
disp('Problem 2(c): Simulation over 10,000 periods');

T = 10000;           % Number of simulation periods
k_sim   = zeros(T,1);
z_sim   = zeros(T,1);   % To record actual z levels
z_ind_sim = zeros(T,1); % To record shock state indices

% Initialize simulation at the nonstochastic steady state and median shock.
[~, mid_ind] = min(abs(k_grid - k_ss));
k_sim(1)    = k_grid(mid_ind);
z_ind_sim(1)= ceil(N_z/2);
z_sim(1)    = z_vals(z_ind_sim(1));

% Precompute cumulative probabilities for the shock transition matrix.
cumP = cumsum(P,2);

rng(123);  % Set seed for reproducibility.
for t = 2:T
    % Determine next shock state using the Markov chain.
    r = rand;
    current_state = z_ind_sim(t-1);
    next_state = find(cumP(current_state,:) >= r, 1, 'first');
    z_ind_sim(t) = next_state;
    z_sim(t) = z_vals(next_state);
    % Find nearest index for current capital.
    [~, ind_k] = min(abs(k_grid - k_sim(t-1)));
    % Use the policy function corresponding to the current shock state.
    k_sim(t) = policy_k(ind_k, current_state);
end

% Plot the simulated capital path.
figure;
plot(1:T, k_sim, 'b', 'LineWidth',1.5);
xlabel('Time Period');
ylabel('Capital, k');
title('Simulated Capital Stock Path');
grid on;

% Plot a histogram of simulated capital.
figure;
histogram(k_sim, 50, 'Normalization','pdf','FaceColor',[0.7 0.3 0.3]);
xlabel('Capital, k');
ylabel('Density');
title('Histogram of Simulated Capital');
grid on;

