% Parameters
beta = 0.98;
sigma = 2;
y = 1; % Fixed endowment
tau = 0.0234;
gamma = 0.02;
theta_1 = 0.74;
theta_2 = 1.36;

% Transition matrix
% Transition prob starting on theta_1
p11 = 0.73;
p21 = 1 - p11;
% Transition probabilities starting in theta_2
p22 = 0.73;
p12 = 1 - p22;

Np = 1000; % Size of the grid
mp_min = y + tau; % CIA constraint (m' = y + tau)
mp_max = 2; % Maximum money level (optional)

% Create the grid of m
mgrid = linspace(mp_min, mp_max, Np)';

% Maximum number of iterations
Nk = 1000;

% Distribution parameters
% Finer grid for the distribution (twice as many points)
Dp = 2*Np-1;
dist_mgrid = linspace(mp_min, mp_max, Dp)';

Dk = 1000; % Maximum number of distribution iterations

% Create a transition matrix
QQ = [p11 p21; p12 p22];
nshocks = 2;

% Create the ergodic distribution
prob = QQ - eye(nshocks);
prob(:, nshocks) = ones(nshocks, 1);
a = zeros(nshocks, 1);
a(nshocks) = 1;
epr = linsolve(prob', a);


%% Case A Implementation
% Initialize distributions for both cases
F1_a = zeros(Dp, 1);
F2_a = zeros(Dp, 1);

% Create uniform distribution as initial guess
unif_dist = (dist_mgrid - dist_mgrid(1)) / (dist_mgrid(end) - dist_mgrid(1));

% Initial distributions
F1_a = epr(1) * unif_dist;
F2_a = epr(2) * unif_dist;

% Pre-compute the inverse policy functions for efficiency
% We need to find g^-1(m') for each m' in our grid
g1_unique = unique(g_pol_1, 'stable');
g2_unique = unique(g_pol_2, 'stable');

% Find the corresponding m values for each unique g value
m1_values = zeros(size(g1_unique));
m2_values = zeros(size(g2_unique));

for i = 1:length(g1_unique)
    idx = find(g_pol_1 == g1_unique(i), 1, 'first');
    m1_values(i) = mgrid(idx);
end

for i = 1:length(g2_unique)
    idx = find(g_pol_2 == g2_unique(i), 1, 'first');
    m2_values(i) = mgrid(idx);
end

% Case A iterations
disp('Starting Case A iterations...');
for iter = 1:Dk
    F1_a_new = zeros(Dp, 1);
    F2_a_new = zeros(Dp, 1);
    
    for j = 1:Dp
        mp = dist_mgrid(j);
        
        % Find g^-1(mp) for each shock
        % For low shock (θ₁)
        if mp <= g_pol_1(1)
            F1_a_new(j) = 0;
        elseif mp >= g_pol_1(end)
            F1_a_new(j) = F1_a(end);
        else
            % Find g^-1(mp, θ₁)
            m_value = interp1(g1_unique, m1_values, mp, 'linear');
            % Find Ψₛ(g^-1(mp, θ₁), θ₁)
            F1_a_new(j) = interp1(dist_mgrid, F1_a, m_value, 'linear');
        end
        
        % For high shock (θ₂)
        if mp <= g_pol_2(1)
            F2_a_new(j) = 0;
        elseif mp >= g_pol_2(end)
            F2_a_new(j) = F2_a(end);
        else
            % Find g^-1(mp, θ₂)
            m_value = interp1(g2_unique, m2_values, mp, 'linear');
            % Find Ψₛ(g^-1(mp, θ₂), θ₂)
            F2_a_new(j) = interp1(dist_mgrid, F2_a, m_value, 'linear');
        end
    end
    
    % Check convergence
    diff1 = max(abs(F1_a_new - F1_a));
    diff2 = max(abs(F2_a_new - F2_a));
    
    % Update distributions
    F1_a = F1_a_new;
    F2_a = F2_a_new;
    
    if mod(iter, 100) == 0
        disp(['Iteration ' num2str(iter) ', Max diff: ' num2str(max(diff1, diff2))]);
    end
    
    if max(diff1, diff2) < 1e-8
        disp(['Case A converged after ' num2str(iter) ' iterations']);
        break;
    end
    
    if iter == Dk
        disp('Case A reached maximum iterations without converging');
    end
end

%% Case B Implementation
%% Case B Implementation (Improved)
% Initialize distributions for Case B
F1_b = zeros(Dp, 1);
F2_b = zeros(Dp, 1);

% Initial distributions (same as Case A)
F1_b = epr(1) * unif_dist;
F2_b = epr(2) * unif_dist;

% Case B iterations
disp('Starting Case B iterations...');
for iter = 1:Dk
    F1_b_new = zeros(Dp, 1);
    F2_b_new = zeros(Dp, 1);
    
    for j = 1:Dp
        mp = dist_mgrid(j);
        
        % The recursive equation for Case B is:
        % Ψₛ₊₁(m',θⱼ) = ∑ Q(θ,θⱼ)Ψₛ(g⁻¹(m',θ),θ)
        
        % Low shock (j=1)
        
        % Contribution from transition θ₁ → θ₁
        if mp <= g_pol_1(1)
            val1_1 = 0;
        elseif mp >= g_pol_1(end)
            val1_1 = p11 * F1_b(end);  % If mp is beyond the range, use the last value
        else
            % Find g^-1(mp, θ₁)
            m_value = interp1(g1_unique, m1_values, mp, 'linear');
            % Find Ψₛ(g^-1(mp, θ₁), θ₁) and multiply by p11
            val1_1 = p11 * interp1(dist_mgrid, F1_b, m_value, 'linear');
        end
        
        % Contribution from transition θ₂ → θ₁
        if mp <= g_pol_2(1)
            val1_2 = 0;
        elseif mp >= g_pol_2(end)
            val1_2 = p12 * F2_b(end);  % If mp is beyond the range, use the last value
        else
            % Find g^-1(mp, θ₂)
            m_value = interp1(g2_unique, m2_values, mp, 'linear');
            % Find Ψₛ(g^-1(mp, θ₂), θ₂) and multiply by p12
            val1_2 = p12 * interp1(dist_mgrid, F2_b, m_value, 'linear');
        end
        
        val1 = val1_1 + val1_2;
        F1_b_new(j) = val1;
        
        % High shock (j=2)
        
        % Contribution from transition θ₁ → θ₂
        if mp <= g_pol_1(1)
            val2_1 = 0;
        elseif mp >= g_pol_1(end)
            val2_1 = p21 * F1_b(end);  % If mp is beyond the range, use the last value
        else
            % Find g^-1(mp, θ₁)
            m_value = interp1(g1_unique, m1_values, mp, 'linear');
            % Find Ψₛ(g^-1(mp, θ₁), θ₁) and multiply by p21
            val2_1 = p21 * interp1(dist_mgrid, F1_b, m_value, 'linear');
        end
        
        % Contribution from transition θ₂ → θ₂
        if mp <= g_pol_2(1)
            val2_2 = 0;
        elseif mp >= g_pol_2(end)
            val2_2 = p22 * F2_b(end);  % If mp is beyond the range, use the last value
        else
            % Find g^-1(mp, θ₂)
            m_value = interp1(g2_unique, m2_values, mp, 'linear');
            % Find Ψₛ(g^-1(mp, θ₂), θ₂) and multiply by p22
            val2_2 = p22 * interp1(dist_mgrid, F2_b, m_value, 'linear');
        end
        
        val2 = val2_1 + val2_2;
        F2_b_new(j) = val2;
    end
    
    % Check convergence
    diff1 = max(abs(F1_b_new - F1_b));
    diff2 = max(abs(F2_b_new - F2_b));
    
    % Update distributions
    F1_b = F1_b_new;
    F2_b = F2_b_new;
    
    if mod(iter, 100) == 0
        disp(['Iteration ' num2str(iter) ', Max diff: ' num2str(max(diff1, diff2))]);
    end
    
    if max(diff1, diff2) < 1e-8
        disp(['Case B converged after ' num2str(iter) ' iterations']);
        break;
    end
    
    if iter == Dk
        disp('Case B reached maximum iterations without converging');
    end
end

%% Plotting
% Plot results for Case A
figure;
subplot(2,1,1);
plot(mgrid, g_pol_1, 'b-', 'LineWidth', 1.5); hold on;
plot(mgrid, mgrid, 'k--', 'LineWidth', 1); % 45-degree line
grid on;
xlabel('Money (m)');
ylabel('Money next period (m'')');
legend('Low shock (θ₁)', '45-degree line');
title('Policy Function for Low Shock (θ₁)');

subplot(2,1,2);
plot(dist_mgrid, F1_a, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Money (m)');
ylabel('Ψ(m,θ₁)');
title('Stationary Distribution for Low Shock (θ₁) - Case A');

figure;
subplot(2,1,1);
plot(mgrid, g_pol_2, 'r-', 'LineWidth', 1.5); hold on;
plot(mgrid, mgrid, 'k--', 'LineWidth', 1); % 45-degree line
grid on;
xlabel('Money (m)');
ylabel('Money next period (m'')');
legend('High shock (θ₂)', '45-degree line');
title('Policy Function for High Shock (θ₂)');

subplot(2,1,2);
plot(dist_mgrid, F2_a, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Money (m)');
ylabel('Ψ(m,θ₂)');
title('Stationary Distribution for High Shock (θ₂) - Case A');


% Plot results for Case B
figure;
subplot(2,1,1);
plot(mgrid, g_pol_1, 'b-', 'LineWidth', 1.5); hold on;
plot(mgrid, mgrid, 'k--', 'LineWidth', 1); % 45-degree line
grid on;
xlabel('Money (m)');
ylabel('Money next period (m'')');
legend('Low shock (θ₁)', '45-degree line');
title('Policy Function for Low Shock (θ₁)');

subplot(2,1,2);
plot(dist_mgrid, F1_b, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Money (m)');
ylabel('Ψ(m,θ₁)');
title('Stationary Distribution for Low Shock (θ₁) - Case B');

figure;
subplot(2,1,1);
plot(mgrid, g_pol_2, 'r-', 'LineWidth', 1.5); hold on;
plot(mgrid, mgrid, 'k--', 'LineWidth', 1); % 45-degree line
grid on;
xlabel('Money (m)');
ylabel('Money next period (m'')');
legend('High shock (θ₂)', '45-degree line');
title('Policy Function for High Shock (θ₂)');

subplot(2,1,2);
plot(dist_mgrid, F2_b, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Money (m)');
ylabel('Ψ(m,θ₂)');
title('Stationary Distribution for High Shock (θ₂) - Case B');


% Plot distributions
figure;
subplot(2,1,1);
plot(dist_mgrid, F1_a, 'b-', 'LineWidth', 1.5); hold on;
plot(dist_mgrid, F1_b, 'b--', 'LineWidth', 1.5);
grid on;
xlabel('Money (m)');
ylabel('Ψ(m,θ₁)');
legend('Case A', 'Case B');
title('Comparison of Stationary Distributions for Low Shock (θ₁)');

subplot(2,1,2);
plot(dist_mgrid, F2_a, 'r-', 'LineWidth', 1.5); hold on;
plot(dist_mgrid, F2_b, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Money (m)');
ylabel('Ψ(m,θ₂)');
legend('Case A', 'Case B');
title('Comparison of Stationary Distributions for High Shock (θ₂)');