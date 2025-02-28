%% ---------------------- Part (a): VFI with Variable Labor (σ = δ = 1) ----------------------
clear all; close all; clc;

% The planner's problem:
%
%   V(k) = max_{c,k',ℓ} [ u(c,ℓ) + β V(k') ]
%
% subject to:
%   
%   c + k' = k^α ℓ^(1-α) + (1-δ)k
%
% Here, we assume σ = 1 and δ = 1.
%
% Utility function (log form for σ = 1):
%
%   u(c,ℓ) = γ log(c) + (1-γ) log(1-ℓ)
%
% We solve this problem using Value Function Iteration (VFI).

% ---------------------- Parameters ----------------------
alpha = 0.36;
beta  = 0.95;
delta = 1;       % as given in part (a)
sigma = 1;       % logarithmic utility limit
gamma = 0.5;     % weight on consumption (arbitrary)
lbar  = 1/3;     % steady-state labor assumption

% ---------------------- Grid Setup ----------------------
Np    = 500; 
k_min = 0.1;
k_max = 20;
kgrid = linspace(k_min, k_max, Np)';

% Nl    = 100;  
% l_min = 0.01;
% l_max = 0.99;
% lgrid = linspace(l_min, l_max, Nl)';
% 
% % Initial guess: Set to zero instead of trying to approximate steady-state utility.
% V = zeros(Np,1);
% 
% % Preallocate arrays
% policy_k  = zeros(Np,1);
% policy_l  = zeros(Np,1);
% policy_c  = zeros(Np,1);
% V_new     = zeros(Np,1);
% 
% % ---------------------- VFI Parameters ----------------------
% tol      = 1e-5;   % Convergence tolerance
% maxiter  = 500;    % Reduce max iterations
% 
% disp('Running Value Function Iteration for Part (a)...');
% 
% for iter = 1:maxiter
%     V_old = V;
% 
%     for ik = 1:Np
%         k = kgrid(ik);
%         max_val = -Inf;
%         best_kp = k_min; % Default to smallest possible k'
%         best_l = l_min;
%         best_c = 1e-6; % Ensure nonzero consumption
% 
%         % Compute next-period wealth for all choices of labor
%         for il = 1:Nl
%             ell = lgrid(il);
%             output = k^alpha * ell^(1-alpha);
% 
%             % Feasible consumption range (avoid negative values)
%             kprime_max = min(k_max, output); % k' cannot exceed available wealth
%             kprime_grid = linspace(k_min, kprime_max, Np);
% 
%             % Compute consumption
%             c_vals = output - kprime_grid;
%             valid = c_vals > 0; % Ensure c > 0
% 
%             if any(valid)
%                 util_vals = gamma * log(c_vals(valid)) + (1 - gamma) * log(1 - ell);
% 
%                 % Interpolate V(k')
%                 V_future = interp1(kgrid, V_old, kprime_grid(valid), 'linear', 'extrap');
% 
%                 % Compute Bellman RHS
%                 val = util_vals + beta * V_future;
% 
%                 % Find the max value
%                 [max_v, idx] = max(val);
% 
%                 % Assign the single best k', ℓ, c
%                 if max_v > max_val
%                     max_val = max_v;
%                     best_kp = kprime_grid(valid); % Select only the best index
%                     best_l = ell;
%                     best_c = c_vals(valid);
% 
%                     % Extract the best scalar values
%                     best_kp = best_kp(idx);
%                     best_c = best_c(idx);
%                 end
%             end
%         end
% 
%         % Store best choices
%         V_new(ik) = max_val;
%         policy_k(ik) = best_kp;
%         policy_l(ik) = best_l;
%         policy_c(ik) = best_c;
%     end
% 
%     % Check convergence
%     diff = max(abs(V_new - V_old));
%     V = V_new;
% 
%     if mod(iter,10) == 0
%         fprintf('Iteration %d, diff = %.8f\n', iter, diff);
%     end
% 
%     if diff < tol
%         fprintf('Convergence achieved in %d iterations (Part a).\n', iter);
%         break;
%     end
% end
% 
% % ---------------------- Plot Policy Functions ----------------------
% figure;
% subplot(3,1,1);
% plot(kgrid, policy_c, 'b', 'LineWidth', 1.5);
% title('Part (a): Consumption Policy Function');
% xlabel('Capital, k'); ylabel('Consumption, c'); grid on;
% 
% subplot(3,1,2);
% plot(kgrid, policy_k, 'r', 'LineWidth', 1.5);
% title('Part (a): Next–Period Capital Policy Function');
% xlabel('Capital, k'); ylabel('k'''); grid on;
% 
% subplot(3,1,3);
% plot(kgrid, policy_l, 'k', 'LineWidth', 1.5);
% title('Part (a): Labor Policy Function');
% xlabel('Capital, k'); ylabel('Labor, ℓ'); grid on;
% 
% % ---------------------- Plot Final Value Function ----------------------
% figure;
% plot(kgrid, V, 'LineWidth', 1.5);
% title('Part (a): Final Value Function V(k)');
% xlabel('Capital, k'); ylabel('V(k)'); grid on;
% 
% fprintf('VFI Completed.\n');


%% ---------------------- Part (c): Steady–State Calibration ---------------------- 
% Given:
%   α = 0.36, δ = 0.03, β = 0.98, σ = 2, and l̄ = 1/3.
% We solve for:
%   k̄, c̄, and γ using steady-state conditions.

alpha = 0.36; 
delta = 0.03; 
beta  = 0.98; 
sigma = 2;
lbar  = 1/3;

% Solve Euler equation for k̄
kbar = fzero(@(k) beta*(alpha*k^(alpha-1)*lbar^(1-alpha)+1-delta)-1, 10);

% Compute steady-state consumption
cbar = kbar^alpha * lbar^(1-alpha) - delta*kbar;

% Calibrate γ from intratemporal condition
rhs = (1-alpha)*kbar^alpha * lbar^(-alpha);
gamma_cal = 1 / (1 + rhs*(1-lbar)/cbar);

fprintf('\nPart (c): Steady State Results\n');
fprintf('   k̄ = %.4f\n', kbar);
fprintf('   c̄ = %.4f\n', cbar);
fprintf('   Calibrated γ = %.4f\n', gamma_cal);


%% ---------------------- Part (d): Orthogonal Collocation ----------------------
clear all; close all; clc;

% Approximate the policy functions for next–period capital and labor:
%
%   k' ≡ ĝ(k; a) = Σ_{i=1}^{n} a_i φ_i(k)
%   ℓ ≡ ĥ(k; b) = Σ_{i=1}^{n} b_i ψ_i(k)
%
% with basis functions:
%
%   φ_i(k) = k·T_{i-1}(x),   ψ_i(k) = T_{i-1}(x),
%
% where T_{i-1}(x) are Chebyshev polynomials, and 
% x = (2*k - (k_min+k_max))/(k_max-k_min) transforms k into [-1, 1].
%
% We solve for the coefficient vectors a and b such that the Euler 
% and intratemporal conditions (expressed via the residual operator) vanish.

% ---------------------- Parameters ----------------------
alpha     = 0.36; 
delta     = 0.03;
beta      = 0.98;
sigma     = 2;
lbar      = 1/3;
gamma_cal = 0.48; 

% Compute steady-state capital from the Euler condition:
kbar = fzero(@(k) beta*(alpha*k^(alpha-1)*lbar^(1-alpha)+1-delta)-1, 10);
cbar = kbar^alpha * lbar^(1-alpha) - delta*kbar;

% Collocation Domain ----------------------
n       = 8;           % number of basis functions
k_min   = 0.3 * kbar;  % lower bound for collocation domain
k_max   = 1.2 * kbar;  % upper bound for collocation domain
j       = (1:n)';
x_nodes = cos((2*j-1)*pi/(2*n));         % Chebyshev nodes in [-1,1]
k_nodes = ((k_max - k_min)/2)*x_nodes + (k_max+k_min)/2;  % transform to k-space

% Basis Functions ----------------------
% Chebyshev polynomial basis functions.
T   = @(i,x) cos((i-1).*acos(x));  % T_{i-1}(x)
phi = @(i,k) k .* T(i, (2*k - (k_min+k_max))/(k_max-k_min));
psi = @(i,k) T(i, (2*k - (k_min+k_max))/(k_max-k_min));

% ---------------------- Initial Guess ----------------------
% A good initial guess is to set a1 ≈ 1 (so that g(k) ≈ k) and b1 ≈ lbar.
X0 = [1; zeros(n-1,1); lbar; zeros(n-1,1)];  % X = [a; b]

% ---------------------- Solve Collocation System ----------------------
options = optimoptions('fsolve','Display','iter','TolFun',1e-8);
[X_sol, ~, ~, ~] = fsolve(@(X) collocation_resid(X, k_nodes, alpha, beta, delta, sigma, gamma_cal, k_min, k_max), X0, options);

% Extract coefficients
a = X_sol(1:n);
b = X_sol(n+1:2*n);

% ---------------------- Evaluate Approximated Policy Functions ----------------------
Nfine = 100;
k_fine = linspace(k_min, k_max, Nfine)';
g_k  = zeros(size(k_fine));
h_k  = zeros(size(k_fine));
for i = 1:Nfine
    k_val = k_fine(i);
    temp_phi = zeros(n,1);
    temp_psi = zeros(n,1);
    for j = 1:n
        temp_phi(j) = phi(j, k_val);
        temp_psi(j) = psi(j, k_val);
    end
    g_k(i) = a' * temp_phi;
    h_k(i) = b' * temp_psi;
end
% Consumption policy: c = k^α ℓ^(1-α) + (1-δ)*k - g(k)
c_policy = k_fine.^alpha .* h_k.^(1-alpha) + (1-delta)*k_fine - g_k;

% ---------------------- Plot Policy Functions ----------------------
figure;
subplot(3,1,1);
plot(k_fine, g_k, 'b','LineWidth',1.5);
title('Part (d): Approximated Policy for Next–Period Capital');
xlabel('Capital, k'); ylabel('k'''); grid on;

subplot(3,1,2);
plot(k_fine, h_k, 'r','LineWidth',1.5);
title('Part (d): Approximated Policy for Labor');
xlabel('Capital, k'); ylabel('Labor, ℓ'); grid on;

subplot(3,1,3);
plot(k_fine, c_policy, 'k','LineWidth',1.5);
title('Part (d): Approximated Policy for Consumption');
xlabel('Capital, k'); ylabel('Consumption, c'); grid on;

% ---------------------- Collocation Residual Function ----------------------
function F = collocation_resid(X, k_nodes, alpha, beta, delta, sigma, gamma_cal, k_min, k_max)
% X is a vector of length 2*n: the first n entries are the coefficients a for 
% the capital policy, and the next n entries are the coefficients b for the labor policy.
n = length(X) / 2;
a = X(1:n);
b = X(n+1:2*n);

F = zeros(2*n, 1);

% Define basis functions in this function as well, using the provided k_min and k_max.
T   = @(i,x) cos((i-1).*acos(x));
phi = @(i,k) k .* T(i, (2*k - (k_min+k_max))/(k_max-k_min));
psi = @(i,k) T(i, (2*k - (k_min+k_max))/(k_max-k_min));

for j = 1:n
    k = k_nodes(j);
    % Compute the basis functions at the current state k.
    temp_phi = zeros(n,1);
    temp_psi = zeros(n,1);
    for i = 1:n
        temp_phi(i) = phi(i, k);
        temp_psi(i) = psi(i, k);
    end
    kprime = a' * temp_phi;   % g(k)
    ell    = b' * temp_psi;    % h(k)
    
    % Current consumption: c = k^α * ell^(1-α) + (1-δ)*k - g(k)
    c = k^alpha * ell^(1-alpha) + (1-delta)*k - kprime;
    
    % Evaluate future policies at kprime.
    temp_phi_future = zeros(n,1);
    temp_psi_future = zeros(n,1);
    for i = 1:n
        temp_phi_future(i) = phi(i, kprime);
        temp_psi_future(i) = psi(i, kprime);
    end
    kprime_future = a' * temp_phi_future;   % g(k')
    ell_future    = b' * temp_psi_future;     % h(k')
    
    % Future consumption: c' = (kprime)^α * (ell_future)^(1-α) + (1-δ)*kprime - g(k')
    c_future = kprime^alpha * ell_future^(1-alpha) + (1-delta)*kprime - kprime_future;
    
    % Compute marginal utilities for current consumption:
    % u(c,ℓ) = [ c^(gamma_cal) (1-ℓ)^(1-gamma_cal) ]^(1-sigma) / (1-sigma)
    % Thus, u_c = [c^(gamma_cal)*(1-ℓ)^(1-gamma_cal)]^(-sigma)*gamma_cal*c^(gamma_cal-1)*(1-ℓ)^(1-gamma_cal)
    %      u_ell = -[c^(gamma_cal)*(1-ℓ)^(1-gamma_cal)]^(-sigma)*(1-gamma_cal)*c^(gamma_cal)*(1-ℓ)^(-gamma_cal-1)
    uc = (c^(gamma_cal) * (1-ell)^(1-gamma_cal))^(-sigma) * gamma_cal * c^(gamma_cal-1) * (1-ell)^(1-gamma_cal);
    ul = - (c^(gamma_cal) * (1-ell)^(1-gamma_cal))^(-sigma) * (1-gamma_cal) * c^(gamma_cal) * (1-ell)^(-gamma_cal-1);
    
    % Marginal utility for future consumption:
    uc_future = (c_future^(gamma_cal) * (1-ell_future)^(1-gamma_cal))^(-sigma) * gamma_cal * c_future^(gamma_cal-1) * (1-ell_future)^(1-gamma_cal);
    
    % Euler equation residual:
    %   u_c(c,ell) - beta * u_c(c',ell')*[alpha * kprime^(alpha-1)*ell_future^(1-alpha) + 1-delta] = 0
    RE = uc - beta * uc_future * (alpha * kprime^(alpha-1) * ell_future^(1-alpha) + 1-delta);
    
    % Intratemporal (labor) residual:
    %   u_ell(c,ell) + u_c(c,ell)*(1-alpha)*k^alpha*ell^(-alpha) = 0
    RL = ul + uc*(1-alpha)*k^alpha*ell^(-alpha);
    
    F(j)   = RE;
    F(j+n) = RL;
end
end
