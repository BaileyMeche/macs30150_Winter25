%%%%%%%%%%%%% ----- FILE: pset_q2 -------%%%%%%%%%%
%%%%%% Problem 2: The Neoclassical Growth Model II

% maximize v(k) = \max_{c,k'} \left[ u(c) + \beta v(k') \right]
% subject to k' + c = F(k, 1    ) + (1 - \delta)k
% u(c) = \frac{c^{1-\sigma}}{1-\sigma}
%%% (Part a): Find the Euler equation.
% c_t^{-2}=\beta c_{t+1}^{-2}[\alpha k_{t+1}^{\,\alpha-1} + 1-\delta]

%%% (Part b) Find the steady state of the model.
% \bar{c}=\bar{k}^{\alpha}+(1-\delta)\,\bar{k}-\bar{k}

%%% (Part c): Find a linear approximation around the steady state of both the Euler equation and the feasibility constraint
% \begin{pmatrix}
% \hat{k}_{t+1}\\[6pt]
% \hat{c}_{t+1}
% \end{pmatrix}
% =
% \begin{pmatrix}
% \frac{1}{\beta} & -1 \\ \frac{c^{-\sigma} \alpha(\alpha-1)k^{\alpha-2} }{\sigma c^{-(1+\delta)}} & 1+ \frac{\beta c^{-\sigma} \alpha(\alpha-1)k^{\alpha-2}}{\sigma c^{-(1+\delta)}}
% \end{pmatrix}
% \begin{pmatrix}
% \hat{k}_t\\[4pt]
% \hat{c}_t
% \end{pmatrix}

% Parameters
beta = 0.95;
alpha = 0.36;
delta = 0.025;
sigma = 2;

% Find steady state k and c first
k_ss = ((1/beta - (1-delta))/(alpha))^(1/(alpha-1));
c_ss = k_ss^alpha - delta*k_ss;

% Calculate derivatives at steady state
f_prime = alpha * k_ss^(alpha-1);
f_double_prime = alpha * (alpha-1) * k_ss^(alpha-2);
u_double_prime = -sigma * c_ss^(-sigma-1);
u_prime = c_ss^(-sigma);

% Construct matrix A
A = [1/beta, -1;
     -(u_prime*f_double_prime/u_double_prime), 1 + (beta*u_prime*f_double_prime/u_double_prime)];

%%% (Part e)
% Find eigenvalues and eigenvectors
[V, D] = eig(A);

% Policy functions coefficients
V_inv = inv(V);
v21 = V_inv(2,1);
v22 = V_inv(2,2);
consumption_policy = v21/v22;
capital_policy = min(diag(D)); % lambda1

% Simulate transition dynamics (Part e)
T = 100;
k0 = 0.1 * k_ss;

k_path = zeros(T+1, 1);
c_path = zeros(T+1, 1);
k_path(1) = k0;

for t = 1:T
    k_dev = k_path(t) - k_ss;
    c_path(t) = c_ss + consumption_policy * k_dev;
    k_path(t+1) = k_ss + capital_policy * k_dev;
end
c_path(T+1) = c_ss + consumption_policy * (k_path(T+1) - k_ss);

% Figure for Part (e)
% Transition simulation using policy functions
figure;
subplot(2,1,1);
plot(0:T, k_path, 'r', 'DisplayName', 'Capital Path (Part e)');
title('Capital Path (Part e: Policy Function Simulation)');
ylabel('Capital');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T, c_path, 'r', 'DisplayName', 'Consumption Path (Part e)');
title('Consumption Path (Part e: Policy Function Simulation)');
ylabel('Consumption');
xlabel('Time');
legend show;
grid on;

% Run Dynare nonlinear simulation (Part f)
dynare PQuestion2F.mod noclearall;

% Extract Dynare results
dynare_k = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_c = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';

% Figure for Part (f)
% Compare MATLAB Part (e) with Dynare nonlinear results
figure;
subplot(2,1,1);
plot(0:T, k_path, 'r', 'DisplayName', 'MATLAB Capital Path (Part e)');
hold on;
plot(0:(length(dynare_k)-1), dynare_k, 'b--', 'DisplayName', 'Dynare Capital Path (Part f)');
title('Capital Path (Part f: Nonlinear Simulation Comparison)');
ylabel('Capital');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T, c_path, 'r', 'DisplayName', 'MATLAB Consumption Path (Part e)');
hold on;
plot(0:(length(dynare_c)-1), dynare_c, 'b--', 'DisplayName', 'Dynare Consumption Path (Part f)');
title('Consumption Path (Part f: Nonlinear Simulation Comparison)');
ylabel('Consumption');
xlabel('Time');
legend show;
grid on;

% Run Dynare linear simulation (Part g)
dynare PQuestion2G.mod noclearall;

% Extract Dynare results
dynare_k_linear = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_c_linear = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';

% Figure for Part (g)
% Compare MATLAB Part (e) with Dynare linear results
figure;
subplot(2,1,1);
plot(0:T, k_path, 'r', 'DisplayName', 'MATLAB Capital Path (Part e)');
hold on;
plot(0:(length(dynare_k_linear)-1), dynare_k_linear, 'b--', 'DisplayName', 'Dynare Capital Path (Part g)');
title('Capital Path (Part g: Linear Approximation Comparison)');
ylabel('Capital');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T, c_path, 'r', 'DisplayName', 'MATLAB Consumption Path (Part e)');
hold on;
plot(0:(length(dynare_c_linear)-1), dynare_c_linear, 'b--', 'DisplayName', 'Dynare Consumption Path (Part g)');
title('Consumption Path (Part g: Linear Approximation Comparison)');
ylabel('Consumption');
xlabel('Time');
legend show;
grid on;

fprintf('Policy functions:\n');
fprintf('c_t - c = %.4f(k_t - k)\n', consumption_policy);
fprintf('k_{t+1} - k = %.4f(k_t - k)\n', capital_policy);
