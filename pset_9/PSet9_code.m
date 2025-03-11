% PSet 9

%% Policy functions VFI

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
mp_max = 2;  % Maximum money level (optional)

% Create the grid of a
mgrid = linspace(mp_min, mp_max, Np)';

% Maximum number of iterations
Nk = 1000;

% consumption has a direct correspondence with the other variables 
% endowment is fixed
c0 = (mgrid ./ (1 + gamma)) + y + tau - mgrid';

% If any of the above are negative, replace with zero
c0(c0 < 0) = 0;

% Define the matrices of returns
ut_c = u_t(sigma, c0);

% Definition of placeholders for Value Functions
g_1 = zeros(Np, Nk);
g_2 = zeros(Np, Nk);

% Initial guess for value function (V(m',theta'))
g_1(:,1) = real(u_t(sigma, mgrid)); 
g_2(:,1) = real(u_t(sigma, mgrid));

% Iterations to fill vf with the maximum
for i = 1:Nk-1
    % Compute value functions
    % VF = theta * u(c) + beta * sum[VF(m',theta')*Q(theta, theta')]
    % Where VF(m',theta') is g_1 or g_2, respectively
    VF_1 = theta_1 * ut_c + (beta * p11 * g_1(:,i)' + beta * p12 * g_2(:,i)');
    VF_2 = theta_2 * ut_c + (beta * p21 * g_1(:,i)' + beta * p22 * g_2(:,i)');

    % Solve for the next period value function
    [g_1(:, i+1), mind_1] = max(VF_1'); 
    [g_2(:, i+1), mind_2] = max(VF_2');

    % Check for convergence
    vdiff = [g_1(:, i+1); g_2(:, i+1)] - [g_1(:, i); g_2(:, i)];
    if max(abs(vdiff)) <= 0.001 * (1 - beta)
        disp("Convergence achieved");
        conv = i;
        break;
    end
end

% Computing the policy functions

% Define the assets placeholders
mpol_1 = zeros(Np, 1);
mpol_2 = zeros(Np, 1);

for i = 1:Np
    mpol_1(i) = mgrid(mind_1(i));
    mpol_2(i) = mgrid(mind_2(i));
end

% Define the consumption 
cpol_1 = (mgrid ./ (1 + gamma)) + y + tau - mpol_1;
cpol_2 = (mgrid ./ (1 + gamma)) + y + tau - mpol_2;

% Plotting value functions
figure;
plot(mgrid, g_1(:, conv), 'b'); hold on;
plot(mgrid, g_2(:, conv), 'r');
grid;
xlabel('Money (m)');
ylabel('Value Function');
title('Value Functions');
legend({'Low theta', 'High theta'}, 'location', 'southeast');
hold off;

% Plotting policy functions for money
figure; 
plot(mgrid, mpol_1); hold on;
plot(mgrid, mpol_2); hold on;
plot(mgrid, mgrid)
grid;
xlabel('a');
ylabel("Assets next period (a')");
title("Policy functions for a'");
legend({'Low theta', 'High theta', '45'}, 'location', 'southeast');

% Plotting policy functions for consumption
figure; 
plot(mgrid, cpol_1); hold on;
plot(mgrid, cpol_2);
grid;
xlabel('a');
ylabel('Consumption');
title('Policy functions for consumption');
legend({'Low theta', 'High theta'}, 'location', 'southeast');

% Define the utility function
function result = u_t(sigma, c)
    % Ensure no log of negative or zero values
    c(c <= 0) = 1e-6; 
    if abs(sigma - 1) == 0
        result = log(c); % Log utility for sigma = 1
    else
        result = (c.^(1 - sigma)) / (1 - sigma); % CRRA utility
    end
end


%% Policy functions EGM

clear; clc; close all;

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
mp_max = 2;  % Maximum money level (optional)

% Create the grid of m
mgrid = linspace(mp_min, mp_max, Np)';

% Maximum number of iterations
Nk = 1000;

% Initial guesses of m
g0_1 = mgrid;
g0_2 = mgrid;

% Preallocate memory for future iterations
g_1 = zeros(Np, Nk);
g_2 = zeros(Np, Nk);

% First column is the initial guess
g_1(:,1) = g0_1;
g_2(:,1) = g0_2;

for j = 1:Nk-1
    uprime_1 = theta_1 * ((mgrid ./ (1 + gamma)) + y + tau - g_1(:,j)) .^ (-sigma);
    uprime_2 = theta_2 * ((mgrid ./ (1 + gamma)) + y + tau - g_2(:,j)) .^ (-sigma);
    
    % Compute expectations under the Markov transition matrix
    EV_1 = ((beta/(1+gamma)) * (p11 * uprime_1' + p12 * uprime_2'))';
    EV_2 = ((beta/(1+gamma)) * (p21 * uprime_1' + p22 * uprime_2'))';
    
    % Compute mstar
    ms_1 = ((EV_1/theta_1) .^ (-1/sigma) + mgrid - y - tau) * (1+gamma);
    ms_2 = ((EV_2/theta_2) .^ (-1/sigma) + mgrid - y - tau) * (1+gamma); 
    
    for i=1:Np
        if mgrid(i) <= ms_1(1)
            g_1(i,j+1) = mgrid(1);
        elseif mgrid(i) >= ms_1(Np)
            g_1(i, j+1) = mgrid(Np);
        else
            g_1(i, j+1) = interp1(ms_1, mgrid, mgrid(i), 'linear', 'extrap');
        end
    end

    for i=1:Np
        if mgrid(i) <= ms_2(1)
            g_2(i,j+1) = mgrid(1);
        elseif mgrid(i) >= ms_2(Np)
            g_2(i, j+1) = mgrid(Np);
        else
            g_2(i, j+1) = interp1(ms_2, mgrid, mgrid(i), 'linear', 'extrap');
        end
    end

    % Check convergence
    gdiff = max(abs([g_1(:,j+1); g_2(:,j+1)] - [g_1(:,j); g_2(:,j)]));
    if gdiff <= 1e-10
        disp("Convergence achieved after " + num2str(j) + " iterations.");
        g_pol_1 = g_1(:, j+1);
        g_pol_2 = g_2(:, j+1);
        break;
    end
end

% If convergence isn't achieved, store final iteration values
if gdiff > 1e-10
    disp("Maximum iterations reached without convergence.");
    g_pol_1 = g_1(:, end);
    g_pol_2 = g_2(:, end);
end


% Compute policy functions for consumption
c_pol_1 = (mgrid / (1 + gamma)) + y + tau - g_pol_1; 
c_pol_2 = (mgrid / (1 + gamma)) + y + tau - g_pol_2; 


% Plotting money policy function
figure;
plot(mgrid, g_pol_1, 'r', 'LineWidth', 2); hold on;
plot(mgrid, g_pol_2, 'b', 'LineWidth', 2);
plot(mgrid, mgrid, 'k--', 'LineWidth', 1); % 45-degree line
xlabel('Assets');
ylabel('Assets t+1');
legend('Low theta', 'High theta', '45-degree line');
title('Policy function: Assets');

% Plotting consumption policy function
figure;
plot(mgrid, c_pol_1, 'r', 'LineWidth', 2); hold on;
plot(mgrid, c_pol_2, 'b', 'LineWidth', 2);
xlabel('Assets');
ylabel('Consumption');
legend('Low theta', 'High theta');
title('Policy function: Consumption');


%% Computing distributions a)

% Parameters are the same for the policy functions

% Distribution parameters
% New grid for the distribution
Dp = 2*Np-1;
dist_mgrid = linspace(mp_min, mp_max, Dp)';

Dk = 1000; % Maximum number of distribution iterations

% Other options: normal dist
unif_dist = (dist_mgrid - dist_mgrid(1)) / (dist_mgrid(Dp) - dist_mgrid(1));


% Create a transition matrix with the prob
QQ = [p11 p21; p12 p22];
nshocks = 2;

% Create the ergodic distribution (cn 7)
% The ergodic distribution (or stationary distribution) of a Markov chain 
% represents the long-run probability of being in each state when the system 
% reaches equilibrium. It tells us how often, on average, we would expect the 
% system to be in each state after a long enough period.

% Linearization of the system of equations (eye is the identity matrix of
% size nshocks)
% pi1 = 0.73 * pi1 + 0.27 * pi2
% pi2 = 0.27 * pi1 + 0.73 * pi2

prob = QQ - eye(nshocks);
% [(0.73 - 1)  0.27]
% [0.27  (0.73 - 1)]

% but this has infinite solutions, hence we first enforce pi1 + pi2 = 1
% by enforcing the last column to be 1 (normalization condition)
% This is the same as normalizing the eigenvector to sum to 1 in tauchen's
% [-0.27  1]
% [ 0.27  1]

prob(:, nshocks) = ones(nshocks, 1);

% On the other side of the equality, we would have a = [0 1] to solve
% for the ergodic distribution over states that remain unchanged over the
% long run and sum up to 1

a = zeros(nshocks, 1);
a(nshocks) = 1;

% Ergodic distribution
epr = linsolve(prob', a);

% Placeholders for my distributions 
% initial_dist = (psi_0(m1) + [psi_0(mM) - psi_0(m1)] * unif_dist) * epr
% but since psi_0(mM)=1 and psi_0(m1)=0
F1 = zeros(Dp, Dk);
F1(:,1) = epr(1) * unif_dist; %epr(1) is for the ergodic dist for theta_1

F2 = zeros(Dp, Dk);
F2(:,1) = epr(2) * unif_dist; %epr(2) is for the ergodic dist for theta_2


% Iterations for convergence
for i = 1:Dk-1
    for j=1:Dp

        % For low shock

        if dist_mgrid(j) < g_pol_1(1)
            val_1 = 0;
        elseif dist_mgrid(j) > g_pol_1(Np)
            val_1 = epr(1);
        else
            % This is to obtain the unique values of the policy
            % function and extract the index for each one of these values.
            [g1r, ind_1] = unique(g_pol_1, "last");

            % For the first interpolation I want an interpolated value for dist_mgrid(j) 
            % using the function values (gr1) of mgrid(ind_1).
            % For the second interpolation I want an interpolated value of
            % the distribution function F1(:,i) at that query point (new dgrid).

            % The reason g1r is X and mgrid(ind_1) is Y is because I want
            % the value in mgrid (inverse policy function used)

            % Otherwise dist_mgrid is X and F1(:,i) is Y because I
            % want the value of the distribution such that I can use it
            % as my new value of the function in i+1.

            val_1 = interp1(dist_mgrid, F1(:,i), interp1(g1r, mgrid(ind_1), dist_mgrid(j), "linear"), "linear");
               
        end

        % For high shock
        
        if dist_mgrid(j) <= g_pol_2(1)
            val_2 = 0;
        elseif dist_mgrid(j) >= g_pol_2(Np)
            val_2 = epr(2);
        else
            [g2r, ind_2] = unique(g_pol_2, "last");

            val_2 = interp1(dist_mgrid, F2(:,i), interp1(g2r, mgrid(ind_2), dist_mgrid(j), "linear"), "linear");
        end
        
        % New distribution
        F1(j, i+1) = val_1;
        F2(j, i+1) = val_2;
        
    end

    diff1 = max(abs(F1(:,i+1)-F1(:,i)));
    diff2 = max(abs(F2(:,i+1)-F2(:,i)));

    F1_dist = F1(:,i+1);
    F2_dist = F2(:,i+1);

    if diff1 <= 1e-10
        disp("convergence achieved of F1");
        disp(i);
        break
    end

    if diff2 <= 1e-10
        disp("convergence achieved of F2");
        disp(i);
        break
    end
end




% Plot distribution functions when they converge
figure; 
subplot(2,1,1);
plot(mgrid, g_pol_1); hold on;
plot(mgrid, mgrid); % 45-degree line
grid on;
xlabel('Assets');
ylabel('Assets t+1');
legend('Low theta', '45-degree line');
title('Policy function: Assets');

subplot(2,1,2);
plot(dist_mgrid, F1(:,1)); hold on;
plot(dist_mgrid, F1(:,2)); hold on;
plot(dist_mgrid, F1(:,3)); hold on;
plot(dist_mgrid, F1(:,4)); hold on;
plot(dist_mgrid, F1_dist);
grid on;
xlabel('dgrid');
ylabel('Psi(m, theta_1)');
title('Distribution functions for low shock');
legend({'Psi_0', 'Psi_1', 'Psi_2', 'Psi_3', 'Last distribution'}, 'location', 'southeast');

figure; 
subplot(2,1,1); 
plot(mgrid, g_pol_2); hold on;
plot(mgrid, mgrid); % 45-degree line
grid on;
xlabel('Assets');
ylabel('Assets t+1');
legend('High theta', '45-degree line');
title('Policy function: Assets');

subplot(2,1,2);
plot(dist_mgrid, F2(:,1)); hold on;
plot(dist_mgrid, F2(:,2)); hold on;
plot(dist_mgrid, F2(:,3)); hold on;
plot(dist_mgrid, F2(:,4)); hold on;
plot(dist_mgrid, F2_dist);
grid on;
xlabel('dgrid');
ylabel('Psi(m, theta_2)');
title('Distribution functions for high shock');
legend({'Psi_0', 'Psi_1', 'Psi_2', 'Psi_3', 'Last distribution'}, 'location', 'southeast');


%% Computing distributions b)

% Parameters are the same for the policy functions

% Distribution parameters
% New grid for the distribution
Dp = 2*Np-1;
dist_mgrid = linspace(mp_min, mp_max, Dp)';

Dk = 1000; % Maximum number of distribution iterations

% Other options: normal dist
unif_dist = (dist_mgrid - dist_mgrid(1)) / (dist_mgrid(Dp) - dist_mgrid(1));


% Create a transition matrix with the prob
QQ = [p11 p21; p12 p22];
nshocks = 2;

% Create the ergodic distribution (cn 7)
% The ergodic distribution (or stationary distribution) of a Markov chain 
% represents the long-run probability of being in each state when the system 
% reaches equilibrium. It tells us how often, on average, we would expect the 
% system to be in each state after a long enough period.

% Linearization of the system of equations (eye is the identity matrix of
% size nshocks)
% pi1 = 0.73 * pi1 + 0.27 * pi2
% pi2 = 0.27 * pi1 + 0.73 * pi2

prob = QQ - eye(nshocks);
% [(0.73 - 1)  0.27]
% [0.27  (0.73 - 1)]

% but this has infinite solutions, hence we first enforce pi1 + pi2 = 1
% by enforcing the last column to be 1 (normalization condition)
% This is the same as normalizing the eigenvector to sum to 1 in tauchen's
% [-0.27  1]
% [ 0.27  1]

prob(:, nshocks) = ones(nshocks, 1);

% On the other side of the equality, we would have a = [0 1] to solve
% for the ergodic distribution over states that remain unchanged over the
% long run and sum up to 1

a = zeros(nshocks, 1);
a(nshocks) = 1;

% Ergodic distribution
epr = linsolve(prob', a);

% Placeholders for my distributions 
% initial_dist = (psi_0(m1) + [psi_0(mM) - psi_0(m1)] * unif_dist) * epr
% but since psi_0(mM)=1 and psi_0(m1)=0
F1 = zeros(Dp, Dk);
F1(:,1) = epr(1) * unif_dist; %epr(1) is for the ergodic dist for theta_1

F2 = zeros(Dp, Dk);
F2(:,1) = epr(2) * unif_dist; %epr(2) is for the ergodic dist for theta_2


% Iterations for convergence
for i = 1:Dk-1
    for j=1:Dp

        % For low shock

        if dist_mgrid(j) < g_pol_1(1)
            val_1 = 0;
        elseif dist_mgrid(j) > g_pol_1(Np)
            val_1 = epr(1);
        else
            % This is to obtain the unique values of the policy
            % function and extract the index for each one of these values.
            [g1r, ind_1] = unique(g_pol_1, "last");

            % For the first interpolation I want an interpolated value for dist_mgrid(j) 
            % using the function values (gr1) of mgrid(ind_1).
            % For the second interpolation I want an interpolated value of
            % the distribution function F1(:,i) at that query point (new dgrid).

            % The reason g1r is X and mgrid(ind_1) is Y is because I want
            % the value in mgrid (inverse policy function used)

            % Otherwise dist_mgrid is X and F1(:,i) is Y because I
            % want the value of the distribution such that I can use it
            % as my new value of the function in i+1.

            val_1 = interp1(dist_mgrid, F1(:,i), interp1(g1r, mgrid(ind_1), dist_mgrid(j), "linear"), "linear");
               
        end

        % For high shock
        
        if dist_mgrid(j) <= g_pol_2(1)
            val_2 = 0;
        elseif dist_mgrid(j) >= g_pol_2(Np)
            val_2 = epr(2);
        else
            [g2r, ind_2] = unique(g_pol_2, "last");

            val_2 = interp1(dist_mgrid, F2(:,i), interp1(g2r, mgrid(ind_2), dist_mgrid(j), "linear"), "linear");
        end
        
        % New distribution
        F1(j, i+1) = val_1 * QQ(1,1) + val_1 * QQ(2,1);
        F2(j, i+1) = val_2 * QQ(1,2) + val_2 * QQ(2,2);
        
    end

    diff1 = max(abs(F1(:,i+1)-F1(:,i)));
    diff2 = max(abs(F2(:,i+1)-F2(:,i)));

    F1_dist = F1(:,i+1);
    F2_dist = F2(:,i+1);

    if diff1 <= 1e-10
        disp("convergence achieved of F1");
        disp(i);
        break
    end

    if diff2 <= 1e-10
        disp("convergence achieved of F2");
        disp(i);
        break
    end
end




% Plot distribution functions when they converge
figure; 
subplot(2,1,1);
plot(mgrid, g_pol_1); hold on;
plot(mgrid, mgrid); % 45-degree line
grid on;
xlabel('Assets');
ylabel('Assets t+1');
legend('Low theta', '45-degree line');
title('Policy function: Assets');

subplot(2,1,2);
plot(dist_mgrid, F1(:,1)); hold on;
plot(dist_mgrid, F1(:,2)); hold on;
plot(dist_mgrid, F1(:,3)); hold on;
plot(dist_mgrid, F1(:,4)); hold on;
plot(dist_mgrid, F1_dist);
grid on;
xlabel('dgrid');
ylabel('Psi(m, theta_1)');
title('Distribution functions for low shock');
legend({'Psi_0', 'Psi_1', 'Psi_2', 'Psi_3', 'Last distribution'}, 'location', 'southeast');

figure; 
subplot(2,1,1); 
plot(mgrid, g_pol_2); hold on;
plot(mgrid, mgrid); % 45-degree line
grid on;
xlabel('Assets');
ylabel('Assets t+1');
legend('High theta', '45-degree line');
title('Policy function: Assets');

subplot(2,1,2);
plot(dist_mgrid, F2(:,1)); hold on;
plot(dist_mgrid, F2(:,2)); hold on;
plot(dist_mgrid, F2(:,3)); hold on;
plot(dist_mgrid, F2(:,4)); hold on;
plot(dist_mgrid, F2_dist);
grid on;
xlabel('dgrid');
ylabel('Psi(m, theta_2)');
title('Distribution functions for high shock');
legend({'Psi_0', 'Psi_1', 'Psi_2', 'Psi_3', 'Last distribution'}, 'location', 'southeast');


%% Computing the equilibrium

clear; clc; close all;

% Same information as before but we don't have tau

% Parameters
beta = 0.98;
sigma = 2;
y = 1; % Fixed endowment
%tau = 0.0234; we don't have this
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

% Transition matrix
QQ = [p11 p21; p12 p22];

% Creating the placeholder for money supply and the first guess
TT = 200;
m_supply = zeros(TT,1);
m_supply(1) = 1.5; % Arbitrary, first guess

for s=1:TT+1
    
    % Use the guess in monery supply to get Tau
    tau = m_supply(s) * gamma / (1 + gamma);

    Np = 1000; % Size of the grid

    mp_min = y + tau; % CIA constraint (m' = y + tau)
    mp_max = 6;  % Maximum money level (optional)
    
    % Create the grid of m
    mgrid = linspace(mp_min, mp_max, Np)';
    
    % Maximum number of iterations
    Nk = 1000;
    

    % POLICY FUNCTIONS


    % Initial guesses of m
    g0_1 = mgrid;
    g0_2 = mgrid;
    
    % Preallocate memory for future iterations
    g_1 = zeros(Np, Nk);
    g_2 = zeros(Np, Nk);
    
    % First column is the initial guess
    g_1(:,1) = g0_1;
    g_2(:,1) = g0_2;
    
    for j = 1:Nk-1
        uprime_1 = theta_1 * ((mgrid ./ (1 + gamma)) + y + tau - g_1(:,j)) .^ (-sigma);
        uprime_2 = theta_2 * ((mgrid ./ (1 + gamma)) + y + tau - g_2(:,j)) .^ (-sigma);
        
        % Compute expectations under the Markov transition matrix
        EV_1 = ((beta/(1+gamma)) * (p11 * uprime_1' + p12 * uprime_2'))';
        EV_2 = ((beta/(1+gamma)) * (p21 * uprime_1' + p22 * uprime_2'))';
        
        % Compute mstar
        ms_1 = ((EV_1/theta_1) .^ (-1/sigma) + mgrid - y - tau) * (1+gamma);
        ms_2 = ((EV_2/theta_2) .^ (-1/sigma) + mgrid - y - tau) * (1+gamma); 
        
        for i=1:Np
            if mgrid(i) <= ms_1(1)
                g_1(i,j+1) = mgrid(1);
            elseif mgrid(i) >= ms_1(Np)
                g_1(i, j+1) = mgrid(Np);
            else
                g_1(i, j+1) = interp1(ms_1, mgrid, mgrid(i), 'linear', 'extrap');
            end
        end
    
        for i=1:Np
            if mgrid(i) <= ms_2(1)
                g_2(i,j+1) = mgrid(1);
            elseif mgrid(i) >= ms_2(Np)
                g_2(i, j+1) = mgrid(Np);
            else
                g_2(i, j+1) = interp1(ms_2, mgrid, mgrid(i), 'linear', 'extrap');
            end
        end
    
        % Check convergence
        gdiff = max(abs([g_1(:,j+1); g_2(:,j+1)] - [g_1(:,j); g_2(:,j)]));
        if gdiff <= 1e-10
            disp("Convergence achieved after " + num2str(j) + " iterations.");
            g_pol_1 = g_1(:, j+1);
            g_pol_2 = g_2(:, j+1);
            break;
        end
    end
    
    % If convergence isn't achieved, store final iteration values
    if gdiff > 1e-10
        disp("Maximum iterations reached without convergence.");
        g_pol_1 = g_1(:, end);
        g_pol_2 = g_2(:, end);
    end
     
    

    % DISTRIBUTIONS 

    % Distribution parameters
    % New grid for the distribution
    Dp = 2*Np-1;
    dist_mgrid = linspace(mp_min, mp_max, Dp)';
    
    Dk = 10000; % Maximum number of distribution iterations
    
    % Other options: normal dist
    unif_dist = (dist_mgrid - dist_mgrid(1)) / (dist_mgrid(Dp) - dist_mgrid(1));
    
    
    % Create a transition matrix with the prob
    QQ = [p11 p21; p12 p22];
    nshocks = 2;
    
    % Create the ergodic distribution (cn 7)
    % The ergodic distribution (or stationary distribution) of a Markov chain 
    % represents the long-run probability of being in each state when the system 
    % reaches equilibrium. It tells us how often, on average, we would expect the 
    % system to be in each state after a long enough period.
    
    % Linearization of the system of equations (eye is the identity matrix of
    % size nshocks)
    % pi1 = 0.73 * pi1 + 0.27 * pi2
    % pi2 = 0.27 * pi1 + 0.73 * pi2
    
    prob = QQ - eye(nshocks);
    % [(0.73 - 1)  0.27]
    % [0.27  (0.73 - 1)]
    
    % but this has infinite solutions, hence we first enforce pi1 + pi2 = 1
    % by enforcing the last column to be 1 (normalization condition)
    % This is the same as normalizing the eigenvector to sum to 1 in tauchen's
    % [-0.27  1]
    % [ 0.27  1]
    
    prob(:, nshocks) = ones(nshocks, 1);
    
    % On the other side of the equality, we would have a = [0 1] to solve
    % for the ergodic distribution over states that remain unchanged over the
    % long run and sum up to 1
    
    a = zeros(nshocks, 1);
    a(nshocks) = 1;
    
    % Ergodic distribution
    epr = linsolve(prob', a);
    
    % Placeholders for my distributions 
    % initial_dist = (psi_0(m1) + [psi_0(mM) - psi_0(m1)] * unif_dist) * epr
    % but since psi_0(mM)=1 and psi_0(m1)=0
    F1 = zeros(Dp, Dk);
    F1(:,1) = epr(1) * unif_dist; %epr(1) is for the ergodic dist for theta_1
    
    F2 = zeros(Dp, Dk);
    F2(:,1) = epr(2) * unif_dist; %epr(2) is for the ergodic dist for theta_2
    
    
    % Iterations for convergence
    for i = 1:Dk-1
        for j=1:Dp
    
            % For low shock
    
            if dist_mgrid(j) < g_pol_1(1)
                val_1 = 0;
            elseif dist_mgrid(j) > g_pol_1(Np)
                val_1 = epr(1);
            else
                % This is to obtain the unique values of the policy
                % function and extract the index for each one of these values.
                [g1r, ind_1] = unique(g_pol_1, "last");
    
                % For the first interpolation I want an interpolated value for dist_mgrid(j) 
                % using the function values (gr1) of mgrid(ind_1).
                % For the second interpolation I want an interpolated value of
                % the distribution function F1(:,i) at that query point (new dgrid).
    
                % The reason g1r is X and mgrid(ind_1) is Y is because I want
                % the value in mgrid (inverse policy function used)
    
                % Otherwise dist_mgrid is X and F1(:,i) is Y because I
                % want the value of the distribution such that I can use it
                % as my new value of the function in i+1.
    
                val_1 = interp1(dist_mgrid, F1(:,i), interp1(g1r, mgrid(ind_1), dist_mgrid(j), "linear"), "linear");
                   
            end
    
            % For high shock
            
            if dist_mgrid(j) <= g_pol_2(1)
                val_2 = 0;
            elseif dist_mgrid(j) >= g_pol_2(Np)
                val_2 = epr(2);
            else
                [g2r, ind_2] = unique(g_pol_2, "last");
    
                val_2 = interp1(dist_mgrid, F2(:,i), interp1(g2r, mgrid(ind_2), dist_mgrid(j), "linear"), "linear");
            end
            
            % New distribution
            F1(j, i+1) = val_1 * QQ(1,1) + val_2 * QQ(2,1);
            F2(j, i+1) = val_1 * QQ(1,2) + val_2 * QQ(2,2);
            
        end
    
        diff1 = max(abs(F1(:,i+1)-F1(:,i)));
        diff2 = max(abs(F2(:,i+1)-F2(:,i)));
    
        F1_dist = F1(:,i+1);
        F2_dist = F2(:,i+1);
    
        if diff1 <= 1e-10
            disp("convergence achieved of F1");
            disp(i);
            break
        end
    
        if diff2 <= 1e-10
            disp("convergence achieved of F2");
            disp(i);
            break
        end
    end


    % EQUILIBRIUM

    dopt = zeros(Dp, nshocks);
    
    dopt(:,1) = F1_dist;
    dopt(:,2) = F2_dist;
    
    % Added the two optimal distribution functions (make them a single
    % distribution)
    FO = sum(dopt, 2);
    % Average money holdings (Real balances = money demand)
    % Cobbweb way of finding equilibrium
    RB = (FO(2:Dp) - FO(1:Dp-1))' * ((dist_mgrid(2:Dp) + dist_mgrid(1:Dp-1))/2) + (FO(1) * dist_mgrid(1));
    
    % This is why the initial for loop is up to TT+1 (money supply)
    % money supply = money demand
    m_supply(s+1) = RB;

    if abs(m_supply(s+1) - m_supply(s)) <= 1e-10
        disp("money supply star found")
        disp(m_supply(s));
        disp(s+1);
        break;
    end
    
end




% Plotting money policy function
figure;
plot(mgrid, g_pol_1, 'r'); hold on;
plot(mgrid, g_pol_2, 'b');
plot(mgrid, mgrid, 'k--'); % 45-degree line
xlabel('Money');
ylabel('Money t+1');
legend('Low theta', 'High theta', '45-degree line');
title('Policy function: Money');

% Plotting distribution function
figure;
plot(dist_mgrid, FO);
grid on;
xlabel('Money');
ylabel('Psi(m, theta)');
title('Distribution functions');
legend({'Distribution'}, 'location', 'southeast');


% Compute policy functions for consumption
c_pol_1 = (mgrid / (1 + gamma)) + y + tau - g_pol_1; 
c_pol_2 = (mgrid / (1 + gamma)) + y + tau - g_pol_2;  

% Interpolated consumption defined on the grid dgrid
icons_1 = interp1(mgrid, c_pol_1, dist_mgrid);
icons_2 = interp1(mgrid, c_pol_2, dist_mgrid);

agicons_1 = (dopt(2:Dp, 1) - dopt(1:Dp-1, 1))' * ((icons_1(2:Dp) + icons_1(1:Dp-1))/2) + (dopt(1, 1) * icons_1(1));
agicons_2 = (dopt(2:Dp, 2) - dopt(1:Dp-1, 2))' * ((icons_2(2:Dp) + icons_2(1:Dp-1))/2) + (dopt(1, 2) * icons_2(1));

%Aggregate consumption
agicons = agicons_1 + agicons_2;

% Plotting consumption function
figure;
plot(mgrid, c_pol_1, 'r'); hold on;
plot(mgrid, c_pol_2, 'b');
xlabel('Money');
ylabel('Consumption');
legend('Low theta', 'High theta');
title('Policy function: Consumption');