close all; clear all; 

% PSet 8 - VFI & Endogenous Grid Method
clear; clc; close all;

% Parameters & Grid Setup
% Economic parameters
price    = 1.0124;
discRate = 0.99322;
rho      = 1.5;

% Endowment levels
end_high = 1;    % High endowment
end_low  = 0.1;  % Low endowment

% Markov transition probabilities
% When starting with high endowment
pi_hh = 0.925;
pi_hl = 1 - pi_hh;
% When starting with low endowment
pi_ll = 0.5;
pi_lh = 1 - pi_ll;

% Asset grid parameters
gridSize = 1000;
amin     = -2; % Borrowing limit
amax     = 4;  % Maximum asset level
assetGrid = linspace(amin, amax, gridSize)';

% Iteration settings
maxIter_VFI = 1000;

%% Part I: Value Function Iteration (VFI)

% Build consumption matrices (rows: current asset, columns: chosen asset)
cons_low  = assetGrid + end_low - price .* assetGrid';
cons_high = assetGrid + end_high - price .* assetGrid';

% Replace negative consumption values with zero
cons_low(cons_low < 0)   = 0;
cons_high(cons_high < 0) = 0;

% Compute immediate returns using our CRRA utility function
R_low  = myUtility(rho, cons_low);
R_high = myUtility(rho, cons_high);

% Preallocate value function matrices for both endowment states
VF_low  = zeros(gridSize, maxIter_VFI);
VF_high = zeros(gridSize, maxIter_VFI);

% Use utility of current asset levels as the initial guess
VF_low(:,1)  = real(myUtility(rho, assetGrid));
VF_high(:,1) = real(myUtility(rho, assetGrid));

% Main iteration loop for VFI
for it = 1:maxIter_VFI-1
    % Compute expected continuation values (row vectors)
    EV_low  = discRate * (pi_ll * VF_low(:,it)' + pi_lh * VF_high(:,it)');
    EV_high = discRate * (pi_hl * VF_low(:,it)' + pi_hh * VF_high(:,it)');
    
    % Choose the asset index that maximizes the sum of immediate and continuation values
    [VF_low(:, it+1), ind_low]  = max((R_low + EV_low)');
    [VF_high(:, it+1), ind_high] = max((R_high + EV_high)');
    
    % Check convergence: if the change in value functions is small enough, exit loop.
    diffVF = [VF_low(:,it+1); VF_high(:,it+1)] - [VF_low(:,it); VF_high(:,it)];
    if max(abs(diffVF)) <= 0.001 * (1 - discRate)
        disp("Convergence achieved in VFI.");
        break;
    end
end

% Derive asset policy functions using the maximizing indices
policyAsset_low  = zeros(gridSize, 1);
policyAsset_high = zeros(gridSize, 1);
for j = 1:gridSize
    policyAsset_low(j)  = assetGrid(ind_low(j));
    policyAsset_high(j) = assetGrid(ind_high(j));
end

% Derive consumption policy functions
policyCons_low  = assetGrid + end_low  - price * policyAsset_low;
policyCons_high = assetGrid + end_high - price * policyAsset_high;

% Plot the value functions
figure;
plot(assetGrid, VF_low(:, it+1), 'b', 'LineWidth', 1.5); hold on;
plot(assetGrid, VF_high(:, it+1), 'r', 'LineWidth', 1.5);
grid on; xlabel('Assets (a)'); ylabel('Value Function');
title('Value Functions');
legend({'Low Endowment', 'High Endowment'}, 'location', 'southeast');
hold off;

% Plot the asset policy functions
figure;
plot(assetGrid, policyAsset_low, 'LineWidth', 1.5); hold on;
plot(assetGrid, policyAsset_high, 'LineWidth', 1.5);
plot(assetGrid, assetGrid, '--', 'LineWidth', 1);
grid on; xlabel('Current Assets (a)'); ylabel('Next Period Assets (a'')');
title('Asset Policy Functions');
legend({'Low Endowment', 'High Endowment', '45-degree Line'}, 'location', 'southeast');
hold off;

% Plot the consumption policy functions
figure;
plot(assetGrid, policyCons_low, 'LineWidth', 1.5); hold on;
plot(assetGrid, policyCons_high, 'LineWidth', 1.5);
grid on; xlabel('Assets (a)'); ylabel('Consumption');
title('Consumption Policy Functions');
legend({'Low Endowment', 'High Endowment'}, 'location', 'southeast');
hold off;

%% Part II: Endogenous Grid Method (EGM)

% Clear variables that will be re-used in the next section, but preserve parameters
clearvars -except price discRate rho end_high end_low pi_hh pi_hl pi_ll pi_lh gridSize amin amax assetGrid

% EGM iteration parameters
tolEGM = 1e-10;
maxIter_EGM = 1000;

% Initial guesses for the future asset policies
futureAsset_low  = assetGrid;
futureAsset_high = assetGrid;

% Preallocate arrays for asset policy iteration
A_low  = zeros(gridSize, maxIter_EGM);
A_high = zeros(gridSize, maxIter_EGM);
A_low(:,1)  = futureAsset_low;
A_high(:,1) = futureAsset_high;

% Preallocate arrays for consumption policies (not used until convergence)
C_low  = zeros(gridSize, maxIter_EGM);
C_high = zeros(gridSize, maxIter_EGM);

for k = 1:maxIter_EGM-1
    % Compute marginal utilities from consumption implied by current asset policies
    margUtil_low  = (assetGrid + end_low - price * A_low(:,k)).^(-rho);
    margUtil_high = (assetGrid + end_high - price * A_high(:,k)).^(-rho);
    
    % Calculate expected marginal utilities
    expMargUtil_low  = (pi_ll * margUtil_low + pi_lh * margUtil_high);
    expMargUtil_high = (pi_hl * margUtil_low + pi_hh * margUtil_high);
    
    % Update consumption policies using the Euler equation inversion
    newCons_low  = (discRate / price * expMargUtil_low).^(-1/rho);
    newCons_high = (discRate / price * expMargUtil_high).^(-1/rho);
    
    % Implied new asset policies
    newA_low  = newCons_low  + price * assetGrid - end_low;
    newA_high = newCons_high + price * assetGrid - end_high;
    
    % Update policies using interpolation back onto the original grid
    A_low(:,k+1)  = interp1(newA_low, assetGrid, assetGrid, 'linear', 'extrap');
    A_high(:,k+1) = interp1(newA_high, assetGrid, assetGrid, 'linear', 'extrap');
    
    % Impose borrowing constraints
    A_low(:,k+1)  = max(amin, min(amax, A_low(:,k+1)));
    A_high(:,k+1) = max(amin, min(amax, A_high(:,k+1)));
    
    % Convergence check
    diffA = max(abs([A_low(:,k+1); A_high(:,k+1)] - [A_low(:,k); A_high(:,k)]));
    if diffA <= tolEGM
        disp("EGM convergence achieved after " + num2str(k) + " iterations.");
        finalA_low  = A_low(:, k+1);
        finalA_high = A_high(:, k+1);
        break;
    end
end

if diffA > tolEGM
    disp("Maximum iterations reached without convergence in EGM.");
    finalA_low  = A_low(:, end);
    finalA_high = A_high(:, end);
end

% Compute final consumption policies from converged asset policies
finalC_low  = assetGrid + end_low  - price * finalA_low;
finalC_high = assetGrid + end_high - price * finalA_high;

% Plot the asset policy functions from EGM
figure;
plot(assetGrid, finalA_low, 'LineWidth', 2); hold on;
plot(assetGrid, finalA_high, 'LineWidth', 2);
plot(assetGrid, assetGrid, '--', 'LineWidth', 1);
xlabel('Assets'); ylabel('Next Period Assets');
legend('Low Endowment', 'High Endowment', '45-degree Line');
title('EGM: Asset Policy Functions');
grid on;
hold off;

% Plot the consumption policy functions from EGM
figure;
plot(assetGrid, finalC_low, 'LineWidth', 2); hold on;
plot(assetGrid, finalC_high, 'LineWidth', 2);
xlabel('Assets'); ylabel('Consumption');
legend('Low Endowment', 'High Endowment');
title('EGM: Consumption Policy Functions');
grid on;
hold off;

%% Local Utility Function Definition
function u = myUtility(rho_val, cons)
    % Ensure consumption is strictly positive (avoid log(0) or negative values)
    cons(cons <= 0) = 1e-6;
    if abs(rho_val - 1) < 1e-6
        u = log(cons);
    else
        u = (cons.^(1 - rho_val)) / (1 - rho_val);
    end
end
