% a) Construct a finite 5-state approximation for the Markov transition 
% of the process lnzt, using Tauchen’s method. Use 3 standard deviations 
% to approximate the lower and upper end of the grid (λ = 3). 
% Find the ergodic distribution for this process.
clc; 
clear;

% Parameters
N = 5; % 5- state approx 
m = 3;
alpha = 0.36;
beta = 0.98; 
delta = 0.03; 
gamma = 2;
sigma = 0.1;
rho = 0.7; % zt+1 = rho*zt + e (so rho is the lag term on old shock)

% mean:
mu = 0; 

% discretized state space Z and transition matrix Q
[Z, transition_matrix] =  tauchen(N, mu, rho, sigma, m);

transition_matrix_transpose = transpose(transition_matrix);
[V,D] = eig(transition_matrix_transpose);
pi_vector_invariant_dist = V(:, 1); %choose one that corresponds to eignevalue = 1

% Plotting to check:
pi_vector_norm = pi_vector_invariant_dist/sum(pi_vector_invariant_dist); %to make it sum up to 1 
plot(Z, pi_vector_norm);

% b) With the discretized process for zt found in question a), solve the 
% model using VFI. Construct the grid with 1000 points with the lower 
% bound for the capital stock of 0.01k ̄ and the upper bound 4k ̄, 
% where k ̄ is the non-stochastic steady state of the capital stock. 
% Find the value functions and the policy functions. 

% In addition to the above parameters:

Np = 1000; % # of points for evaluation
Nk = 5000; % max number of iterations allowed

kss = (alpha*beta/ (1-(1-delta)*beta))^(1/(1-alpha)); % capital, non-stochastic SS
kgrid = linspace (0.01*kss, 4*kss, Np)'; % grid for approximation w/ lower bound 0.01*kss and upper bound 4*kss

%NOTE TO SELF: kgrid represents the ALL possible values of capital that can
%happen 

z = exp(Z); % shocks from ln z from Tauchen method

cons1 = z(1)*kgrid.^alpha+(1-delta)*kgrid-kgrid';
cons2 = z(2)*kgrid.^alpha+(1-delta)*kgrid-kgrid';
cons3 = z(3)*kgrid.^alpha+(1-delta)*kgrid-kgrid';
cons4 = z(4)*kgrid.^alpha+(1-delta)*kgrid-kgrid';
cons5 = z(5)*kgrid.^alpha+(1-delta)*kgrid-kgrid';

% If any of the above turn out negative, assign them 0 instead:
cons1 (cons1<0) = 0;
cons2 (cons2<0) = 0;
cons3 (cons3<0) = 0;
cons4 (cons4<0) = 0;
cons5 (cons5<0) = 0;

% Different Matrices of returns:

Ret1 = u_t(gamma, cons1);
Ret2 = u_t(gamma, cons2); 
Ret3 = u_t(gamma, cons3);
Ret4 = u_t(gamma, cons4);
Ret5 = u_t(gamma, cons5);

% Different Matrices as placeholders for the Value Functions:
v1 = zeros(Np,Nk); 
v2 = zeros(Np,Nk); 
v3 = zeros(Np,Nk); 
v4 = zeros(Np,Nk); 
v5 = zeros(Np,Nk); 

% Establishing the initial guesses:
v1(:,1) = u_t(gamma, kgrid); 
v2(:,1) = u_t(gamma, kgrid);
v3(:,1) = u_t(gamma, kgrid); 
v4(:,1) = u_t(gamma, kgrid);
v5(:,1) = u_t(gamma, kgrid);

for i=1:Nk-1
    v1(:, i+1) = max((Ret1+beta*transition_matrix(1,1)*v1(:,i)'+ beta*transition_matrix(1,2)*v2(:,i)'+beta*transition_matrix(1,3)*v3(:,i)'+beta*transition_matrix(1,4)*v4(:,i)'+beta*transition_matrix(1,5)*v5(:,i)')')'; 
    v2(:, i+1) = max((Ret2+beta*transition_matrix(2,1)*v1(:,i)'+ beta*transition_matrix(2,2)*v2(:,i)'+beta*transition_matrix(2,3)*v3(:,i)'+beta*transition_matrix(2,4)*v4(:,i)'+beta*transition_matrix(2,5)*v5(:,i)')')'; 
    v3(:, i+1) = max((Ret3+beta*transition_matrix(3,1)*v1(:,i)'+ beta*transition_matrix(3,2)*v2(:,i)'+beta*transition_matrix(3,3)*v3(:,i)'+beta*transition_matrix(3,4)*v4(:,i)'+beta*transition_matrix(3,5)*v5(:,i)')')'; 
    v4(:, i+1) = max((Ret4+beta*transition_matrix(4,1)*v1(:,i)'+ beta*transition_matrix(4,2)*v2(:,i)'+beta*transition_matrix(4,3)*v3(:,i)'+beta*transition_matrix(4,4)*v4(:,i)'+beta*transition_matrix(4,5)*v5(:,i)')')'; 
    v5(:, i+1) = max((Ret5+beta*transition_matrix(5,1)*v1(:,i)'+ beta*transition_matrix(5,2)*v2(:,i)'+beta*transition_matrix(5,3)*v3(:,i)'+beta*transition_matrix(5,4)*v4(:,i)'+beta*transition_matrix(5,5)*v5(:,i)')')'; 
    v = [v1(:,i+1); v2(:,i+1); v3(:, i+1);v4(:, i+1) ;v5(:, i+1)] -[v1(:,i);v2(:, i);v3(:, i) ;v4(:, i);v5(:, i)]; 
        if max(abs(v)) <= 0.001*(1 - beta) % using 0.001*(1 - beta) as epsilon
            % storing convergent value functions & capital index
            [v1s, kindex1] = max((Ret1+beta*transition_matrix(1,1)*v1(:,i)'+ beta*transition_matrix(1,2)*v2(:,i)'+beta*transition_matrix(1,3)*v3(:,i)'+beta*transition_matrix(1,4)*v4(:,i)'+beta*transition_matrix(1,5)*v5(:,i)')'); 
            [v2s, kindex2] = max((Ret2+beta*transition_matrix(2,1)*v1(:,i)'+ beta*transition_matrix(2,2)*v2(:,i)'+beta*transition_matrix(2,3)*v3(:,i)'+beta*transition_matrix(2,4)*v4(:,i)'+beta*transition_matrix(2,5)*v5(:,i)')'); 
            [v3s, kindex3] = max((Ret3+beta*transition_matrix(3,1)*v1(:,i)'+ beta*transition_matrix(3,2)*v2(:,i)'+beta*transition_matrix(3,3)*v3(:,i)'+beta*transition_matrix(3,4)*v4(:,i)'+beta*transition_matrix(3,5)*v5(:,i)')'); 
            [v4s, kindex4] = max((Ret4+beta*transition_matrix(4,1)*v1(:,i)'+ beta*transition_matrix(4,2)*v2(:,i)'+beta*transition_matrix(4,3)*v3(:,i)'+beta*transition_matrix(4,4)*v4(:,i)'+beta*transition_matrix(4,5)*v5(:,i)')'); 
            [v5s, kindex5] = max((Ret5+beta*transition_matrix(5,1)*v1(:,i)'+ beta*transition_matrix(5,2)*v2(:,i)'+beta*transition_matrix(5,3)*v3(:,i)'+beta*transition_matrix(5,4)*v4(:,i)'+beta*transition_matrix(5,5)*v5(:,i)')');

            disp("Convergence!");
            break;
        end
end 


% Storing all of our states' kindex matrices/vectors in one big matrix:
kindex = [kindex1' kindex2' kindex3' kindex4' kindex5'];

% Placeholder for Capital's policy function:
kpolfn = zeros(Np,5);

for j = 1:5
    for i =1:Np
        kpolfn(i, j) =kgrid(kindex(i,j));
    end
end

%NOTE TO SELF: kpolfn represents a specific path of capital that
%materializes for each of the 5 states 

% Storing all value functions matrices in one big matrix:
vals = [v1s' v2s' v3s' v4s' v5s']; 

% Plotting to check:
figure; 
plot(kgrid, vals);
grid;
xlabel('k');
ylabel('values');
title('Value functions');
legend({'v1', 'v2', 'v3', 'v4', 'v5'}, 'location', 'southeast');

figure;
plot(kgrid, [kpolfn kgrid]); %kgrid & kgrid is for 45 degrees line
grid;
xlabel('k');
ylabel('k');
title('Policy functions');
legend({'g1', 'g2', 'g3', 'g4', 'g5', '45 degree'}, 'location', 'southeast');

% c) With the policy functions for capital obtained in b), simulate
% the stock of capital. To do this depart from an arbitrary capital stock
% in the grid and simulate the draws of the Markov chain for the shock. 
% You will need to use a random number generator. 
% Simulate 5000 observations for the capital stock and compute the 
% histogram. This should be close to the ergodic distribution of capital!

observations = 10000;

% First: simulating the Markov Chain Shocks 
z_sim = zeros(observations, 1); %placeholder for simulated markov chain of shocks
z_sim(1) = 5; %arbitrary starting point (from 1-5)

for i = 1:observations-1
    % Get the probability distribution for the current state
    probabilities = transition_matrix(z_sim(i), :); 
    
    % Normalize in case the probabilities do not sum to 1
    probabilities = probabilities / sum(probabilities);
    
    % Generate a random number between 0 and 1
    rand_num = rand;
    
    % Find the next state using inverse transform sampling
    cumulative_probabilities = cumsum(probabilities);  % Cumulative sum of probabilities
    next_state = find(rand_num <= cumulative_probabilities, 1);  % Find the state
    
    % Update the Markov chain state
    z_sim(i+1) = next_state;
end


% Second: simulating the capital 
k_sim = zeros(observations,1); % placeholder for capital 
k_sim(1) = kgrid(100); % some arbitrary capital stock

for i=1:observations-1
    k_sim(i+1) = kpolfn(kgrid==k_sim(i), z_sim(i+1));
end

% find() function finds the index at which the value of k_sim(i) is
% materialized in kgrid, then uses that in addition to the shock in the
% k policy function (1000 x 5 grid) to get the next value of capital.

% Plotting to check;
figure; 
plot(1:observations, k_sim); 
grid;
xlabel('Time');
ylabel('Capital');
title('Simulating Capital');

% Histogram plot

k_average = mean(k_sim(1000:5000)); 

figure;
histogram(k_sim);
xline(k_average);
grid;
xlabel("Capital");
ylabel("Frequency");
title("Capital histogram");

function result = u_t(gamma, c)
    % Check if gamma is close enough to 1 to consider it as 1
    if abs(gamma - 1) < eps
        result = log(c); % Utility for the case gamma is 1
    else
        result = (c.^(1-gamma))./(1-gamma); % CRRA utility for gamma not equal to 1
    end
end