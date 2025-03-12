clc; 
clear;

beta = 0.95;
sigma = 2;
r = 0.2;
Np = 1000;
Nk = 1000;

wgrid = linspace(0.1, 5, Np)';

%% Analytical solution

B = (1-(((1+r)^(1-sigma))*beta)^(1/sigma))^(-sigma);
X = (((1+r)^(1-sigma))*beta*B)^(-1/sigma)/(1+(((1+r)^(1-sigma))*beta*B)^(-1/sigma));

%Analytical solutions of w and c
c_an = X*wgrid;
c_an(c_an<0) = 0;
wprime_an = (wgrid-c_an)*(1+r);

% From the guess and verify method
v_guess = (1/(1-sigma))*B*wgrid .^ (1-sigma);

%% Value Function Iteration 
%with my resource constraint
consu = wgrid - (wgrid'/(1+r));

consu(consu<0) = 0;

Ret = ut(sigma, consu);

V = zeros(Np,Nk);  % Initial guess for VFI

%Im using this guess to start (the first column of the value function)
V(:,1) = ut(sigma, wgrid);

%any v should work but the closer you are to the solution the better

for i=1:Nk-1
    %Bellman equation. v is the same vector but you are evaluating the k_j
    %the last transpose is because Gauss goes col by col
    V(:,i+1) = max((Ret+beta*V(:,i)')')';
        if max(abs(V(:,i+1) - V(:,i))) <= 0.001*(1-beta);
            [value, windex] = max((Ret+beta*V(:,i)')');
            disp(value);
            disp(windex);
            disp("convergence achieved");
            disp(i);
            break;
        end
end

% Compute the policy function
wprime_value = wgrid(windex); % Maps indices to actual capital values

% Extract optimal consumption using correct indexing
c_value = zeros(Np, 1); % Preallocate for efficiency
for i = 1:Nk
    c_value(i) = consu(i,windex(i));
end

%% Plot Value Functions (Part c)
figure;
plot(wgrid, v_guess, 'r', 'LineWidth', 2); hold on;
plot(wgrid, value, 'b--', 'LineWidth', 2);
legend('Analytical', 'Numerical (VFI)');
title('Comparison of Value Functions');

%% Plot Value Functions (Part c)
figure;
plot(wgrid, c_an, 'r', 'LineWidth', 2); hold on;
plot(wgrid, c_value, 'b--', 'LineWidth', 2);
legend('Analytical', 'Numerical (VFI)');
title('Comparison of Consumption Functions');

%% Plot Value Functions (Part c)

line = linspace(0, 5, Np);
y = line;

figure;
plot(wgrid, wprime_an, 'r', 'LineWidth', 2); hold on;
plot(wgrid, wprime_value, 'b--', 'LineWidth', 2); hold on;
plot(wgrid, y, 'black', 'LineWidth', 1)
legend('Analytical', 'Numerical (VFI)');
title('Comparison of Capital Functions');


%% Orthogonal collocation method for Cake-Eating Problem

% Parameters
r = 0.2;  
sigma = 2; 
beta = 0.95;

param = [r, sigma, beta];

% First guess for coefficients
a0 = [0.2; 0.1; 0.01; 0.001; 0.0001; 0.00001; 0.000001];  
n = length(a0);
l = 1:1:n;

% Chebyshev nodes (important to transpose z)
z = cos((2*(1:n)' - 1)*pi/(2*n));

% Define the bounds for wealth w
w_min = min(wgrid);
w_max = max(wgrid);

% Transform Chebyshev nodes to range [w_min, w_max]
w = ((z+1)*(w_max-w_min)/2)+w_min;

%This is phi_i(k)
aa = T(w);

%Now we multiply by a to get the function g(k; a)
bb = g(w, a0);

%Finally, we compute the residual function
cc = R(w, a0, param);

% Solve for coefficients using nonlinear least squares
as = lsqnonlin(@(a) R(w, a, param), a0);

% Evaluate Chebyshev polynomials on wgrid
t = 2*((wgrid - w_min)./(w_max-w_min))-1;
Ts = zeros(length(wgrid), n);
for i = 0:n-1
    Ts(:, i+1) = wgrid .* cos(i*acos(t));
end

% Approximate policy function
gpol = Ts * as;

% Plot results
figure;
plot(wgrid, gpol, 'r', 'LineWidth', 2); hold on;
plot(wgrid, wprime_an, 'b--', 'LineWidth', 2); hold on;
plot(wgrid, y, 'black', 'LineWidth', 1)
xlabel('Wealth w');
ylabel('Next-period wealth w''');
legend('Orthogonal Approximation', 'Analytical Solution');
title('Policy Function for w'' using Collocation Method');
grid on;

% Functions
function [Ts] = T(w)
    w_min = 0.1;
    w_max = 5;

    % Compute transformation into Chebyshev space
    t_cheb = 2*((w-w_min)./(w_max-w_min))-1;

    % Initialize the matrix for Chebyshev polynomials
    Ts = zeros(length(t_cheb), length(t_cheb));

    %Compute Chebyshev polynomials for each w
    for i = 0:length(t_cheb)-1
        Ts(:, i+1) = w .* cos(i * acos(t_cheb));
    end
end


function [g] = g(w, a)
    g = T(w) * a;
end


function [R] = R(w, a, param)
    r = param(1);
    sigma = param(2);
    beta = param(3);

    % Compute the policy function g(w, a)
    gw = g(w, a);
    ggw = g(gw, a); % g(g(w, a), a), the next-period capital policy function
    
    % Compute consumption using feasibility constraint
    c = w - gw / (1 + r);
    
    % Compute the RHS of the Euler equation
    rhs = ( (c.^(-sigma)) / (beta * (1 + r)) ).^(-1/sigma);
    
    % Residual function
    R = gw - rhs - (ggw/(1+r));
end


function[c] = ut(g,b)
%c = zeros(size(b));
    if g< 1 | g>1 
    c = (1/(1-g))*b.^(1-g); 
    else; 
    c = log(b);
    end
end