// %%%%%%%%%%%%%%%%------ FILE: PQuestion3G.mod ---------%%%%%%%%%%%%%%%%%
// Problem 3: Linear Approximation

var k, c, l;                  // State and control variables
parameters beta, alpha, delta, gamma, k_ss, c_ss, l_ss; // Model parameters and steady-state values
predetermined_variables k;    // Capital is predetermined

// Parameter values
alpha  = 0.36;
beta   = 0.95;
delta  = 0.025;
l_ss = 1/3; 

% (i) => alpha * k_ss^(alpha - 1) * l_ss^(1 - alpha) = 1/beta - (1-delta).
lhs = (1/beta) - (1 - delta);
% => k_ss^(alpha-1) = lhs / [ alpha * l_ss^(1-alpha) ]
k_ss = ( lhs / ( alpha * l_ss^(1-alpha) ) )^(1/(alpha-1));  % watch sign of (alpha-1) < 0

% (ii):
c_ss = k_ss^alpha * l_ss^(1 - alpha) - delta * k_ss;

% (iii) => gamma = ( (1-alpha)*k_ss^alpha ) / ( c_ss * l_ss^alpha )
gamma = ( (1 - alpha)*k_ss^alpha ) / ( c_ss * l_ss^alpha );

model;
    // Euler Equation
    (1/c) = beta*(1/c(+1))*( alpha*k(+1)^(alpha-1)*l(+1)^(1-alpha) + 1 - delta );
    
    // Labor Supply
    c = ((1 - alpha) * k^alpha) / (gamma * l^alpha);
    
    // Resource Constraint
    k(+1) + c = k^alpha * l^(1 - alpha) + (1 - delta) * k;
end;

initval;
    k =0.1*k_ss; // Initial capital at steady state
    c = 0.5*c_ss; // Initial consumption at steady state
    l = 0.5*l_ss; // Initial labor at steady state
end;

endval;
    k =k_ss; // Initial capital at steady state
    c =c_ss; // Initial consumption at steady state
    l = l_ss; // Initial labor at steady state
end;

//check; // Check steady-state eigenvalues for stability

resid; // Check residuals for equations//

// Linear Approximation
perfect_foresight_setup(periods=200);
perfect_foresight_solver(linear_approximation);

// Linear results
rplot k;
rplot c;
rplot l;
