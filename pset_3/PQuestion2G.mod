%%%%%%%%%%%%%%%%------ FILE: PQuestion2G.mod ---------%%%%%%%%%%
// Nonlinear simulation
var k, c;
parameters beta, sigma, alpha, delta;
predetermined_variables k;

beta = 0.95;
sigma = 2;
alpha = 0.36;
delta = 0.025;

model;
[name='Euler']
(c^(-sigma)) = beta * (c(+1)^(-sigma)) * (alpha * k(+1)^(alpha - 1) + 1 - delta);
[name='Resource Constraint']
k(+1) + c = k^alpha + (1 - delta) * k;
end;

initval;
    k = 0.1 * ((alpha * beta) / (1 - beta * (1 - delta)))^(1 / (1 - alpha));
    c = 0.1;
end;

endval;
    k = ((alpha * beta) / (1 - beta * (1 - delta)))^(1 / (1 - alpha));
    c = k^alpha - delta * k;
end;

resid;
perfect_foresight_setup(periods = 200);
perfect_foresight_solver;

rplot k;
rplot c;

// Linear approximation
var k, c;
parameters beta, sigma, alpha, delta;
predetermined_variables k;

beta = 0.95;
sigma = 2;
alpha = 0.36;
delta = 0.025;

model;
[name='Euler']
(c^(-sigma)) = beta * (c(+1)^(-sigma)) * (alpha * k(+1)^(alpha - 1) + 1 - delta);
[name='Resource Constraint']
k(+1) + c = k^alpha + (1 - delta) * k;
end;

initval;
    k = 0.1 * ((alpha * beta) / (1 - beta * (1 - delta)))^(1 / (1 - alpha));
    c = 0.1;
end;

endval;
    k = ((alpha * beta) / (1 - beta * (1 - delta)))^(1 / (1 - alpha));
    c = k^alpha - delta * k;
end;

resid;
perfect_foresight_setup(periods = 200);
perfect_foresight_solver(linear_approximation);

rplot k;
rplot c;
