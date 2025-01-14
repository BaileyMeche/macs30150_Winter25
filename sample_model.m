// Define model name
var y c k i h w r z;
varexo epsilon;

// Parameters
parameters beta rho alpha delta theta psi sigma;

// Parameter values
beta    = 0.99;          // Discount factor
rho     = 0.95;          // AR(1) persistence of technology
alpha   = 0.36;          // Capital share
delta   = 0.025;         // Depreciation rate
theta   = 2.0;           // Labor elasticity
psi     = 0.5;           // Coefficient on labor disutility
sigma   = 2.0;           // Risk aversion

// Model equations
model;
    // Euler equation
    1/c = beta * (1/c(+1)) * (1 + r(+1) - delta);

    // Capital accumulation
    k = (1 - delta) * k(-1) + i;

    // Production function
    y = z + alpha * k(-1) + (1 - alpha) * h;

    // Labor market condition
    w = (1 - alpha) * exp(z) * h^(alpha);

    // Wage equation
    w = psi * c * h^theta;

    // Resource constraint
    y = c + i;

    // Technology process
    z = rho * z(-1) + epsilon;
end;

// Initial values
initval;
    y = 1;
    c = 0.8;
    k = 10;
    i = 0.2;
    h = 0.3;
    w = 0.8;
    r = 0.04;
    z = 0;
    epsilon = 0;
end;

// Shocks
shocks;
    var epsilon; stderr 0.01;
end;

// Steady-state computation
steady;
check;

// Simulation and IRFs
stoch_simul(order=1, irf=20, periods=200);