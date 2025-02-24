var k c l z;
varexo e; 
predetermined_variables k;

parameters beta alpha delta rho gamma sigma eta zbar;

// Define parameter values
beta = 0.98;
alpha = 0.36;
delta = 0.025;
rho = 0.95;
//gamma = 0.3851;    // Assuming sigma = 1
zbar = 1;     // Steady-state value of z
sigma = 1;

k_b = (1/3)*(((1/beta)-1+delta)/alpha)^(1/(alpha-1));
c_b = (1/3)^(1-alpha)*k_b^alpha-delta*k_b;
gamma = c_b/((2/3)*(1-alpha)*(3*k_b)^alpha+c_b);



// Model equations 
model;

//c = z*k(-1)^alpha*l^(1-alpha)+(1-delta)*k(-1)-k;
//1 = beta*(c/c(+1))*(alpha*z(+1)*k^(alpha-1)*l(+1)^(1-alpha)+1-delta);
//log(z)=rho*log(z(-1))+e;
//(1-gamma)*c = gamma*(1-l)*z*(1-alpha)*(k(-1)^alpha)*(l^(-alpha));

gamma*(1-l)*(1-alpha)*z*(k/l)^alpha=(1-gamma)*c;

c^(gamma-1)*(1-l)^(1-gamma)*(c^gamma*(1-l)^(1-gamma))^(-sigma) = beta*c(+1)^(gamma-1)*(1-l(+1))^(1-gamma)*(c(+1)^gamma*(1-l(+1))^(1-gamma))^(-sigma)*(alpha*z(+1)*(k(+1)/l(+1))^(alpha-1)+1-delta);

log(z)=rho*log(z(-1))+e;

k(+1)+c = z*k^alpha*l^(1-alpha)+(1-delta)*k;

end;


// Calculate steady state
steady_state_model;
    l = 1/3;
    k = (1/3)*(((1/beta)-1+delta)/alpha)^(1/(alpha-1));
    c = (1/3)^(1-alpha)*k^alpha-delta*k;
    z = zbar;
end;

// Initialize at steady state
initval;
    l = 1/3;
    k = 8.4690;
    c = 0.8565;
    z = zbar;
end;

// Check steady state
steady;

// Set shock properties
shocks;
    var e = 0.1^2;
end;

// Compute policy and transition functions
stoch_simul(irf=200,order=1,pruning);
//  stoch_simul(order=1, irf=40);

