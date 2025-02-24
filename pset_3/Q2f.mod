var k, c;
parameters beta, alpha, delta, sigma;
predetermined_variables k;

beta = 0.95;    // Discount factor
alpha = 0.36;   // Capital share
delta = 0.025;  // Depreciation rate
sigma = 2;      // Risk Aversion

model;

[name='Euler']
c^(-sigma)/(c(+1)^(-sigma)) = beta*(alpha*k(+1)^(alpha-1)+1-delta);

[name='resource constraint']
k(+1)+c=k^(alpha)+(1-delta)*k;

end;

initval;
    k= 0.1*10.991;
    c= 0.1*2.0954;
end;

endval;
    k= 10.991;
    c= 2.0954;
end;


resid;

perfect_foresight_setup(periods=200);

perfect_foresight_solver;
//perfect_foresight_solver(linear_approximation);

//rplot k;
//rplot c;