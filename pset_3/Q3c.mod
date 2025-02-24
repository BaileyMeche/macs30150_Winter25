var k, c, l;
parameters beta, alpha, delta, gamma;

beta = 0.95;    // Discount factor
alpha = 0.36;   // Capital share
delta = 0.025;  // Depreciation rate
gamma = 2.17;

model;
[name='Euler']
c(+1)/c=beta*(alpha*k^(alpha-1)*l^(1-alpha)+1-delta);

[name='resource constraint']
k+c=k(-1)^(alpha)*l^(1-alpha)+(1-delta)*k(-1);

[name='labor']
gamma*c = (1-alpha)*(k(-1)^alpha)*(l^(-alpha));
end;

initval;
    k = 0.1*3.672;
    c = 0.1*0.6989;
    l = 0.1*(1/3);
end;

endval;
    k= 3.672;
    c= 0.6989;   
    l= 1/3;
end;

resid;

perfect_foresight_setup(periods=200);

perfect_foresight_solver;
//perfect_foresight_solver(linear_approximation);

//rplot k;
//rplot c;
//rplot l;