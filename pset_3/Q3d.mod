var k, c, l;
parameters beta, alpha, delta, gamma, k_ss, c_ss, l_ss;
predetermined_variables k;

beta = 0.95;    // Discount factor
alpha = 0.36;   // Capital share
delta = 0.025;  // Depreciation rate

l_ss = 1/3; 

lhs = (1/beta) - (1 - delta);
k_ss = ( lhs / ( alpha * l_ss^(1-alpha) ) )^(1/(alpha-1));
c_ss = k_ss^alpha * l_ss^(1 - alpha) - delta * k_ss;

gamma = ( (1 - alpha)*k_ss^alpha ) / ( c_ss * l_ss^alpha );

model;

[name='Euler']
c(+1)/c=beta*(alpha*k^(alpha-1)*l^(1-alpha)+1-delta);
    
[name='labor']
gamma*c = (1-alpha)*(k(-1)^alpha)*(l^(-alpha));
    
[name='resource constraint']
k(+1)+c=k^(alpha)*l^(1-alpha)+(1-delta)*k;

end;

initval;
    k = 0.1*k_ss;
    c = 0.1*c_ss;
    l = 0.1*l_ss;
end;

endval;
    k= k_ss;
    c= c_ss;   
    l= l_ss;
end;

resid;

perfect_foresight_setup(periods=200);

//perfect_foresight_solver;
perfect_foresight_solver(linear_approximation);

//rplot k;
//rplot c;
//rplot l;