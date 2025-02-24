%Params

alpha = 0.36;
beta = 0.98;
delta = 0.025;
sigma = 1;
rho = 0.95;

save params_model alpha beta delta sigma rho

%Analytical Steady State
k_b = (1/3)*(((1/beta)-1+delta)/alpha)^(1/(alpha-1))
c_b = (1/3)^(1-alpha)*k_b^alpha-delta*k_b
gamma = c_b/((2/3)*(1-alpha)*(3*k_b)^alpha+c_b)


dynare Q2.mod




