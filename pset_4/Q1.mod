var k z c;
varexo e;

parameters beta alpha sigma delta rho;

load params_model;
set_param_value('beta', beta);
set_param_value('alpha', alpha);
set_param_value('sigma', sigma);
set_param_value('delta', delta);
set_param_value('rho', rho);

model;

c = z*k(-1)^alpha+(1-delta)*k(-1)-k;

1 = beta*(c/c(+1))^sigma*(alpha*z(+1)*k^(alpha-1)+1-delta);

log(z) = rho*log(z(-1))+e;

end;

initval;
k = k_b;
c = c_b;
z = 1;
end;

resid;

steady(solve_algo=2);

check;

shocks;
var e = ve^2;
end;

stoch_simul(irf=200,order=1,pruning);