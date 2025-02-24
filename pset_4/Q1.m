%Params

alpha = 0.36;
beta = 0.95;
delta = 0.025;
sigma = 2;
rho = 0.95;

save params_model alpha beta delta sigma rho

%Analytical Steady State
z_b = 1;
k_b = ((1/beta - 1 + delta)/alpha)^(1/(alpha-1));
c_b = k_b^alpha-delta*k_b;

%Solutions at steady state
uc = c_b^(-sigma);
ucc = -sigma*c_b^(-sigma-1);
fk = alpha*k_b^(alpha-1)+1-delta;
fkk = alpha*(alpha-1)*k_b^(alpha-2);
fkz = alpha*k_b^(alpha-1);
fz = k_b^alpha;


A = [1/beta -1
    -((uc*fkk)/ucc) 1+(beta*uc*fkk)/ucc];

B = [fz
    -beta*uc/ucc*(fkk*fz+fkz*rho)]

[V,D] = eig(A);

%Sorting Values and vectors

[E, ind] = sort(diag(D));
DS = diag(E)
VS = V(:,ind)
IV = inv(VS);

C = IV*B

% Calculate policy functions
c_k = -IV(2,1)/IV(2,2)
c_z = -C(2)/((DS(2,2)-rho)*IV(2,2))

k_k = 1/beta + IV(2,1)/IV(2,2)
k_z = fz + C(2)/((DS(2,2)-rho)*IV(2,2))

% Simulate

T = 200;
ve = .1;

et = zeros(T,1);
et(2) = ve;

CT = zeros(T,1);
KT = zeros(T,1);
ZT = zeros(T,1);

for i=1:T-1;
    ZT(i+1) = rho*ZT(i)+et(i+1);
    CT(i+1) = c_k*KT(i) + c_z*ZT(i+1);
    KT(i+1) = k_k*KT(i) + k_z*ZT(i+1);
end;

figure;
subplot(3,1,1);
plot(0:T-1, CT, 'r', 'DisplayName', 'Consumption');
title('Consumption');
ylabel('Consumption');
xlabel('Time');
grid on;

subplot(3,1,2);
plot(0:T-1, KT, 'r', 'DisplayName', 'Capital');
title('Capital');
ylabel('Capital');
xlabel('Time');
grid on;

subplot(3,1,3);
plot(0:T-1, ZT, 'r', 'DisplayName', 'Productivity Level');
title('Productivity Level');
ylabel('Productivity Level');
xlabel('Time');
grid on;