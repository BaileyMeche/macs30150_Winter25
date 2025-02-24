addpath C:\dynare\6.2\matlab

%Params

alpha = 0.36;
beta = 0.95;
delta = 0.025;
sigma = 2;

k_b = 10.991;
c_b = 2.0954;
k0 = 0.1*k_b;

A = [1.053 -1; -0.004736 1.0045];

[V,D] = eig(A);

%Sorting Values and vectors

[E, ind] = sort(diag(D));
DS = diag(E);
VS = V(:,ind);
IV = inv(VS);

T = 200;
kt = zeros(T,1);
kt(1) = k0;
ka = zeros(T,1);
ka(1) = k0;
ct = zeros(T,1);

for i=1:T-1
    kt(i+1) = k_b+0.9558*(kt(i)-k_b);
    %ka(i+1) = alpha*beta*ka(i)^alpha;
    ct(i) = c_b+0.0972*(kt(i)-k_b);
end

%plot(kt)
%title('Capital')
%plot(ct)
%title('Consumption')

figure;
subplot(2,1,1);
plot(0:T-1, kt);
title('Capital Path');
ylabel('k');
xlabel('Time');
grid on;

subplot(2,1,2);
plot(0:T-1, ct);
title('Consumption Path');
ylabel('c');
xlabel('Time');
grid on;

% Run Dynare nonlinear simulation (Part f)
dynare Q2f.mod noclearall;

% Extract Dynare results
dynare_k = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_c = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';

figure;
subplot(2,1,1);
plot(0:T-1, kt, 'r', 'DisplayName', 'Part e Capital');
hold on;
plot(0:(length(dynare_k)-1), dynare_k, 'b', 'DisplayName', 'Dynare Capital Path');
title('Capital Path');
ylabel('k');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T-1, ct, 'r', 'DisplayName', 'Part e Consumption');
hold on;
plot(0:(length(dynare_c)-1), dynare_c, 'b', 'DisplayName', 'Dynare Consumption Path');
title('Consumption Path');
ylabel('c');
xlabel('Time');
legend show;
grid on;


dynare Q2g.mod noclearall;
% Extract Dynare results
dynare_k = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_c = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';

figure;
subplot(2,1,1);
plot(0:T-1, kt, 'r', 'DisplayName', 'Part e Capital');
hold on;
plot(0:(length(dynare_k)-1), dynare_k, '--b', 'DisplayName', 'Dynare Capital Path: Linear Approx');
title('Capital Path');
ylabel('k');
xlabel('Time');
legend show;
grid on;

subplot(2,1,2);
plot(0:T-1, ct, 'r', 'DisplayName', 'Part e Consumption');
hold on;
plot(0:(length(dynare_c)-1), dynare_c, '--b', 'DisplayName', 'Dynare Consumption Path: Linear Approx');
title('Consumption Path');
ylabel('c');
xlabel('Time');
legend show;
grid on;

dynare Q3c.mod noclearall;
% Extract Dynare results
dynare_k = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_c = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';
dynare_l = oo_.endo_simul(strcmp(M_.endo_names, 'l'), :)';

figure;
plot(0:(length(dynare_k)-1), dynare_k, 'b', 'DisplayName', 'Capital');
hold on;
plot(0:(length(dynare_c)-1), dynare_c, 'r', 'DisplayName', 'Consumption');
plot(0:(length(dynare_l)-1), dynare_l, 'g', 'DisplayName', 'labor');
title('Capital Labor and Consumption');
ylabel('Value');
xlabel('Time');
legend show;
grid on;

dynare Q3d.mod noclearall;
% Extract Dynare results
dynare_kL = oo_.endo_simul(strcmp(M_.endo_names, 'k'), :)';
dynare_cL = oo_.endo_simul(strcmp(M_.endo_names, 'c'), :)';
dynare_lL = oo_.endo_simul(strcmp(M_.endo_names, 'l'), :)';

figure;
plot(0:(length(dynare_k)-1), dynare_k, 'r', 'DisplayName', 'Capital');
hold on;
plot(0:(length(dynare_kL)-1), dynare_kL, '--b', 'DisplayName', 'Approx Capital');
title('Capital');
ylabel('Value');
xlabel('Time');
legend show;
grid on;

figure;
plot(0:(length(dynare_c)-1), dynare_c, 'r', 'DisplayName', 'Consumption');
hold on;
plot(0:(length(dynare_cL)-1), dynare_cL, '--b', 'DisplayName', 'Approx Consumption');
title('Consumption');
ylabel('Value');
xlabel('Time');
legend show;
grid on;

figure;
plot(0:(length(dynare_l)-1), dynare_l, 'r', 'DisplayName', 'Labor');
hold on;
plot(0:(length(dynare_lL)-1), dynare_lL, '--b', 'DisplayName', 'Approx Labor');
title('Labor');
ylabel('Value');
xlabel('Time');
legend show;
grid on;