%% Q3

clc;
clear all; close all;

% Number of iterations:
Nk = 5000; 

% Parameters
alpha = 0.36;
beta = 0.95; 
sigma = 2;
delta = 0.03;
%phi = 2;
phi = 2; 

% Steady states:
lss = 1/3;
kss = lss*((beta*alpha)/(1+beta*(delta-1)))^(1/(1-alpha));
css = kss^(alpha) * lss^(1 - alpha) - delta*kss; 

% Chi:
chi = (1-alpha)*kss^alpha*lss^(-alpha)*(kss^alpha*lss^(1-alpha)-delta*kss)^(-sigma)*lss^(-1/phi);

% Setting up kgrid:
kgrid = linspace(0.01*kss, 1.2*kss,500)';

% Setting up lgrid; b/c now we are also varying it
lgrid = linspace(0,1,500)';

%Initial value function
vf = ones(1,500)*ut2(css, lss, sigma, chi, phi);

% Placeholders:
vf_place = zeros(1,500);
kpol = zeros(1,500);
lpol = zeros(1,500);
cpol = zeros(1,500);

difference = 90;

%Solving using VFI
for p = 1: Nk
    for i = 1:500
        k = kgrid(i);
        v_temp = -inf;
        for h = 1:500 
            l = lgrid(h);
            c = ((chi*l^(1/phi))/((1-alpha)*k^alpha*l^(-alpha)))^(1/-sigma);
            kplus1 = k^alpha*l^(1-alpha)+(1-delta)*k-c;
                if kplus1 < 0
                    continue
                else
                    vplus1 = interp1(kgrid,vf,kplus1,"linear", "extrap");
                    v = c^(1-sigma)/(1-sigma)-(chi*phi/(1+phi))*l^(1+1/phi)+beta*vplus1;
                         if v > v_temp 
                            v_temp = v;
                            k_new = kplus1;
                            l_new = l;
                            c_new = c;
                         else
                             continue
                         end
                
                end 
       end

    vf_place(i) = v_temp;
    kpol(i) = k_new;
    lpol(i) = l_new;
    cpol(i) = c_new;

    end

difference = max(abs(vf_place-vf));
vf = vf_place;

if difference < 0.0001
    disp("Convergence Achieved")
    break
end 

end

% Plotting
figure; 
plot([1:500], vf);
xlabel('Time');
ylabel('Value');
title ('Value Function, phi = 2');

figure;
plot([1:500], kpol);
xlabel('Time');
ylabel('Capital');
title ('Capital Policy Function, phi = 2');

figure;
plot([1:500], cpol);
xlabel('Time');
ylabel('Consumption');
title ('Consumption Policy Function, phi = 2');

figure;
plot([1:500], lpol);
xlabel('Time');
ylabel('Labor');
title ('Labor Policy Function, phi = 2');


function u = ut2(c, l, sigma, chi, phi)
    if abs(sigma - 1) < eps
        cutility = log(c);
    else
        cutility = (c.^(1 - sigma))/(1 - sigma);
    end

    ldisutility = chi*phi*(l.^(1 + 1 / phi))/(1 + phi);
    % Total utility
    u = cutility - ldisutility;
end
%%
figure;
plot(kgrid, lpol, 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Labor Supply (l)');
title('Labor vs Capital Policy Function, \phi = 2');
grid on;