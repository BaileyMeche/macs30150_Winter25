% Zeros of 6th degree chebychev polynomial
chebyRoots = zeros(6,1);

for k = 1:6
    chebyRoots(k) = -cos(((2*k-1)*pi)/12);
end

% Rescale to our range

for k = 1:6
    chebyRoots(k) = (chebyRoots(k)+1)*(100)/2;
end

% Use our equation for R(t;a)


% R(t,a)
R = @(t,a) sum((1:6).*a.*(t.^(0:5)))-0.3*(0.01+sum(a.*(t.^(1:6))))^0.36+0.1*(0.01+sum(a.*t.^(1:6)));

% Need something for guess for solver
aguess = [0.6,0.5,0.4,0.3,0.02,0.001];

% Only want to solve for and vary a
fa = @(a) real(f(a,R,chebyRoots));
as = fsolve(fa, aguess);

% Ok so now to get the approximate capital solution

k_hat = @(t) real(0.01+sum(as.*(t.^(1:6))))

%fplot(k_hat,[0,100])
%title('Approximate Solution')

%And here's the analytical solution

k_anal = @(t) real(((0.01^0.64-3)*exp(-0.064*t)+3)^(1/0.64))


%fplot(k_anal,[0,100])
%fplot(k_hat,[0,100])
%title('Analytical Solution')

figure;
fplot(k_anal, [1,100], 'r', 'DisplayName', 'Analytical Solution');
hold on;
fplot(k_hat, [1,100], '--b', 'DisplayName', 'Approximate Solution');
title('Solow');
ylabel('k');
xlabel('Time');
legend show;
grid on;

% Construct set of functions to solve
function fun = f(a,g,roots)
    fun = zeros(6,1);
    for num = 1:6
        fun(num) = g(roots(num),a);
    end
end