% Define the symbolic variables
syms a11 a12 a21 a22 lambda

%%% PART 1
fprintf('Question 1: Solve the general matrix\n')

% Define the matrix
A = [a11, a12; a21, a22];

% Compute the characteristic polynomial
charPoly = det(A - lambda * eye(2));

% Solve for eigenvalues
eigenvalues = solve(charPoly, lambda);

% Compute eigenvectors for each eigenvalue
eigenvectors = [];
for i = 1:length(eigenvalues)
    eigenvectors = [eigenvectors, null(A - eigenvalues(i) * eye(2))];
end

% Display eigenvalues in LaTeX
disp('Eigenvalues:');
for i = 1:length(eigenvalues)
    disp(latex(eigenvalues(i)));
end

% Display eigenvectors in LaTeX
disp('Eigenvectors:');
for i = 1:length(eigenvectors)
    disp(latex(eigenvectors(:, i)));
end

%%% PART 2
fprintf('\nQuestion 2: a_21 = 0\n')

% Define the matrix
A = [a11, a12; 0, a22];

% Compute the characteristic polynomial
charPoly = det(A - lambda * eye(2));

% Solve for eigenvalues
eigenvalues = solve(charPoly, lambda);

% Compute eigenvectors for each eigenvalue
eigenvectors = [];
for i = 1:length(eigenvalues)
    eigenvectors = [eigenvectors, null(A - eigenvalues(i) * eye(2))];
end

% Display eigenvalues in LaTeX
disp('Eigenvalues:');
for i = 1:length(eigenvalues)
    disp(latex(eigenvalues(i)));
end

% Display eigenvectors in LaTeX
disp('Eigenvectors:');
for i = 1:length(eigenvectors)
    disp(latex(eigenvectors(:, i)));
end

%%% PART 3
fprintf('\nQuestion 3: a_12 = 0\n')

% Define the matrix
A = [a11, 0; a21, a22];

% Compute the characteristic polynomial
charPoly = det(A - lambda * eye(2));

% Solve for eigenvalues
eigenvalues = solve(charPoly, lambda);

% Compute eigenvectors for each eigenvalue
eigenvectors = [];
for i = 1:length(eigenvalues)
    eigenvectors = [eigenvectors, null(A - eigenvalues(i) * eye(2))];
end

% Display eigenvalues in LaTeX
disp('Eigenvalues:');
for i = 1:length(eigenvalues)
    disp(latex(eigenvalues(i)));
end

% Display eigenvectors in LaTeX
disp('Eigenvectors:');
for i = 1:length(eigenvectors)
    disp(latex(eigenvectors(:, i)));
end