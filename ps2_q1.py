import sympy as sp

### PART 1: Solve the general matrix:

# Define the symbols
a11, a12, a21, a22, lam = sp.symbols('a11 a12 a21 a22 lam')
x, y = sp.symbols('x y')  # For eigenvectors

# Define the matrix A
A = sp.Matrix([[a11, a12],
            [0, a22]])

# Compute the characteristic polynomial det(A - lambda*I) = 0
I = sp.eye(2)  # Identity matrix of size 2x2
char_poly = (A - lam * I).det()

# Solve the characteristic equation for eigenvalues
eigenvalues = sp.solve(sp.Eq(char_poly, 0), lam)

# Solve for eigenvectors for each eigenvalue
eigenvectors = {}
for eigenvalue in eigenvalues:
    null_space = (A - eigenvalue * I).nullspace()
    eigenvectors[eigenvalue] = null_space

# Generate LaTeX strings for the eigenvalues
eigenvalue_latex = [sp.latex(ev) for ev in eigenvalues]

# Generate LaTeX strings for eigenvectors
eigenvector_latex = {sp.latex(ev): [sp.latex(v) for v in vectors] for ev, vectors in eigenvectors.items()}

# Display eigenvalues and eigenvectors in LaTeX
print('Question 1:', '\n',eigenvalue_latex, '\n', eigenvector_latex)

