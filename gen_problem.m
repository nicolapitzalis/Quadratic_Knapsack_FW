function problem = gen_problem(n, actv, density, cond_num, lmin, lmax, umin, umax, seed)

rng(seed)

% generate l, u- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
problem.l = lmin * ones( n, 1 ) + ( lmax - lmin ) * rand( n, 1 );
problem.u = umin * ones( n , 1 ) + ( umax - umin ) * rand( n , 1 );

% generate Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Max eigenvalue is 1, control the smallest eigenvalue based on the condition number
max_eig = 1;  % Fix max eigenvalue to 1 for simplicity
min_eig = max_eig / cond_num;  % Compute smallest eigenvalue based on condition number
eigenvalues = linspace(min_eig, max_eig, n);  % Linearly spaced eigenvalues

% Create a random orthogonal matrix
[U, ~] = qr(randn(n)); 

% Construct Q with controlled eigenvalues: Q = U * D * U^T
D = diag(eigenvalues);
problem.Q = U * D * U';

% Apply sparsity pattern (if density < 1), keep dense if density = 1
% Apply sparsity pattern (if density < 1), keep dense if density = 1
if density < 1
    % Create a random mask with the given density
    mask = rand(n) < density;
    % Ensure symmetry of the mask
    mask = triu(mask, 0) + triu(mask, 1)';
    % Apply mask to make problem.Q sparse while keeping symmetry
    problem.Q = problem.Q .* mask;
end

% generate q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
z = zeros(n, 1);

% outb( i ) = true if z( i ) will be out of the bounds
outb = rand(n, 1) <= actv;

% 50/50 chance of being left of lb or right of ub
lr = rand(n, 1) <= 0.5;
l = outb & lr;
r = outb & ~lr;

% a random amount left of the lb (0)
z(l) = problem.l(l) .* (rand(sum(l), 1) - 1);

% a random amount right of the ub (u)
z(r) = problem.u(r) .* (1 + rand(sum(r), 1));

outb = ~outb;  % entries that will be inside the bound
% pick at random in [0, u]
z(outb) = rand(sum(outb), 1) .* problem.u(outb);

problem.q = -2 * problem.Q * z;

% generate ax=b - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
problem.a = rand(n, 1);
x0 = problem.l + (problem.u - problem.l) .* rand(size(problem.l));
problem.b = problem.a' * x0;

end
