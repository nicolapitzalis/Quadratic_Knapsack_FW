function problem = gen_problem(n, actv, density, max_eig, lmin, lmax, umin, umax, seed)

rng(seed)

% generate l, u- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
problem.l = lmin * ones( n, 1 ) + ( lmax - lmin ) * rand( n, 1 );
problem.u = umin * ones( n , 1 ) + ( umax - umin ) * rand( n , 1 );

% generate Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rc = rand(n, 1) * max_eig;
problem.Q = sprandsym(n, density, rc);

% generate q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
z = zeros( n , 1 );

% outb( i ) = true if z( i ) will be out of the bounds
outb = rand( n , 1 ) <= actv;

% 50/50 chance of being left of lb or right of ub
lr = rand( n , 1 ) <= 0.5;
l = outb & lr;
r = outb & ~ lr;

% a random amount left of the lb (0)
z( l ) = - rand( sum( l ) , 1 ) .* problem.u( l );

% a random amount right of the ub (u)
z( r ) = problem.u( r ) .* ( 1 + rand( sum( r ) , 1 ) );

outb = ~ outb;  % entries that will be inside the bound
% pick at random in [ 0 , u ]
z( outb ) = rand( sum( outb ) , 1 ) .* problem.u( outb );

problem.q = - 0.5 * problem.Q * z;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% generate ax=b - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
problem.a = rand(n, 1);
x0 = problem.l + (problem.u - problem.l) .* rand(size(problem.l));
problem.b = problem.a' * x0;
% r = rand();
% problem.b = r * (problem.a'*problem.l) + (1-r) * (problem.a'*problem.u);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end