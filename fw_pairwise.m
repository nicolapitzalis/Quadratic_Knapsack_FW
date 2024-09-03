function [val, varargout] =  fw_pairwise(BCQP, varargin)

% ----- Inputs -----
if ~ isstruct( BCQP )
   error( 'BCQP not a struct' );
end

if ~ isfield( BCQP , 'Q' ) || ~ isfield( BCQP , 'q' ) || ...
   ~ isfield( BCQP , 'u' )
   error( 'BCQP not a well-formed struct' );
end

if ~ isreal( BCQP.Q ) || ~ isreal( BCQP.q ) || ~ isreal( BCQP.u )
   error( 'BCQP not a well-formed struct' );
end

n = size( BCQP.q , 1 );
if size( BCQP.q , 2 ) ~= 1 || ~ isequal( size( BCQP.Q ) , [ n n ] ) || ...
   ~ isequal( size( BCQP.u ) , [ n 1 ] )
   error( 'BCQP not a well-formed struct' );
end

if ~isempty( varargin )
   eps = varargin{ 1 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-6;
end

if length( varargin ) > 1
   MaxIter = round( varargin{ 2 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
else
   MaxIter = 1000;
end

if length( varargin ) > 2
   verbose = varargin{ 3 };
   if ~ islogical( verbose ) || ~ isscalar( verbose )
      error( 'Verbose is not a logical scalar' );
   end
else
   verbose = false;
end
% ------------------


% ----- Initialization -----
time_start = tic;
log_gap = zeros(MaxIter + 1, 1);

x0 = LP(BCQP, zeros(n, 1));
x = x0;  % Initial solution

bestlb = - inf;  % Best lower bound so far (= none, really)

S = x0;  % Initialize active set
beta = 1;  % Initial coefficient of x0

fprintf('Frank-Wolfe method\n');
if verbose
    fprintf('iter\tf(x)\t\tlb\t\tgap\n\n');
end

i = 1;
% --------------------------


% ----- Main loop -----
while true
   
   % Compute the function value and gradient
   val = x' * BCQP.Q * x + BCQP.q' * x;
   g = 2 * BCQP.Q * x + BCQP.q;

   % Compute FW direction argmin_s {<g, s>: as >= b, l <= a <= u}
   s = LP(BCQP, g);

   % Compute Away direction argmax_v {<g, v>, v in S}
   [~, index_v] = max(g' * S); 
   v = S(:, index_v);

   % Pairwise direction
   d = s - v;
   step_size_max = beta(index_v);

   % Compute the gap
   [bestlb, gap] = get_gap(val, g, s, x, bestlb);
   log_gap(i) = gap;

   % Print iteration log
   if verbose
      fprintf('%4d\t%1.8e\t%1.8e\t%1.4e\n', i, val, bestlb, gap);
   end

   % Check stopping criteria
   if gap <= eps
      status = 'optimal';
      break;
   end
   if i > MaxIter
      status = 'stopped';
      break;
   end

   % Line search to find the step size
   step_size = line_search(BCQP, d, g, step_size_max);
   x = x + step_size * d;
   i = i + 1;

   % Update the active set S and the coefficients beta
   [found, index_s] = ismember(s', S', 'rows');

   if ~found
      S = [S, s];
      beta = [beta; step_size];
   else
      beta(index_s) = beta(index_s) + step_size;
   end

   [~, index_v] = ismember(v', S', 'rows');
   beta(index_v) = beta(index_v) - step_size;

   if beta(index_v) == 0 % step_size_max - step_size = 0 --> alpha_v = 0 --> remove v from S 
      S(:, index_v) = []; 
      beta(index_v) = [];
   end

end

time_elapsed = toc(time_start);
% ---------------------


% ----- Outputs -----
fprintf( 'stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n' , ...
         i , status , val , gap );

log = struct('gaps', log_gap, 'solution', x, 'iteration', i, 'time', time_elapsed);

if nargout > 1
   varargout{ 1 } = log;
end

if nargout > 2
   varargout{ 2 } = status;
end
% ---------------------

end 



