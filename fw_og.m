function [ v , varargout ] =  fw_og( BCQP , varargin )

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

time_start = tic;
log_gap = zeros(MaxIter + 1, 1);

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
x0 = LP(BCQP, zeros(n, 1));
x = x0;  % start from the middle of the box
check_feasibility(x, BCQP.l, BCQP.u, BCQP.a, BCQP.b);
bestlb = - inf;  % best lower bound so far (= none, really)

fprintf( 'Frank-Wolfe method\n');
if verbose
   fprintf( 'iter\tf(x)\t\tlb\t\tgap\n\n' );
end

i = 1;

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
   v = x' * BCQP.Q * x + BCQP.q' * x;
   g = 2 * BCQP.Q * x + BCQP.q;


   y = LP(BCQP, g);


   [bestlb, gap] = get_gap(v, g, y, x, bestlb);
   
   log_gap(i) = gap;
   
   if verbose
      fprintf( '%4d\t%1.8e\t%1.8e\t%1.4e\n' , i , v , bestlb , gap );
   end

   if gap <= eps
      status = 'optimal';
      break;
   end
   
   if i > MaxIter
      status = 'stopped';
      break;
   end
   % fixed step size
   
   d = y - x;
   alpha = line_search(BCQP, d, g, 1);
   % alpha = 2 / (i + 2);
   x = x + alpha * d;
   i = i + 1;

end

time_elapsed = toc(time_start);

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n' , ...
         i , status , v , gap );

log = struct('gaps', log_gap, 'solution', x, 'iteration', i, 'time', time_elapsed);

if nargout > 1
   varargout{ 1 } = log;
end

if nargout > 2
   varargout{ 2 } = status;
end

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



