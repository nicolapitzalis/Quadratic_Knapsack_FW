function [ v , varargout ] =  FWBCQP( BCQP , varargin )

Verbose = true;      % if we log every iteration

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
   t = varargin{ 1 };
   if ~ isreal( t ) || ~ isscalar( t )
      error( 't is not a real scalar' );
   end
   if t < 0 || t > 1
      error( 't is not in [0, 1]' );
   end
else
   t = 0;
end

if length( varargin ) > 1
   eps = varargin{ 2 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-6;
end

if length( varargin ) > 2
   MaxIter = round( varargin{ 3 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
else
   MaxIter = 1000;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
x0 = LP(BCQP, zeros(n, 1));
x = x0;  % start from the middle of the box

bestlb = - inf;  % best lower bound so far (= none, really)

fprintf( 'Frank-Wolfe method\n');
if Verbose
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

   if Verbose
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

   d = y - x;
   alpha = line_search(BCQP, d, g, 1);
   x = x + alpha * d;
   i = i + 1;

end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n' , ...
         i , status , v , gap );

if nargout > 1
   varargout{ 1 } = x;
end

if nargout > 2
   varargout{ 2 } = status;
end

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



