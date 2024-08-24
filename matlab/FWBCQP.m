function [ v , varargout ] =  FWBCQP( BCQP , varargin )

%function [ v , x , status ] = FWBCQP( BCQP , t , eps , MaxIter )
%
% Apply the (possibly, stabilized) Franke-Wolfe algorithm with exact line
% search to the convex Box-Constrained Quadratic program
%
%  (P) min { (1/2) x^T * Q * x + q * x : 0 <= x <= u }
%
% encoded in the structure BCQP.
%
% Input:
%
% - BCQP, the structure encoding the BCQP to be solved within its fields:
%
%   = BCQP.Q: n \times n symmetric positive semidefinite real matrix
%
%   = BCQP.q: n \times 1 real vector
%
%   = BCQP.u: n \times 1 real vector > 0
%
% - t (real scalar scalar, optional, default value 0): if the stablized
%   version of the approach is used, then the new point is chosen in the
%   box of relative sixe t around the current point, i.e., the component
%   x[ i ] is allowed to change by not more than plus or minus t * u[ i ].
%   if t = 0, then the non-stabilized version of the algorithm is used.
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the relative gap
%   between the value of the best primal solution (the current one) and the
%   value of the best lower bound obtained so far is less than or equal to
%   eps
%
% - MaxIter (integer scalar, optional, default value 1000): the maximum
%   number of iterations
%
% Output:
%
% - v (real scalar): the best function value found so far (possibly the
%   optimal one)
%
% - x ([ n x 1 ] real column vector, optional): the best solution found so
%   far (possibly the optimal one)
%
% - status (string, optional): a string describing the status of the
%   algorithm at termination, with the following possible values:
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the norm of the gradient at x
%     is less than the required threshold
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the bast solution found so far, but not
%     necessarily the optimal one
%
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 13-12-22
 Version 1.21
 Copyright Antonio Frangioni
 =======================================
%}

Interactive = false;  % if we pause at every iteration

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

x = BCQP.u / 2;  % start from the middle of the box

bestlb = - inf;  % best lower bound so far (= none, really)

fprintf( 'Frank-Wolfe method\n');
if Verbose
   fprintf( 'iter\tf(x)\t\tlb\t\tgap\n\n' );
end

i = 1;

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ Interactive
   tic;
end

while true
    % compute function value and direction = - gradient - - - - - - - - - -

    v = 0.5 * x' * BCQP.Q * x + BCQP.q' * x;
    g = BCQP.Q * x + BCQP.q;

    % solve min { < g , y > : 0 <= y <= u }
    y = linear_subproblem(BCQP, g);
%     y = zeros( n , 1 );
%     ind = g < 0;
%     y( ind ) = BCQP.u( ind );
    
    % compute the lower bound: remember that the first-order approximation
    % is f( x ) + g ( y - x )
    lb = v + g' * ( y - x );
    if lb > bestlb
       bestlb = lb;
    end

    % compute the relative gap
    gap = ( v - bestlb ) / max( [ abs( v ) 1 ] );

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    if Verbose
       fprintf( '%4d\t%1.8e\t%1.8e\t%1.4e\n' , i , v , bestlb , gap );
    end
   
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if gap <= eps      % look, ma! the *true* stopping criterion, for once!
       status = 'optimal';
       break;
    end
    
    if i > MaxIter
       status = 'stopped';
       break;
    end

    % in the stabilized case, restrict y in the box - - - - - - - - - - - -

    if t > 0
       y = max( x - t * BCQP.u , min( x + t * BCQP.u , y ) );
    end
    
    % compute step size - - - - - - - - - - - - - - - - - - - - - - - - - -
    % we are taking direction d = y - x and y is feasible, hence the
    % maximum stepsize is 1

    d = y - x;

    % compute optimal unbounded stepsize:
    % min (1/2) ( x + a d )' * Q * ( x + a d ) + q' * ( x + a d ) =
    %     (1/2) a^2 ( d' * Q * d ) + a d' * ( Q * x + q ) [ + const ]
    %
    % ==> a = - d' * ( Q * x + q ) / d' * Q * d
    %
    den = d' * BCQP.Q * d;
    
    if den <= 1e-16  % d' * Q * d = 0  ==>  f is linear along d
       alpha = 1;    % just take the maximum possible stepsize
    else
       % optimal unbounded stepsize restricted to max feasible step
       alpha = min( [ ( - g' * d ) / den , 1 ] );
    end
    
    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    x = x + alpha * d;
  
    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    i = i + 1;

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n' , ...
         i , status , v , gap );
if ~ Interactive
   toc
end

if nargout > 1
   varargout{ 1 } = x;
end

if nargout > 2
   varargout{ 2 } = status;
end

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



