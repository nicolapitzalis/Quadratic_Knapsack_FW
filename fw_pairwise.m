function [ val , varargout ] =  fw_pairwise( BCQP , varargin )

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

if length( varargin ) > 3
   verbose = varargin{ 4 };
   if ~ islogical( verbose ) || ~ isscalar( verbose )
      error( 'Verbose is not a logical scalar' );
   end
else
   verbose = false;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
x0 = LP(BCQP, zeros(n, 1));
x = x0;  % start from the middle of the box

bestlb = - inf;  % best lower bound so far (= none, really)

S = {x0};  % Initial active set
alpha = 1;  % Initial step_size for x0

fprintf( 'Frank-Wolfe method\n');
if verbose
   fprintf( 'iter\tf(x)\t\tlb\t\tgap\n\n' );
end

i = 1;

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
   val = x' * BCQP.Q * x + BCQP.q' * x;
   g = 2 * BCQP.Q * x + BCQP.q;

   % Compute FW direction
   s = LP(BCQP, g);

   % Compute Away direction
   % disp(i)
   
   % iterate over the active set S
   max_val = -inf;
   for j = 1:length(S)
      val_away = g' * S{j};
      if val_away > max_val
         max_val = val_away;
         v_index = j;
      end
   end
   % [~, v_index] = max(g' * cell2mat(S)); % argmax_{v_k in S} <g, v_k> TODO: change
   v = S{v_index};

   % Pair-wise step
   d = s - v;
   index_v = cellfun(@(x_) isequal(x_, v), S); % find index of v in S
   alpha_v = alpha(index_v);
   step_size_max = alpha_v;

   [bestlb, gap] = get_gap(val, g, s, x, bestlb);

   if verbose
      fprintf( '%4d\t%1.8e\t%1.8e\t%1.4e\n' , i , val , bestlb , gap );
   end

   if gap <= eps
      status = 'optimal';
      break;
   end
   
   if i > MaxIter
      status = 'stopped';
      break;
   end

   step_size = line_search(BCQP, d, g, step_size_max);
   x = x + step_size * d;
   i = i + 1;

   % Update the active set S and coefficients alpha
   index_s = cellfun(@(x_) isequal(x_, s), S); % find index of s in S
      
   if all(index_s == 0) % s not in S
      % add s to S
      S = [S, s];
      alpha = [alpha, step_size];
   end
   
   index_v = cellfun(@(x_) isequal(x_, v), S); % find index of v in S
   alpha(index_s) = alpha(index_s) + step_size;
   alpha(index_v) = alpha(index_v) - step_size;

   if alpha(index_v) == 0 % step_size_max - step_size = 0 --> alpha_v = 0 --> remove v from S 
      S(index_v) = [];
      alpha(index_v) = [];
   end
   

end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n' , ...
         i , status , val , gap );

if nargout > 1
   varargout{ 1 } = x;
end

if nargout > 2
   varargout{ 2 } = status;
end

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



