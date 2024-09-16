function [val, varargout] =  fw_pairwise(P, varargin)

   % profile on
   % ----- Inputs -----
   if ~isstruct(P)
      error('P not a struct');
   end

   if ~isfield(P, 'Q') || ~ isfield(P, 'q') || ...
      ~isfield(P, 'u') || ~isfield(P, 'l')
      error('P not a well-formed struct');
   end

   if ~ isreal(P.Q) || ~isreal(P.q) || ~isreal(P.u) || ~isreal(P.l)
      error( 'P not a well-formed struct' );
   end

   n = size(P.q, 1);
   if size(P.q, 2) ~= 1 || ~isequal(size(P.Q), [n n]) || ...
      ~isequal(size(P.u), [n 1]) || ~isequal(size(P.l), [n 1])
      error('P not a well-formed struct');
   end

   if ~isempty(varargin)
      eps = varargin{1};
      if ~isreal(eps) || ~isscalar(eps)
         error('eps is not a real scalar');
      end
   else
      eps = 1e-6;
   end

   if length(varargin) > 1
      MaxIter = round( varargin{2} );
      if ~ isscalar(MaxIter)
         error('MaxIter is not an integer scalar');
      end
   else
      MaxIter = 1000;
   end

   if length(varargin) > 2
      verbose = varargin{3};
      if ~islogical(verbose) || ~isscalar(verbose)
         error('Verbose is not a logical scalar');
      end
   else
      verbose = false;
   end
   % ------------------


   % ----- Initialization -----
   time_start = tic;
   log_gap = zeros(MaxIter + 1, 1);

   x0 = LP(P, zeros(n, 1));
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
      val = x' * P.Q * x + P.q' * x;
      g = 2 * P.Q * x + P.q;

      % Compute FW direction argmin_s {<g, s>: as >= b, l <= a <= u}
      s = LP(P, g);

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
      step_size = line_search(P, d, g, step_size_max);
      x = x + step_size * d;
      

      % Update the active set S and the coefficients beta
      [found, index_s] = ismembertol(s', S', 'ByRows', true); 

      if ~found
         S = [S, s];
         beta = [beta; step_size];
      else
         beta(index_s) = beta(index_s) + step_size;
      end

      [~, index_v] = ismembertol(v', S', 'ByRows', true);
      beta(index_v) = beta(index_v) - step_size;

      if beta(index_v) <= 1e-8 % step_size_max - step_size = 0 --> alpha_v = 0 --> remove v from S 
         S(:, index_v) = []; 
         beta(index_v) = [];
      end
      
      i = i + 1;
   end
   % ---------------------


   % ----- Outputs -----
   time_elapsed = toc(time_start);
   
   fprintf('stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n' , i, status, val, gap);

   log = struct('gaps', log_gap, 'iteration', i, 'time', time_elapsed);

   if nargout > 1
      varargout{1} = x;
   end

   if nargout > 2
      varargout{2} = status;
   end

   if nargout > 3
      varargout{3} = log;
   end
   % ---------------------

end 



