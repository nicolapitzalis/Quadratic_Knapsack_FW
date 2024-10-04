function [val, varargout] =  fw_og(P, varargin)

   % The fw_og function implements the Original Frank-Wolfe (FW) algorithm for solving a convex quadratic 
   % programming problem in the form:
   % min_x  x'Qx + q'x
   % subject to:
   % l <= x <= u   and   a'x >= b
   % The algorithm iteratively minimizes a quadratic objective function over a feasible set 
   % constrained by lower and upper bounds. The function tracks the duality gap and provides options for detailed 
   % logging and theoretical bound computation.

   % Syntax:
   % [val, x, status, log] = fw_og(P, epsilon, max_iter, verbose, bound_flag)

   % Inputs:
   % - P: (struct) A structure representing the problem instance with the following fields:
   %   - Q: (matrix) The quadratic matrix Q (should be symmetric and positive semi-definite).
   %   - q: (vector) The linear term in the objective function.
   %   - l: (vector) Lower bounds for the variables.
   %   - u: (vector) Upper bounds for the variables.
   % - epsilon: (float, optional) Stopping criterion for the duality gap. If the duality gap is smaller than epsilon, the algorithm terminates. Default is 1 x 10^{-6}.
   % - max_iter: (integer, optional) Maximum number of iterations allowed. Default is 1000.
   % - verbose: (boolean, optional) Whether to display detailed information during each iteration. Default is `true`.
   % - bound_flag: (boolean, optional) If `true`, the function computes theoretical bounds for the algorithm’s convergence. Default is `false`.

   % Outputs:
   % - val: (float) The objective function value at the final solution.
   % - x: (vector, optional) The solution vector at the final iteration.
   % - status: (string, optional) The termination status of the algorithm. It can be either `'optimal'` if the solution meets the stopping criterion, or `'stopped'` if the maximum number of iterations is reached.
   % - log: (struct, optional) A structure containing logs of the optimization process with the following fields:
   %   - gaps: (vector) The dual gap at each iteration.
   %   - bounds: (vector) The computed theoretical bounds (if `bound_flag` is `true`).
   %   - iteration: (integer) The number of iterations performed.
   %   - time: (float) Time taken to complete the optimization process.
   %
   % Example Usage:
   % P.Q = sprandsym(10, 0.3);  
   % P.q = rand(10, 1);        
   % P.l = zeros(10, 1);       
   % P.u = ones(10, 1);   
   % [val, x, status, log] = fw_og(P, 1e-6, 500, true, true);
   %
   % Notes:
   % - The function uses LP to solve linear subproblems, which minimizes a linear function over the feasible region 
   % defined by the contraints.
   % - The duality gap provides a measure of how close the current solution is to optimality.
   % - If bound_flag is set, theoretical convergence bounds are computed at each iteration based on the curvature of the problem.

   % See Also:
   % - `LP`: Solves the linear subproblem within the Frank-Wolfe method.
   % - `line_search`: Determines the optimal step size for updating the solution.
   % - `curvature_function`: Computes the curvature of the problem for theoretical bound computation.

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
      max_iter = round( varargin{2} );
      if ~ isscalar(max_iter)
         error('max_iter is not an integer scalar');
      end
   else
      max_iter = 1000;
   end

   if length(varargin) > 2
      verbose = varargin{3};
      if ~islogical(verbose) || ~isscalar(verbose)
         error('Verbose is not a logical scalar');
      end
   else
      verbose = true;
   end
   
   if length(varargin) > 3
       bound_flag = varargin{4};
   else
       bound_flag = false;
   end
   % ------------------

   % ----- Initialization -----
   time_start = tic;
   log_gap = zeros(max_iter + 1, 1);
   log_bound = zeros(max_iter + 1, 1);

   x0 = LP(P, zeros(n, 1));
   x = x0;  % start from the middle of the box

   bestlb = - inf;  % best lower bound so far (= none, really)

   theoretical_bound = NaN;
   if bound_flag
      curvature = curvature_function(P);
      gamma = 27/8;
   end

   fprintf( 'Frank-Wolfe method\n');
   if verbose
      fprintf( 'iter\tf(x)\t\tlb\t\tgap\t\ttheoretical_bound\n\n' );
   end

   i = 1;
   % --------------------------
   

   % ----- Main loop -----
   while true
      
      % Compute the function value and gradient
      val = x' * P.Q * x + P.q' * x;
      g = 2 * P.Q * x + P.q;

      % argmin_s {<g, s>: as >= b, l <= a <= u}
      y = LP(P, g);

      % Direction
      d = y - x;

      % Compute the gap
      [bestlb, gap] = get_gap(val, g, y, x, bestlb);
      log_gap(i) = gap;
      
      if bound_flag
        theoretical_bound = (2 * gamma * curvature)/(i+2);
      end
      log_bound(i) = theoretical_bound;

      % Print iteration log
      if verbose
         fprintf('%4d\t%1.8e\t%1.8e\t%1.4e\t%1.4e\n', i, val, bestlb, gap, theoretical_bound);
      end

      % Check stopping criteria
      if gap <= eps
         status = 'optimal';
         break;
      end
      if i > max_iter
         status = 'stopped';
         break;
      end

      % Line search to find the step size
      alpha = line_search(P, d, g, 1);
      x = x + alpha * d;
      i = i + 1;

   end
   % ---------------------


   % ----- Outputs -----
   time_elapsed = toc(time_start);
   
   fprintf('stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n', i, status, val, gap);

   log = struct('gaps', log_gap, 'bounds', log_bound, 'iteration', i, 'time', time_elapsed);

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



