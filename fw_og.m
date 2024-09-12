function [val, varargout] =  fw_og(P, varargin)

   % ----- Inputs -----
   if ~isstruct(P)
      error('P not a struct');
   end

   if ~isfield(P, 'Q') || ~ isfield(P, 'q') || ...
      ~isfield(P, 'u')
      error('P not a well-formed struct');
   end

   if ~ isreal(P.Q) || ~isreal(P.q) || ~isreal(P.u)
      error( 'P not a well-formed struct' );
   end

   n = size(P.q, 1);
   if size(P.q, 2) ~= 1 || ~isequal(size(P.Q), [n n]) || ...
      ~isequal(size(P.u), [n 1])
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
      solution = varargin{3};
      if ~isreal(eps) || ~isscalar(eps)
         error('Solution is not a real scalar');
      end
   else
      solution = NaN;
   end

   if length(varargin) > 3
      verbose = varargin{4};
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
   log_primal_error = zeros(MaxIter + 1, 1);

   x0 = LP(P, zeros(n, 1));
   x = x0;  % start from the middle of the box

   bestlb = - inf;  % best lower bound so far (= none, really)

   fprintf( 'Frank-Wolfe method\n');
   if verbose
      fprintf( 'iter\tf(x)\t\tlb\t\tgap\t\tprimal\t\ttheoretical_bound\n\n' );
   end

   i = 1;
   % --------------------------

   curvature = curvature_function(P);
   gamma = 27/8;

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

      % Compute the primal error
      if ~isnan(solution)
         primal_error = abs(val - solution);
         log_primal_error(i) = primal_error;
      end
      
      theoretical_bound = (2 * gamma * curvature)/(i+2);
      % Print iteration log
      if verbose
         fprintf('%4d\t%1.8e\t%1.8e\t%1.4e\t%1.4e\t%1.4e\n', i, val, bestlb, gap, primal_error, theoretical_bound);
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
      alpha = line_search(P, d, g, 1);
      x = x + alpha * d;
      i = i + 1;

   end
   % ---------------------


   % ----- Outputs -----
   time_elapsed = toc(time_start);
   
   fprintf('stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e, primal = %1.4e\n' , ...
            i, status, val, gap, primal_error);

   log = struct('gaps', log_gap, 'primals', log_primal_error, 'iteration', i, 'time', time_elapsed);

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



