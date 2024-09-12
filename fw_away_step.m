function [val, varargout] =  fw_away_step(P, varargin)

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
      MaxIter = round(varargin{2});
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
   x = x0;  % Initial solution

   bestlb = - inf;  % Best lower bound so far (= none, really)

   S = x0;  % Initialize active set
   beta = 1;  % Initial coefficient of x0

   fprintf('Frank-Wolfe method\n');
   if verbose
      fprintf('iter\tf(x)\t\tlb\t\tgap\t\tprimal\n\n');
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
      d_fw = s - x;

      % Compute Away direction argmax_v {<g, v>, v in S}
      [~, index_v] = max(g' * S); 
      v = S(:, index_v);
      d_away = x - v;

      % Choose between Frank-Wolfe and Away direction
      if -g' * d_fw >= -g' * d_away || size(S, 2) == 1
         d = d_fw;
         step_size_max = 1;
      else
         d = d_away;
         beta_v = beta(index_v);
         step_size_max = beta_v / (1 - beta_v);
      end

      % Compute the gap a.k.a. the dual error
      [bestlb, gap] = get_gap(val, g, s, x, bestlb);
      log_gap(i) = gap;

      % Compute the primal error
      if ~isnan(solution)
         primal_error = abs(val - solution);
         log_primal_error(i) = primal_error;
      end

      % Print iteration log
      if verbose
         fprintf('%4d\t%1.8e\t%1.8e\t%1.4e\t%1.4e\n', i, val, bestlb, gap, primal_error);
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
      i = i + 1;

      % Update the active set S and the coefficients beta
      if d == d_fw
         % Frank-Wolfe step case
         beta = (1 - step_size) * beta; % Update beta
         
         [found, index_s] = ismember(s', S', 'rows');
         % Update beta_s
         if ~found  % s not in S
            beta_s = step_size;   
         else
            beta_s = beta(index_s);
            beta(index_s) = (1 - step_size) * beta_s + step_size;
         end

         if step_size < step_size_max && ~found
            S = [S, s];  % Add s to S % consider preallocating 1000 iterations
            beta = [beta; beta_s]; % Add beta_s to beta
         end
         if step_size == step_size_max
            S = s;  % Reset S to only contain s
            beta = beta_s;
         end
      
         else
            % Away step case
            beta = (1 + step_size) * beta;
            beta(index_v) = (1 + step_size) * beta_v - step_size;

            if step_size == step_size_max
                  S(:, index_v) = [];  % Remove column corresponding to v
                  beta(index_v) = [];
            end
      end

   end
   % ---------------------


   % ----- Outputs -----
   time_elapsed = toc(time_start);

   fprintf( 'stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e, primal = %1.4e\n' , ...
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



