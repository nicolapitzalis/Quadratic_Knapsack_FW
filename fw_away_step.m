function [val, varargout] =  fw_away_step(P, varargin)

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
      MaxIter = round(varargin{2});
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
      verbose = true;
   end

   if length(varargin) > 3
      is_theo_bound = varargin{4};
  else
      is_theo_bound = false;
  end
   % ------------------


   % ----- Initialization -----
   time_start = tic;
   log_gap = zeros(MaxIter + 1, 1);
   log_bound = zeros(MaxIter + 1, 1);

   x0 = LP(P, zeros(n, 1));
   x = x0;  % Initial solution

   bestlb = - inf;  % Best lower bound so far (= none, really)

   % Initialize map for fast lookup
   S_map = containers.Map();  % Create an empty map
   S_keys = {};               % To store string keys representing vectors in S

   % Convert the initial point x0 to a key and add to the map
   key_x0 = mat2str(x0', 16);  % Convert vector to string with high precision
   S_map(key_x0) = 1;          % Map the string to index 1
   S_keys{1} = key_x0;         % Store the key
   S = x0;                     % Initialize active set S with x0
   beta = 1;                   % Initialize beta with coefficient of x0

   fprintf('Frank-Wolfe method\n');
   if verbose
      fprintf('iter\tf(x)\t\tlb\t\tgap\n\n');
   end

   i = 1;
   % --------------------------

   if is_theo_bound
      [~, golden_val, ~, ~] = quadprog(P.Q, P.q, [], [], [], [], P.l, P.u);
      mu = min(eig(P.Q));
      L = max(eig(P.Q));
      delta = computePyramidalWidth(P.a, P.b, P.l, P.u, 1000);
      D = solve_diameter_polytope(P);
      rho = (mu / 4*L) * (delta / D)^2;
   else
      theoretical_bound = NaN;
   end

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

      if is_theo_bound && i > 1
         h = (1-rho) * h;
      elseif is_theo_bound && i == 1
          h = abs(val - golden_val);
      end
      if is_theo_bound
          if (L*(D^2))/2 < h
            theoretical_bound = h + L*(D^2)/2;
         else
            theoretical_bound = D * sqrt(2 * L * h);
          end      
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
      if i > MaxIter
         status = 'stopped';
         break;
      end

      % Line search to find the step size
      step_size = line_search(P, d, g, step_size_max);
      x = x + step_size * d;

      % Update the active set S and the coefficients beta
      if d == d_fw
         % Frank-Wolfe step case
         beta = (1 - step_size) * beta; % Update beta
         
         % [found, index_s] = ismembertol(s', S', 'ByRows', true);
         % % [found, index_s] = ismember(s', S', 'rows');
         % % Update beta_s
         % if ~found  % s not in S
         %    beta_s = step_size;   
         % else
         %    beta_s = beta(index_s);
         %    beta(index_s) = (1 - step_size) * beta_s + step_size;
         % end
         % Convert the new vector s to a key
         key_s = mat2str(s', 16);  % Convert the vector to a string with high precision
         found = isKey(S_map, key_s);
         
         if ~found
            beta_s = step_size;
         else
            index_s = S_map(key_s);
            beta_s = beta(index_s);
            beta(index_s) = (1 - step_size) * beta_s + step_size;
         end

         if step_size < step_size_max && ~found
            S = [S, s];  % Add the new direction to S
            beta = [beta; step_size];  % Add the new coefficient to beta
            S_map(key_s) = size(S, 2);  % Map the string key to the new index in S
            S_keys{end + 1} = key_s;    % Store the key
         end

         if abs(step_size - step_size_max) <= 1e-8 % TODO: this tolerance as global variable or input
            % Reset S to only contain s
            S = s;  % Reset S to only contain s
            beta = beta_s;
            S_map = containers.Map();  % Create an empty map
            S_keys = {};               % To store string keys representing vectors in S
            S_map(key_s) = 1;          % Map the string to index 1
            S_keys{1} = key_s;         % Store the key
         end
      
      else
         % Away step case
         % Convert vector v to a string key
         key_v = mat2str(v', 16);  % Convert the vector to a string with high precision

         beta = (1 + step_size) * beta;
         beta(index_v) = (1 + step_size) * beta_v - step_size;

         if abs(step_size - step_size_max) <= 1e-8
               S(:, index_v) = [];  % Remove column corresponding to v
               beta(index_v) = [];
               % Remove the entry from the map and adjust indices
               remove(S_map, key_v);  % Remove the vector from the map
               S_keys(index_v) = [];  % Remove the key from the key list
   
               % Update the map for all remaining keys in S
               for j = index_v:length(S_keys)
                     S_map(S_keys{j}) = j;  % Adjust the index in the map
               end
         end
      end

      i = i + 1;
   end
   % ---------------------


   % ----- Outputs -----
   time_elapsed = toc(time_start);
   
   fprintf('stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n', i, status, val, gap);

   log = struct('gaps', log_gap, 'bounds', log_bound, 'iteration' , i, 'time', time_elapsed);

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



