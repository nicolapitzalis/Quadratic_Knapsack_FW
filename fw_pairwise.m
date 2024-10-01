function [val, varargout] =  fw_pairwise(P, varargin)

    profile on
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
       verbose = false;
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
    
    theoretical_bound = NaN;
   if bound_flag
      % Compute the theoretical bound
      [~, golden_val, ~, ~] = quadprog(P.Q, P.q, [], [], [], [], P.l, P.u);
      mu = min(eig(P.Q));
      L = max(eig(P.Q));
      delta = compute_pyramidal_width(P.a, P.b, P.l, P.u, 1000);
      D = compute_diameter_polytope(P);
      rho = (mu / 4*L) * (delta / D)^2; % Compute rho
   end

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

       % Compute the theoretical bound
      if bound_flag 
         if i > 1
            h = (1-rho) * h;
         else 
            h = abs(val - golden_val);
         end
      
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
       if i > max_iter
          status = 'stopped';
          break;
       end
 
       % Line search to find the step size
       step_size = line_search(P, d, g, step_size_max);
       x = x + step_size * d;
       
       % Convert the new vector s to a key
       key_s = mat2str(s', 16);  % Convert the vector to a string with high precision
 
       % Update the active set S and the coefficients beta
       if ~isKey(S_map, key_s)  % Check if key_s is in the map
          S = [S, s];  % Add the new direction to S
          beta = [beta; step_size];  % Add the new coefficient to beta
          S_map(key_s) = size(S, 2);  % Map the string key to the new index in S
          S_keys{end + 1} = key_s;    % Store the key
       else
          index_s = S_map(key_s);  % Get the index of the existing direction
          beta(index_s) = beta(index_s) + step_size;  % Update the coefficient
       end
 
       % Convert vector v to a string key
       key_v = mat2str(v', 16);  % Convert the vector to a string with high precision
 
       % Find the index of v using the map
       if isKey(S_map, key_v)  % Check if the key exists in the map
          index_v = S_map(key_v);  % Get the index of v in S
          beta(index_v) = beta(index_v) - step_size;  % Update the coefficient
 
          % Check if the coefficient is below the tolerance and should be removed
          if beta(index_v) <= 1e-16
             S(:, index_v) = [];  % Remove the column corresponding to v
             beta(index_v) = [];  % Remove the corresponding beta value
 
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
    
    fprintf('stop: %d iter, status = %s, fbest = %1.8e, gap = %1.4e\n' , i, status, val, gap);
 
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