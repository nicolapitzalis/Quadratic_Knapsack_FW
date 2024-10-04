function pyramidalWidth = compute_pyramidal_width(a, b, l, u, num_samples)

    % The compute_pyramidal_width function approximates the pyramidal width of the polytope defined by the box 
    % constraints and the inequality constraint a'x >= b.
    % The pyramidal width is an important measure in optimization that helps describe the geometry of the feasible region. 
    % It is calculated by sampling random directions and computing the width of the polytope along these directions.
    % This function provides an approximation of the pyramidal width by solving linear programs (LPs) in random directions.
    %
    % Syntax:
    % pyramidalWidth = compute_pyramidal_width(a, b, l, u, num_samples)
    %
    % Inputs:
    % - a: (vector) Coefficients of the inequality constraint.
    % - b: (scalar) The right-hand side of the inequality constraint a'x >= b.
    % - l: (vector) Lower bounds for the variables.
    % - u: (vector) Upper bounds for the variables.
    % - num_samples: (integer) The number of random directions to sample for the approximation.
    %
    % Outputs:
    % - pyramidalWidth: (scalar) The approximated pyramidal width of the polytope. If no valid width is found, an error is raised.
    %
    % Example Usage:
    % a = [0.3; 0.5];
    % b = 1;
    % l = [0; 0];
    % u = [1; 1];
    % num_samples = 1000;
    % pyramidalWidth = compute_pyramidal_width(a, b, l, u, num_samples);
    %
    % Notes:
    % - The pyramidal width is computed by solving two LPs for each sampled direction: one to maximize and one to minimize the dot product in the given direction.
    % - This function returns the minimum width across all sampled directions, providing an approximation based on the chosen number of samples.

    n = length(l);  % Dimension of the problem
    pyramidalWidth = Inf;  % Initialize to a large number

    % Sample random directions
    directions = randn(num_samples, n);  % Generate random directions
    directions = directions ./ vecnorm(directions, 2, 2);  % Normalize directions

    % Loop through each sampled direction
    for i = 1:num_samples
        r = directions(i, :);  % Current direction
        
        options = optimoptions('linprog', 'Display', 'none');

        % Solve LP to maximize the dot product in the direction of r
        [x_max, fval_max] = linprog(-r, -a', -b, [], [], l, u, options);  % Maximize r'x
        if isempty(x_max) || a' * x_max < b
            continue;  % Skip if no feasible solution or doesn't satisfy a'x >= b
        end

        % Solve LP to minimize the dot product in the direction of r
        [x_min, fval_min] = linprog(r, -a', -b, [], [], l, u, options);  % Minimize r'x
        if isempty(x_min) || a' * x_min < b
            continue;  % Skip if no feasible solution or doesn't satisfy a'x >= b
        end

        % Compute the directional width
        directional_width = fval_max - fval_min;

        % Update the minimum pyramidal directional width
        if directional_width < pyramidalWidth
            pyramidalWidth = directional_width;
        end
    end

    if pyramidalWidth == Inf
        error('No valid pyramidal width found.');
    end
end
