function pyramidalWidth = computePyramidalWidth(a, b, l, u, num_samples)
    % INPUTS:
    % A, b: Linear inequality constraints (A * x <= b)
    % l, u: Box constraints (l <= x <= u)
    % a: Coefficients of the linear inequality (a'x >= b)
    % num_samples: Number of random directions to sample
    %
    % OUTPUT:
    % pyramidalWidth: The computed pyramidal width

    n = length(l);  % Dimension of the problem
    pyramidalWidth = Inf;  % Initialize to a large number

    % Sample random directions
    directions = randn(num_samples, n);  % Generate random directions
    directions = directions ./ vecnorm(directions, 2, 2);  % Normalize directions

    % Loop through each sampled direction
    for i = 1:num_samples
        r = directions(i, :);  % Current direction

        % Solve LP to maximize the dot product in the direction of r
        [x_max, fval_max] = linprog(-r, -a', -b, [], [], l, u);  % Maximize r'x
        if isempty(x_max) || a' * x_max < b
            continue;  % Skip if no feasible solution or doesn't satisfy a'x >= b
        end

        % Solve LP to minimize the dot product in the direction of r
        [x_min, fval_min] = linprog(r, -a', -b, [], [], l, u);  % Minimize r'x
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
    else
        fprintf('Pyramidal width: %.4f\n', pyramidalWidth);
    end
end
