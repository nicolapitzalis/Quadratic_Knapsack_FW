function alpha = line_search(problem, d, g, alpha_max)

  den = 2 * d' * problem.Q * d;

  if den <= 1e-16  % d' * Q * d = 0  ==>  f is linear along d
    alpha = alpha_max;    % just take the maximum possible stepsize
  else
    % optimal unbounded stepsize restricted to max feasible step
    alpha = min( [ ( - g' * d ) / den , alpha_max ] );

  end
end