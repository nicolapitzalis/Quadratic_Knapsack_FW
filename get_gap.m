function [bestlb, gap] = get_gap(v, g, y, x, bestlb)
  
  % The `get_gap` function calculates the duality gap for a given solution in the context of a convex optimization problem.
  % The duality gap is a measure of how far the current solution is from optimality. This function also tracks the best lower bound encountered so far.
  %
  % Syntax:
  % [bestlb, gap] = get_gap(v, g, y, x, bestlb)
  %
  % Inputs:
  % - v: (scalar) The current objective function value at the current point.
  % - g: (vector) The gradient of the objective function at the current point.
  % - y: (vector) The solution to the linear subproblem (e.g., Frank-Wolfe direction).
  % - x: (vector) The current point in the optimization process.
  % - bestlb: (scalar) The best lower bound encountered so far.
  %
  % Outputs:
  % - bestlb: (scalar) The updated best lower bound after comparing it with the current computed lower bound.
  % - gap: (scalar) The relative duality gap, defined as the difference between the objective value and the best lower bound, normalized by the objective value or 1 (whichever is greater).

  % compute new lower bound
  lb = v + g' * ( y - x );

  % update if best
  if lb > bestlb
    bestlb = lb;
  end

  % compute the relative gap
  gap = ( v - bestlb ) / max( [ abs( v ) 1 ] );

end