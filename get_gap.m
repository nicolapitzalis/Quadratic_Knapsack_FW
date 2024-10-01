function [bestlb, gap] = get_gap(v, g, y, x, bestlb)
    
  % compute new lower bound
  lb = v + g' * ( y - x );

  % update if best
  if lb > bestlb
    bestlb = lb;
  end

  % compute the relative gap
  gap = ( v - bestlb ) / max( [ abs( v ) 1 ] );