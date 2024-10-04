function diameter = compute_diameter_polytope(problem)
    
    % The compute_diameter_polytope function calculates an approximation of the diameter of the polytope defined by the box constraints of 
    % the optimization problem. The diameter is approximated as the Euclidean norm of the difference between 
    % the upper bound vector u and the lower bound vector l.
    % 
    % Syntax:
    % diameter = compute_diameter_polytope(problem)
    %
    % Inputs:
    % - problem: (struct) A structure representing the optimization problem. It must contain the following fields:
    %   - l: (vector) Lower bounds for the variables.
    %   - u: (vector) Upper bounds for the variables.
    %
    % Outputs:
    % - diameter: (scalar) The Euclidean distance between the lower and upper bounds, which represents the diameter of the polytope.
    %
    % Example Usage:
    % problem.l = [0; 0];
    % problem.u = [1; 1];
    % diameter = compute_diameter_polytope(problem);

    diameter = norm(problem.u - problem.l);
return;