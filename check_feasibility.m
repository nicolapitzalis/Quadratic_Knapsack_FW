function is_feasible = check_feasibility(x, l, u, a, b)

    % The `check_feasibility` function checks whether a given point x satisfies both the box constraints and an inequality constraint of the form a'x >= b. 
    % It returns a boolean indicating whether the point is feasible or not.
    % The function checks both conditions and returns `true` only if both are satisfied. If the point is feasible, 
    % it prints a message confirming this. If the point is not feasible, it prints which constraint(s) were violated.
    %
    % Syntax:
    % is_feasible = check_feasibility(x, l, u, a, b)
    %
    % Inputs:
    % - x: (vector) The point to be checked for feasibility.
    % - l: (vector) Lower bounds for the box constraints.
    % - u: (vector) Upper bounds for the box constraints.
    % - a: (vector) Coefficients of the inequality constraint.
    % - b: (scalar) The right-hand side of the inequality constraint a'x >= b.
    %
    % Outputs:
    % - is_feasible: (boolean) Returns true if the point satisfies both the box and inequality constraints, otherwise returns `false`.

    % Check box constraints: l <= x <= u
    box_constraints_satisfied = all(x >= l) && all(x <= u);
    
    % Check inequality constraint: a'x >= b
    inequality_constraint_satisfied = dot(a, x) >= b;
    
    % Determine if both constraints are satisfied
    is_feasible = box_constraints_satisfied && inequality_constraint_satisfied;
    
    % Display result
    if is_feasible
        disp('The point satisfies both the box constraints and the inequality constraint.');
    else
        disp('box constraint: ');
        disp(box_constraints_satisfied);
        disp('inequality constraints: ');
        disp(inequality_constraint_satisfied);
    end
end
