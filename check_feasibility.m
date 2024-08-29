function is_feasible = check_feasibility(x, l, u, a, b)
    % x: the point to check
    % l: the lower bounds vector
    % u: the upper bounds vector
    % a: the coefficients vector for the inequality constraint
    % b: the scalar bound for the inequality constraint a'x >= b

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
