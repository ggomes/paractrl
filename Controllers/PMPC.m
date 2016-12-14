%% These functions represents a parameterized MPC controller


function [p_CMPC, fval, exitflag, output] = PMPC ( Np, Nc, x, P_start, dist, A, b, Aeq, beq, lb, ub, xmax)

    options = optimoptions('Algorithm','sqp');
    
    [P_opt, fval, exitflag, output] = fmincon(@(P) ...
                                        objfun_CMPC (Np, x, P, dist), ... 
                                        P_start, ...
                                        A, b, Aeq, beq, lb, ub, ...
                                        nonlcon_PMPC (Np, Nc, x, P, dist, xmax), ...
                                        options);
                             
                                    
    if ( exitflag == 0 )
        % maximum number of iterations exceeded
        [P_opt, fval, exitflag, output] = fmincon(@(P) ...
                                        objfun_PMPC (Np, x, P, dist), ... 
                                        P_opt, ...
                                        A, b, Aeq, beq, lb, ub, ...
                                        nonlcon_PMPC (Np, Nc, x, P, dist, xmax), ...
                                        options);
    end;
        
    if ( exitflag >= 0 )                                 
        p_CMPC = P_opt;
    else
        p_CMPC = [];   % failure
    end;
    
    % later on we are going to add a time budget; in case the controller
    % finds a local optimum within the budget, it may use other random
    % starting points to search the region for more possible local optima

end


%% Defining the objective function for CMPC

% the objective function should be defined across the prediction horizon

function J = objfun_PMPC (Np, x, P, dist)

    J = 0;
    
    for i = 1 : Np
        if ( i <= Nc )
           p = P(i, :);
        else
           p = P(Nc,:);
        end;
        x = PMPC_model (x, p, dist(i,:));
        J = J + fP (x, p);               
    end

end


function eval_f = fP (x, p)

    % this function should be modified later
    eval_f = abs(x) + abs (p);

end


%% Defining the constraints for CMPC

function [c , ceq] = nonlcon_PMPC (Np, Nc, x, P, dist, xmax)

    n = size(1, xmax);
    c = NaN(n, Np);
    ceq = [];
    
    for i = 1 : Np
        if ( i <= Nc )
           p = P(i, :);
        else
           p = P(Nc, :);
        end;
        x = CMPC_model (x, p, dist(i,:));
        c(:, i) = max(x - xmax(:, i), 0);
    end
    
 end

