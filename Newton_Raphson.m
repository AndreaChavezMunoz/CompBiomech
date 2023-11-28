function [x] = Newton_Raphson(H, Ri, Ro, lambda, Pi, materialParameters,x0)

%=========================================================================
% Newton Raphson Scheme, X_n = X_(n-1) - f/f'
% 
% H is the handle to the equilibrium equation
% x0 is the initial guess for the outer radius
%
% x is the approximate solution of the Newton Raphson method (in this case,
% ro)
%=========================================================================

%% Iteration Constraints
tol = 0.001;    % tolerance criteria to stop scheme
maxiter = 1000;  % maximum number of iterations allowed before stopping

%% Newton Raphson Iteration
if false  % check if the number of equations match the number of unknowns
    disp('1Error: number of unknowns does not coincide with number of equations')
    return
    
else
    
    % first perform Newton Raphson for the initial guess 
    
    % evaluate the function, f, with the initial guess, x0
    f = H(Ri, Ro, lambda, Pi, materialParameters, x0); % call: equilibrium_r_or_loaded_tutorial
    
    iter = 0; % initial iteration
    h = 0.1; % step size to perform the derivative (forward difference)    
    
    xinc = x0; % initial guess of the root
    xinc = xinc + h; % current root + derivative step
    f_prime0 = ((H(Ri, Ro, lambda, Pi, materialParameters, xinc) - f)) / h; % compute the derivative of f, between x0 and xinc
    u0 = f/f_prime0; % divide f by its derivative, f/f'
    x = x0 - u0; % update guess for ro: X = X0 - f/f' (Newton-Raphson Method)
    

    % Now, continue iterating until either a solution is found (within the
    % defined tolerance) or you have reached the maximum number of
    % iterations
    
    while (iter < maxiter) % iterate until one of the criteras are met
        
        iter = iter+1; % current iteration number        
        f = H(Ri, Ro, lambda, Pi, materialParameters, x); % evaluate f at x
        xinc = x; % current root guess
        xinc = xinc + h; % current root + derivative step
        f_prime = ((H(Ri, Ro, lambda, Pi, materialParameters, xinc) - f)) / h; % compute the derivative of f, between x and xinc
        u = f/f_prime; % divide f by its derivative, f/f'
        x = x - u; % update guess for ro: X_n = X_(n-1) - f/f' (Newton-Raphson Method)
        
    end
end

%% Outputs

% Stop when we meet the max number of iterations
if iter == maxiter 
    disp ('---')
    disp ('Maximum number of iterations reached')
    disp ('---')
    
% or stop when solution is found: the error < desired tolerance     
elseif norm(f) <= tol
    
    disp ('---')
    disp ('Minimum tolereance reached')
    disp ('---')
    
end
iter % number of iterations you performed
f % function evaluated 
x % root, ro
end