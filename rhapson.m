function [x] = rhapson(x0, g)
%{
x0 is the initial guess
g is the equilibrium function to be minimized that returns the value to be
minimized in the position (1,1) of result
%}

% Iteration Constraints
tol = 0.01;    % tolerance criteria to stop scheme
maxiter = 100;  % maximum number of iterations allowed before stopping

% Newton Raphson Iteration 
% THIS NEEDS TO BE UPDATED, IDK IF THIS IS THE RIGHT WAY OF TESTING 
% IF THE CODE WOULD BREAK
if length([1]) ~= length(x0)  % check if the number of equations match the number of unknowns
    disp('1Error: number of unknowns does not coincide with number of equations')
    return
    
else
    
    % first perform Newton Raphson for the initial guess 
    
    % evaluate the function, f, with the initial guess, x0
    f = g(x0); % call: equilibrium
    
    iter = 0; % initial iteration
    h = 0.01; % step size to perform the derivative (forward difference)    
    
    xinc = x0; % initial guess of the root
    xinc = xinc + h; % current root + derivative step
    f_prime0 = (g(xinc)-f)/h; % compute the derivative of f, between x0 and xinc
    u0 = f/f_prime0; % divide f by its derivative, f/f'
    x = x0 - u0; % update guess for ro: X = X0 - f/f' (Newton-Raphson Method)
    

    % Now, continue iterating until either a solution is found (within the
    % defined tolerance) or you have reached the maximum number of
    % iterations
    
    while iter < maxiter || x < tol % iterate until one of the criteras are met
        
        
        iter = iter+1; % current iteration number        
        f = g(x); % evaluate f at x
        xinc = x; % current root guess
        xinc = xinc + h; % current root + derivative step
        f_prime = (g(xinc)-f)/h; % compute the derivative of f, between x and xinc
        u = f/f_prime; % divide f by its derivative, f/f'
        
        x = x-u; % update guess for ro: X_n = X_(n-1) - f/f' (Newton-Raphson Method)
        
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
g % function evaluated 
x % root, ro
end