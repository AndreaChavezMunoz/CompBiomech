function [output] = equilibrium(Ri, Ro, lambda, Pi, materialParameters)
%{
Input Parameters:
Ri: Reference inner radius.
Ro: Reference outer radius.
lambda: Deformation parameter.
Pi: Applied pressure.
materialParameters: Structure containing material parameters needed for the constitutive model.

Output:
output: A 1x2 array containing the computed pressure and force.

%}

% calculate the current inner radius (using the reference as a guess)
ri = sqrt(Ro.^2-1./lambda*(Ro^2-Ri^2)); 

a = Ri;   % lower limit of the independent variable a
b = Ro;   % upper limit of the independent variable b
P = 0;    % initial value for the pressure integral
F = 0;    % initial value for the force integral

n = 1000;    % number of spatial steps
h = (b - a)/n;  % spatial step size, based on n and the bounds [a,b]

% for loop going from inner to outer radius
for index1 = 0:n-1 
     
    % First position
    R1 = Ri+index1*h; % reference position

    % calculate the current radius by mapping from the reference configuration
    r1 = (ri ^ 2 + (1/lambda) * (R1 ^ 2 - Ri ^ 2)) ^ .5;
    F_vector = [(1/lambda) * (R1/r1), r1/R1, lambda];

    % define the the deformation gradient tensor (use the diag function)
    F1 = diag(F_vector); 

    % calculate sigma_extra using the constitutive model
    sigma_extra1 = Constitutive_model(F1, materialParameters); 
    % evaluate the function inside the pressure integral, at R1
    pressure1 = ((1/r1) * (sigma_extra1(2,2) - sigma_extra1(1, 1))) * (R1/(lambda * r1)); 
    % evaluates the function inside the force integral, at R1
    force1 = (2*sigma_extra1(3, 3) + sigma_extra1(2, 2) + sigma_extra1(1, 1)) * (R1/(lambda * r1)) * r1; 

    % Second position
    R2 = R1+h; 
    % calculate the radius by mapping from the reference configuration
    r2 = (ri ^ 2 + (1/lambda) * (R2 ^ 2 - Ri ^ 2)) ^ .5; 

    % define the the deformation gradient tensor (use the diag function)
    F_vector = [(1/lambda) * (R2/r2), r2/R2, lambda];
    F2 = diag(F_vector); 
    % calculate sigma_extra using the constitutive model
    sigma_extra2 = Constitutive_model(F2, materialParameters); 
    % evaluate the function inside the integral, at R2
    pressure2 = ((1/r2) * (sigma_extra2(2,2) - sigma_extra2(1, 1))) * (R2/(lambda * r2)); 
    % evaluates the function inside the force integral, at R2
    force2 = (2*sigma_extra2(3, 3) + sigma_extra2(2, 2) + sigma_extra2(1, 1)) * (R2/(lambda * r2)) * r2; 

    % Calculate the pressure integral by solving the Trapezoidal Rule
    local_P = h * (pressure1 + pressure2) / 2;
    P = P + local_P; % Adds pressure to total

    % Calculate the force integral by solving the Trapezoidal Rule
    local_F = h * (force1 + force2) / 2;
    F = F + local_F; % Adds force to total

end

pressure = F - Pi; % ultimately fvalue needs to go to zero (within the defined tolerance)
force = F;
output = [pressure, force];

end