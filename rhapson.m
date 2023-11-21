function [output] = rhapson(Ri, Ro, lambda, Pi, materialParameters)

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
                
    R1 = Ri+index1*h; % reference position

    r1 = (ri ^ 2 + (1/lambda) * (R1 ^ 2 - Ri ^ 2)) ^ .5; % calculate the radius by mapping from the reference configuration
    F_vector = [(1/lambda) * (R1/r1), r1/R1, lambda];
    F1 = diag(F_vector); % define the the deformation gradient tensor (use the diag function)

    sigma_extra1 = Constitutive_model(F1, materialParameters); % calculate sigma_extra using the constitutive model
    pressure1 = ((1/r1) * (sigma_extra1(2,2) - sigma_extra1(1, 1))) * (R1/(lambda * r1)); % evaluate the function inside the pressure integral, at R1
    force1 = (2*sigma_extra1(3, 3) + sigma_extra1(2, 2) + sigma_extra1(1, 1)) * (R1/(lambda * r1)) * r1; % evaluates the function inside the force integral, at R1

    R2 = R1+h; % next position
    r2 = (ri ^ 2 + (1/lambda) * (R2 ^ 2 - Ri ^ 2)) ^ .5; % calculate the radius by mapping from the reference configuration

    F_vector = [(1/lambda) * (R2/r2), r2/R2, lambda];
    F2 = diag(F_vector); % define the the deformation gradient tensor (use the diag function)
    sigma_extra2 = Constitutive_model(F2, materialParameters); % calculate sigma_extra using the constitutive model
    pressure2 = ((1/r2) * (sigma_extra2(2,2) - sigma_extra2(1, 1))) * (R2/(lambda * r2)); % evaluate the function inside the integral, at R2
    force2 = (2*sigma_extra2(3, 3) + sigma_extra2(2, 2) + sigma_extra2(1, 1)) * (R2/(lambda * r2)) * r2; % evaluates the function inside the force integral, at R1

    % Calculate the pressure integral by solving the Trapezoidal Rule
    local_P = h * (pressure1 + pressure2) / 2;
    P = P + local_P; 

    % Calculate the force integral by solving the Trapezoidal Rule
    local_F = h * (force1 + force2) / 2;
    F = F + local_F;

end

pressure = F - Pi; % ultimately fvalue needs to go to zero (within the defined tolerance)
force = F;
output = [pressure, force];

end