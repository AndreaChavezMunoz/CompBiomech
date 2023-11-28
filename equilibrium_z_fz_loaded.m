function [z_fz_output] = equilibrium_z_fz_loaded(Ri, Ro, lambda, materialParameters, x0)

ro = x0(1,1); % outer radius

% calculate the current inner radius (using the reference as a guess)
ri = sqrt(ro.^2-1./lambda*(Ro^2-Ri^2)); 

a = Ri;   % lower limit of the independent variable a
b = Ro;   % upper limit of the independent variable b
T2 = 0;    % initial value for the  integral


n = 1000;    % number of spatial steps
h = (b - a)/n;  % spatial step size, based on n and the bounds [a,b]

% for loop going from inner to outer radius
for index1 = 0:n-1 
                
    R1 = Ri+index1*h; % reference position

    r1 = (ri ^ 2 + (1/lambda) * (R1 ^ 2 - Ri ^ 2)) ^ .5; % calculate the radius by mapping from the reference configuration
    F_vector = [lambda, r1/R1, (1/lambda) * (R1/r1)];
    F1 = diag(F_vector); % define the the deformation gradient tensor (use the diag function)

    sigma_extra1 = constitutive_model(F1, materialParameters); % calculate sigma_extra using the constitutive mode
    f1 = (2*sigma_extra1(3, 3) + sigma_extra1(2, 2) + sigma_extra1(1, 1)) * (R1/(lambda * r1^2)); % evaluates the function inside the force integral, at R1

    R2 = R1+h; % next position
    r2 = (ri ^ 2 + (1/lambda) * (R2 ^ 2 - Ri ^ 2)) ^ .5; % calculate the radius by mapping from the reference configuration

    F_vector = [lambda, r2/R2, (1/lambda) * (R2/r2)];
    F2 = diag(F_vector); % define the the deformation gradient tensor (use the diag function)
    sigma_extra2 = constitutive_model(F2, materialParameters); % calculate sigma_extra using the constitutive model
    f2 = (2*sigma_extra2(3, 3) + sigma_extra2(2, 2) + sigma_extra2(1, 1)) * (R2/(lambda * r2^2)); % evaluates the function inside the force integral, at R2

    % Calculate the force integral by solving the Trapezoidal Rule
    T2 = (((f2+f1) / 2) * h) + T2;

end

z_fz_output = T2 ;

end
