function [fvalue] = equilibrium_r_or_loaded(Ri, Ro, lambda, Pi, materialParameters, x0)

%=========================================================================
% RADIAL EQUILIBRIUM - Local
%
% This function calculates: 
%    Pi - int_{ri}^{ro} (tqq-trr)/r dr = 0
%
% x0 is the current solution of the Newton-Raphson method
%=========================================================================

ro = x0(1,1); % outer radius current

ri = sqrt(ro.^2 - 1./lambda*(Ro^2-Ri^2)); % calculate the inner radius

a = Ri;   % lower limit of the independent variable a
b = Ro;   % upper limit of the independent variable b
T1 = 0;    % initial value for the integral (T is the result of the integral)
          
n = 20;    % number of spatial steps
h = (b-a)/n;  % spatial step size, based on n and the bounds [a,b]

for index1 = 0:n-1 % for loop going from inner to outer radius
                
    R1 = Ri + index1*h; % reference position
    r1 = sqrt((ri .^ 2) - ((1 / lambda)*(R1^2 - Ri^2))); % calculate the radius by mapping from the reference configuration
    F1 = diag([R1/(lambda*r1) r1/R1 lambda]); % define the the deformation gradient tensor (use the diag function)
    sigma_extra1 = constitutive_model(F1, materialParameters); % calculate sigma_extra using the constitutive model
    f1 = R1/(lambda * r1^2) .* (sigma_extra1(2, 2) - sigma_extra1(1, 1)); % evaluate the function inside the integral, at R1
    
    R2 = R1 + h; % next position
    r2 = sqrt((ri .^ 2) - ((1 / lambda)*(R2^2 - Ri^2))); % calculate the radius by mapping from the reference configuration
    F2 = diag([R1/(lambda*r1) r1/R1 lambda]); % define the the deformation gradient tensor (use the diag function)
    sigma_extra2 = constitutive_model(F2, materialParameters); % calculate sigma_extra using the constitutive model
    f2 = R2/(lambda * r2^2) .* (sigma_extra2(2, 2) - sigma_extra2(1, 1)); % evaluate the function inside the integral, at R2
    
    % calculate the intergral by solving the Trapezoidal Rule
    T1 = (((f2+f1) / 2) * h) + T1; % remember to add the previous T_(index-1);

end

fvalue = T1 - Pi; % ultimately fvalue needs to go to zero (within the defined tolerance)

end
