
function [stress_rr, stress_theta, stress_zz, r_position] = lagrange_stress(Ri, Ro, lambda, Pi, materialParameters,x0)

%=========================================================================
% RADIAL EQUILIBRIUM - Local
%
% This function calculates: 
%    Pi - int_{ri}^{ro} (tqq-trr)/r dr = 0
%
% x0 is the current solution of the Newton-Raphson method
%=========================================================================

ro = x0(1,1); % Current radius
%ro = Ro; %Why?
ri = sqrt(ro.^2 - 1./lambda*(Ro^2-Ri^2)); % calculate the inner radius

a = Ri;   % lower limit of the independent variable a
b = Ro;   % upper limit of the independent variable b
T1 = 0;    % initial value for the integral (T is the result of the integral)
          
n = 20;    % number of spatial steps
h = (b-a)/n;  % spatial step size, based on n and the bounds [a,b]

% Vector containing lagrange multipliers and stresses
stress_rr = zeros(1, n);
stress_zz = zeros(1, n);
stress_theta = zeros(1, n);
r_position = zeros(1, n);

% Obtains material parameters
c1 = materialParameters(1);
c2 = materialParameters(2);
c3 = materialParameters(3);
c4 = materialParameters(4);
c5 = materialParameters(5);
c6 = materialParameters(6);
c = materialParameters(7);

for index1 = 0:n-1 % for loop going from inner to outer radius
                
    R1 = Ri + index1*h; % reference position
    r1 = sqrt((ri .^ 2) - ((1 / lambda)*(R1^2 - Ri^2))); % calculate the current radius by mapping from the reference configuration
    F1 = diag([lambda r1/R1 R1/(lambda*r1)]); % define the the deformation gradient tensor (use the diag function)
    sigma_extra1 = constitutive_model(F1, materialParameters); % calculate sigma_extra using the constitutive model
    f1 = R1/(lambda * r1^2) .* (sigma_extra1(2, 2) - sigma_extra1(1, 1)); % evaluate the function inside the integral, at R1
    
    R2 = R1 + h; % next position
    r2 = sqrt((ri .^ 2) - ((1 / lambda)*(R2^2 - Ri^2))); % calculate the radius by mapping from the reference configuration
    F2 = diag([lambda r2/R2 R2/(lambda*r2)]); % define the the deformation gradient tensor (use the diag function)
    sigma_extra2 = constitutive_model(F2, materialParameters); % calculate sigma_extra using the constitutive model
    f2 = R2/(lambda * r2^2) .* (sigma_extra2(2, 2) - sigma_extra2(1, 1)); % evaluate the function inside the integral, at R2
    
    % calculate the intergral by solving the Trapezoidal Rule
    T1 = (((f2+f1) / 2) * h) + T1; % remember to add the previous T_(index-1);

    % Calculates cauchy and strain components
    C = F1*F1;
    E = 0.5 * (C - eye(3,3));
    Err = E(1,1);
    Etheta =  E(2,2);
    Ezz = E(3,3);

    % Tries to identify the lagrange multiplier
    Q = c1 * Err^2 + c2 * Etheta^2 + c3 * Ezz^2 + 2*c4*Err*Etheta + 2*c5*Etheta*Ezz + 2*c6*Err*Ezz;

    stress_rr(index1 + 1) = C(1,1) * 0.5 * c * exp(Q) * (2*c1*Err + 2*c4*Etheta + 2*c6*Ezz) + Pi - T1;
    stress_theta(index1 + 1) = C(2,2) * 0.5 * c * exp(Q) * (2*c2*Etheta + 2*c4*Err + 2*c5*Ezz) + Pi - T1;
    stress_zz(index1 + 1) = C(3,3) * 0.5 * c * exp(Q) * (2*c3*Ezz + 2*c5*Etheta + 2*c6*Err) + Pi - T1;

  
    % Save R
    r_position(index1 +1) = r1;

    

end
end
