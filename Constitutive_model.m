
function [sigma_extra] = Constitutive_model(F, parameters)

%=========================================================================
% This function calculates the Cauchy stress extra (sigma_extra) for a  
% given gradient tensor F, Using Neo-Hookean Constitutive model (see Lecture 8)
%=========================================================================
sigma_extra = zeros(3,3);

c1 = parameters(1);
c2 = parameters(2);
c3 = parameters(3);
c4 = parameters(4);
c5 = parameters(5);
c6 = parameters(6);
c = parameters(7);

C = F'*F;
E = 0.5 * (C - eye(3,3));
Err = E(1,1);
Etheta =  E(2,2);
Ezz = E(3,3);

Q = (c1 * Err^2) + (c2 * Etheta^2) + (c3 * Ezz^2) + (2*c4 * Err*Etheta) + (2*c5 * Etheta*Ezz) + (2*c6 * Err*Ezz);

%Cauchy extra stress
sigma_extra(1,1) = (0.5 * c * exp(Q)) * ((2*Err*c1) + (2*c4*Etheta) + (2*c6*Ezz)); %r component
sigma_extra(2,2) = (0.5 * c * exp(Q)) * ((2*c2*Etheta) + (2*c4*Err) + (2*c5*Ezz)); %theta component
sigma_extra(3,3) = (0.5 * c * exp(Q)) * ((2*c3*Ezz) + (2*c5*Etheta) + (2*c6*Err)); %z component

end


