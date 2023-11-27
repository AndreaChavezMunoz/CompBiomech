%=========================================================================
% Lecture 16: In Class Tutorial
%
% This function calculates the radial equilibrium function for an axially 
% stretched and pressurized thick wall vessel and is part of the set of 
% equations you will implement for your vasculature project
%
% Input data: 
%   luminal pressure (Pi), axial stretch (lambdaz_v) 
%   material parameters, radii in ktf (Ri, Ro)
%
% Output data:
%   approximation of the outer radius, ro
%
% The inverse solution of the radial equilibrium involves finding
% the root of the equation:
%   Pi - int_{ri}^{ro} (tqq-trr)/r dr = 0
%=========================================================================
cd 'C:\Users\slee4\OneDrive\Documents\Fall 2023\BIOE 5640\Lecture 16';

close all, clear all, clc
%% Load and Define INPUT DATA

% Control refers to a age-matched healthy model
control_data = load("Input_HypertensionControl_ATA.mat");
test_data = load("Input_Hypertension_ATA.mat");

% Obtains material parameters
material_control = control_data.estimated_parameters;
material_test = test_data.estimated_parameters;

% Loads control parameters
control_lambdaz_exp = control_data.data_kl.lambdaz_exp;
control_Psys_exp = control_data.data_kl.Psys_exp;
control_ref_outer_radius = control_data.data_ktf.or_exp;
control_ref_inner_radius = control_data.data_ktf.ir_exp;
control_thickness = control_data.data_ktf.h_exp;

% Load test parameters
test_lambdaz_exp = test_data.data_kl.lambdaz_exp;
test_Psys_exp = test_data.data_kl.Psys_exp;
test_ref_outer_radius = test_data.data_ktf.or_exp;
test_ref_inner_radius = test_data.data_ktf.ir_exp;
test_thickness = test_data.data_ktf.h_exp;

%% Solve the equilibrium equation
% Equilibrium equation for a pressurized and axially stretched vessel
% Unknown: outer radius, ro
% ro will be estimated as the roots of the equilibrium equation

H=cell(1,2);
H{1,1}=@equilibrium_r_or_loaded; %handle to equilibrium equation 
H{1,2}=@equilibrum_z_fz_loaded;

% Parameter to be estimated in loaded configuration is the outer radius (ro) 
% x0 is the initial guess for ro based on the input experimental data

       

        
% DELIVERABLE - Plot Luminal Pressure vs Outer Diameter
luminalP = 10:2:140;
control_axialStretch = [control_lambdaz_exp, control_lambdaz_exp*1.05, control_lambdaz_exp*0.95];
test_axialStretch = [test_lambdaz_exp, test_lambdaz_exp*1.05, test_lambdaz_exp*0.95];

control_or_est = zeros(1, length(luminalP));
test_or_est = zeros(1, length(luminalP));

control_P_or = figure();
for stretch_index = 1:length(control_axialStretch)
    for P_index = 1:length(luminalP)
        control_or_est(1, P_index) = Newton_Raphson(H{1,1}, control_ref_outer_radius, control_ref_inner_radius, ...
                               control_axialStretch(1, stretch_index), luminalP(1, P_index), material_control, control_ref_outer_radius); % call the Newton Raphson function (Newton_Raphson_tutorial.m)
        plot(luminalP(1,P_index), control_or_est(1,P_index));
        hold on
        
    end

end

test_P_or = figure();
for stretch_index = 1:length(test_axialStretch)
    for P_index = 1:length(luminalP)
        test_or_est(1, P_index) = Newton_Raphson(H{1,1}, test_ref_outer_radius, test_ref_inner_radius, ...
                               test_axialStretch(1, stretch_index), luminalP(1, P_index), material_test, test_ref_outer_radius); % call the Newton Raphson function (Newton_Raphson_tutorial.m)
        plot(luminalP(1,P_index), test_or_est(1,P_index));
        hold on
        
    end

end
        
