%%=========================================================================
% Main Script of Bachelor Thesis of Finn Warlimont
%%=========================================================================

% Dimensions

length_total = 0.400;     %Barbell length [m]
length_middle = 0.200;    %Middle part length [m]
end_length = (length_total-length_middle)/2;


thickness_end = 0.2;    %Thickness of both ends [m]
thickness_neck = 0.1;   %Thickness of thin middle part [m]

interpolation_factors = GeometricParameterization(end_length,thickness_end, thickness_neck);

% Material Constants

rho = 7850;      %Density [kg/m^3]
E = 210e9;       %Young's Modulus [N/m^2]
nu = 0.3;        %Contraction [-]

% Gravity

g = 9.81;       %Gravitational constant [m/s^2]







