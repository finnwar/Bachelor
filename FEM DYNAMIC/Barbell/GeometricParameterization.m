%%Function to find parameters of the function describing the part geometry
% The curve is described by sinusoidal interpolation
function factors = GeometricParameterization(length_end, thickness_end, thickness_middle)
    A = [1 0 1 0 1; 1 sin(length_end) cos(length_end) sin(2*length_end) cos(2*length_end); 
        1 sin(2*length_end) cos(2*length_end) sin(4*length_end) cos(4*length_end);
        0 2 0 4 0;
        0 2*cos(2*length_end) -2*sin(2*length_end) 4*cos(4*length_end) -4*sin(4*length_end)];

    b = [thickness_end/2 (thickness_middle+thickness_end)/4 thickness_middle/2 0 0].';
    
    factors = A\b;

end