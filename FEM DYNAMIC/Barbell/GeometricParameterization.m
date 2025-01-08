%%Function to find parameters of the function describing the part geometry
% The curve is described by a degree 4 polynomial of the form
% f(x)=ax^4+bx^3+cx^2+dx+e
function factors = GeometricParameterization(length_end, thickness_end, thickness_middle)
    A = [0 0 0 0 1;
         0 0 0 1 0;
         (length_end/2)^4 (length_end/2)^3 (length_end/2)^2 (length_end/2) 1;
         length_end^4 length_end^3 length_end^2 length_end 1;
         4*(length_end)^3 3*(length_end)^2 2*(length_end) 1 0];

    b = [thickness_end/2 0 (thickness_end+thickness_end)/2 thickness_middle 0].';
    
    factors = A\b;

end