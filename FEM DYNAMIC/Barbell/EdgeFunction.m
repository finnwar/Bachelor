%%Function describing the upper and lower edge of the part

function y = EdgeFunction(factors, x, length_end, length_middle, thickness_middle,thickness_end)
 total_length = 2*length_end + length_middle;
    if x<0
        error("Out of Bounds")
    if x>(2*length_end+length_middle)
        error("Out of Bounds")
    if (x<(2*length_end+length_middle)) && (x>=(length_end+length_middle))
        y = [factors.'*[(-x+total_length)^4 (-x+total_length)^3 (-x+total_length)^2 (-x+total_length) 1];
             -factors.'*[(-x+total_length)^4 (-x+total_length)^3 (-x+total_length)^2 (-x+total_length) 1]];
    if (x<(length_end+length_middle))&&(x>=length_end)
        y = [thickness_middle/2; -thickness_middle/2];
    if x<length_end
        y = [factors.'*[x^4 x^3 x^2 x 1]; -factors.'*[x^4 x^3 x^2 x 1]];
    if
end