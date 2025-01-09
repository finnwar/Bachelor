%%Function describing the upper and lower edge of the part

function y = EdgeFunction(factors, x, length_end, length_middle, thickness_middle)
 total_length = 2*length_end + length_middle;
    if x<0
        error("Out of Bounds")
    elseif x>(2*length_end+length_middle)
        error("Out of Bounds")
    elseif (x<=(2*length_end+length_middle)) && (x>=(length_end+length_middle))
        y = [factors.'*[1 sin(2*((length_middle+2*length_end)-x)) cos(2*((length_middle+2*length_end)-x)) sin(4*((length_middle+2*length_end)-x)) cos(4*((length_middle+2*length_end)-x))].'
             -factors.'*[1 sin(2*((length_middle+2*length_end)-x)) cos(2*((length_middle+2*length_end)-x)) sin(4*((length_middle+2*length_end)-x)) cos(4*((length_middle+2*length_end)-x))] .'];
    elseif (x<(length_end+length_middle))&&(x>=length_end)
        y = [thickness_middle/2; -thickness_middle/2];
    elseif x<length_end
        y = [factors.'*[1 sin(2*x) cos(2*x) sin(4*x) cos(4*x)].'; -factors.'*[1 sin(2*x) cos(2*x) sin(4*x) cos(4*x)].'];
    end
end