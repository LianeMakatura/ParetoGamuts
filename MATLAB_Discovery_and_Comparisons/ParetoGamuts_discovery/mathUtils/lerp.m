function [val] = lerp(t, min, max)
%LERP linear interpolation, to rescale a parameter into another range
    val = min + t.*(max-min);
end

