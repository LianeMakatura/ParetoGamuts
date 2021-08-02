function [str] = sprintArray(array, datatype)
%SPRINTARRAY print array of any length to a string
%   array - array to print
%   datatype - sprintf string code for the type of the elements in the array (e.g. 'd' for integers)
    element_format = sprintf("%%%s, ", datatype);

    contents = sprintf(element_format, array); % loops through every element in array, stores in same string
    contents = strip(strtrim(contents), 'right', ","); %trim off the trailing whitespace, then comma

    str = "[" + contents + "]";
end

