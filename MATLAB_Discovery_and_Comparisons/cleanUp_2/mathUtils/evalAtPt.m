function [evalMat] = evalAtPt(symbMat, xVals, zVals)
    p = num2cell([xVals, zVals]);
     evalMat = eval(symbMat(p{:}));
end