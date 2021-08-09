classdef cppCaller < handle
    properties  
        % Static
        rD
        rd
        functionID
        
        %Dynamic
        targets
                
        inFileName;
        outFileName;
    end
    methods

    function obj = cppCaller(rD, rd, functionID)
        obj.rD = rD;
        obj.rd = rd;
        obj.functionID = functionID;
        obj.inFileName = ['samplesFile_' num2str(functionID) '.txt'];
        obj.outFileName = ['samplesFile_' num2str(functionID) '_results.txt'];
    end

    % Writes @points to a text file to be read by the executable
    function writeToFile(obj, points)
        if obj.functionID == 3
            points(:,5) = 0.01;
        end
        dlmwrite(obj.inFileName, points, ' '); 
    end

    % Helper function that returns each set of points in the sample
    % results file as a vector
    function contents = readFromFile1(obj)
        contents = dlmread(obj.outFileName);
    end

    % Helper function that returns each set of points in the sample
    % results file as a vector
    function [evaluated, gradient] = readFromFile2(obj)
        contents = dlmread(obj.outFileName);
        evaluated = contents(:,1:obj.rd);
        gradient = contents(:,obj.rd+1:obj.rd + obj.rd*obj.rD);
    end

    % Helper function that returns each set of points in the sample
    % results file as a vector
    function [evaluated, gradient, hess] = readFromFile3(obj)
        contents = dlmread(obj.outFileName);
        evaluated = contents(:,1:obj.rd);
        gradient = contents(:,obj.rd+1:obj.rd + obj.rd*obj.rD);
        hess  = contents(:, obj.rd + obj.rd*obj.rD + 1:obj.rd + obj.rd*obj.rD + obj.rd*obj.rD^2);
    end

    % Calls the evaluation executable
    function callExec1(obj)
        filename = ['cppCalls\perfomanceSampler.exe samplesFile_' num2str(obj.functionID) ' 1 ' num2str(obj.functionID)];
        system(filename);
    end

    % Calls the gradient executable
    function callExec2(obj)
      filename = ['cppCalls\perfomanceSampler.exe samplesFile_' num2str(obj.functionID) ' 2 ' num2str(obj.functionID)];
        system(filename);
    end

    % Calls the hessian executable 
    function callExec3(obj)
      filename = ['cppCalls\perfomanceSampler.exe samplesFile_' num2str(obj.functionID) ' 3 ' num2str(obj.functionID)];
        system(filename);
    end

    % Returns the single-objective scalarization function as a function
    % handle to be passed into fmincon
    function fh = getOptimizationHandle(obj, targets)
        obj.targets = targets;
        fh = @obj.magicFunction;
    end

    % Function representing the single-objective scalarization function
    % for local optimization. Explicitly equal to:
    %  (f_0-target_0)^2 + (f_1-target_1)^2 + ...(f_rd - target_rd)^2
    function [f, gradf] = magicFunction(obj, x)
        %[evalPts, gradPts, hessPts] = obj.evalPointsHess(x);
        [evalPts, gradPts] = obj.evalPointsGradient(x);
        f = sum((evalPts-obj.targets).^2);
        gradf = 2*sum((evalPts-obj.targets).*gradPts, 2);

    end

    % @param
    % points: N x rD matrix representing N sets of inputs
    % :return 
    % evalPoints: N x rd matrix representing the values of the
    % objective functions evaluated at each point        
    function output = evalPoints(obj, points)
        obj.writeToFile(points);
        obj.callExec1();

        output = obj.readFromFile1();
    end

    % @param
    % points: N x rD matrix representing N sets of inputs
    % :return 
    % evalPoints: N x rd matrix representing the values of the
    % objective functions evaluated at each point
    % gradPoints: rD x rd matrix representing the gradients of the
    % objective functions
    function [evalPoints, gradPoints] = evalPointsGradient(obj, points)
        obj.writeToFile(points);
        obj.callExec2();
        [evalPoints, gradPointsFlat] = obj.readFromFile2();
        gradPoints = reshape(gradPointsFlat, [obj.rD, obj.rd]);
    end

    % @param
    % points: N x rD matrix representing N sets of inputs
    % :return 
    % evalPoints: N x rd matrix representing the values of the
    % objective functions evaluated at each point
    % gradPoints: rD x rd matrix representing the gradients of the
    % objective functions
    % hessPoints: rD x rD x rd multi-dimensional array representing the
    % hessians of the objective functions
    function [evalPoints, gradPoints, hessPoints] = evalPointsHessian(obj, points)
        obj.writeToFile(points);
        obj.callExec3();
        [evalPoints, gradPointsFlat, hessPointsFlat] = obj.readFromFile3();
        gradPoints = reshape(gradPointsFlat, [obj.rD, obj.rd]);

        hessPoints = zeros(obj.rD,obj.rD,obj.rd);
        for i = 0:obj.rd-1
            idx = (i*obj.rD^2) + 1;
            hess_i = reshape(hessPointsFlat(idx:idx+obj.rD^2-1),[obj.rD, obj.rD]);
            hessPoints(:,:,i+1) = hess_i;
        end
    end
    end
end