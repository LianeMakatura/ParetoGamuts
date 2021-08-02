function [perf_json, des_json] = write_json(buffers, rD, rA, rd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     num_z_steps = 10;
%     min_z = 0; 
%     max_z = 1;
%     z_step = (max_z - min_z) / (num_z_steps);
%     
%     app_var_values = [];
%     for i=0:num_z_steps
%         app_var_values = min_z + z_step * i;
%     end
    
    % ==== put all (pareto optimal) buffer points in single array
    numBuffs = length(buffers);
    allParetoFronts = [];
    allParetoSets = [];


    z_front = [];
    z_set = [];
    for buffIdx=1:numBuffs
        b = buffers(buffIdx);
        paretoIndices = b.getParetoInd;

        for j = 1:length(paretoIndices)
            entry = b.Buff(paretoIndices(j));
            patch = entry.bestPatch;
            if  patch > 0 
                % Store all the points in a single matrix
                newDes = entry.minD(1:rD);
                newPerf = entry.minF;
                app_val = entry.minD(rD+1:rD+rA);
                z_front = [z_front; newPerf, app_val];
                z_set = [z_set; newDes, app_val];
            end
        end
        app_var_value = entry.minD(rD+1:rD+rA);

    end
%         allParetoFronts = [allParetoFronts; z_front];
%         allParetoSets = [allParetoSets; z_set];
%     if exist('performance.mat', 'file') == 2 %already exists
%         disp('existed!')
%         save('performance.mat', 'z_front', '-append')
%         save('design.mat', 'z_set', '-append')
%         save('app_vals', 'app_var_value', '-append')
%     else
        save('performance.mat', 'z_front')
        save('design.mat', 'z_set')
        save('app_vals', 'app_var_value')
%     end
    
%     perf_json = jsonencode(table(allAppVars,allParetoFronts));
%     des_json = jsonencode(table(allAppVars,allParetoSets));
    
%     fileID = fopen('sumSines_perf.json', 'w');
%     fprintf(fileID, perf_json);
%     fclose(fileID);
%     
%     fileID = fopen('sumSines_des.json', 'w');
%     fprintf(fileID, des_json);
%     fclose(fileID);
end

