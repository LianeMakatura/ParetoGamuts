classdef showFigures < handle
properties
    rD
    rd
    buffer
    functionID
end
methods

    function obj = showFigures(rD, rd, buffer, functionID)
        obj.rD = rD;
        obj.rd = rd;
        obj.buffer = buffer;
        obj.functionID = functionID;
    end
    
    %% Display the samples for this iteration before and after performing local optimization
    %
    %  @param samples: n x rD+rd list of n points
    %  @param projections: same points as @samples after being locally
    %  optimized
    %
    %  samples displayed in red, projections displayed in blue
    function showCurrentIteration(obj, samples, projections)
        figure;
        hold on;
        diff = projections-samples;
        if obj.rd == 2
            scatter(samples(:,obj.rD+1),samples(:,obj.rD+2), 'red');
            scatter(projections(:,obj.rD+1),projections(:,obj.rD+2), 'blue');
            quiver(samples(:,obj.rD+1),samples(:,obj.rD+2),diff(:,obj.rD+1),diff(:,obj.rD+2),0, 'Color', 'g', 'MaxHeadSize', .05);
        else
            scatter3(samples(:,obj.rD+1),samples(:,obj.rD+2), samples(:,obj.rD+3), 'red');
            scatter3(projections(:,obj.rD+1),projections(:,obj.rD+2),projections(:,obj.rD+3), 'blue');
            quiver3(samples(:,obj.rD+1),samples(:,obj.rD+2), samples(:,obj.rD+3),diff(:,obj.rD+1),diff(:,obj.rD+2),diff(:,obj.rD+3),0, 'Color', 'g', 'MaxHeadSize', .05);
        end
        obj.plotGroundTruth();
    end
    
    %% Display the 1/2D image of the patch labels embedded in the buffer
    function showLabeledBuffer(obj)
        B = obj.buffer;
        labels = B.emptyCells*0;
        indices = find(B.emptyCells == 0)';
        for idx = indices
            labels(idx) = B.selectedPatches(idx);
        end
        figure; 
        imagesc(labels);
    end
    
    %% Display the labeled Pareto set (design space) after performing graph-cut
    %  Performs MDS using PCA for D > 3
    function showLabeledParetoSet(obj, keepAllPoints)
        rng(2);
        B = obj.buffer;
        if nargin < 2
          keepAllPoints = true;
        end

        if keepAllPoints
            idx = B.validBufferIndices;
        else
            idx = B.getParetoInd;
        end
              
        numPatches = B.nPatches;
        a = cell(numPatches, 1);

        % Perform multidimensional scaling using PCA
        if obj.rD > 3
            D = zeros(size(idx,1),obj.rD);
            for q = 1:size(idx,1)
                minD = B.Buff(idx(q)).minD;
                D(q,:) = minD;
            end
            [~, reconstructed] = pca(D);
            reduced = reconstructed(:,1:3); 
        end
          
        % Group all points into cells for faster rendering
        for f = 1:length(idx)
            C = B.Buff(idx(f));
            patchIndex = C.bestPatch;
            if obj.rD > 3
                minD = reduced(f,:);
            else
                minD = C.minD;
            end
            if patchIndex > 0
                a{patchIndex} = [a{patchIndex};minD];
            end
        end
          
        % Display points
        colors = rand(B.nPatches, 3);
        figure;
        hold on;
        for t = 1:numPatches
            des = a{t};
            if size(des,1) > 0
                scatter3(des(:,1),des(:,2),des(:,3), 'MarkerEdgeColor',colors(t,:));
            end
        end
    end
    
    %% Display the labeled Pareto front (objective space) after performing graph-cut
    function showLabeledParetoFront(obj, keepAllPoints, new_fig_window, fig_id)
        rng(2);
        B=obj.buffer;
        
        disp(nargin)
        if (nargin <2)
            keepAllPoints = true;
        end
        
        if (nargin <3)
            new_fig_window = true;
        end
        
        idx = B.validBufferIndices;
        if (keepAllPoints)
            idx = B.getParetoInd();
        end
        
        colors = rand(B.nPatches, 3);
        if (new_fig_window)
            figure;
            hold on;        
        elseif exist('fig_handle', 'var')
            figure(fig_id);
        end
        
        % Group all points in individual cells for faster rendering
        a = cell(B.nPatches);
        ids  = zeros(B.nPatches,1);
        for i=idx'
            patchID = B.Buff(i).bestPatch;
            F = B.Buff(i).minF;
            if (patchID > 0)
                ids(patchID) =1;
                a{patchID} =  [a{patchID}; F];
            end
        end        

        % Display result
        for i = find(ids)'
            F = a{i};
            if(obj.rd == 2)
                scatter(F(:,1),F(:,2),'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));   
            else
                scatter3(F(:,1),F(:,2), F(:,3), 'MarkerEdgeColor', colors(i,:));                
            end
        end         
    end
    
    
    %% Display the labeled Pareto front (objective space) after performing graph-cut
    function showLabeledParetoFrontStacked(obj, z, keepAllPoints, new_fig_window)
        rng(2);
        B=obj.buffer;
        
        if (nargin <3)
            keepAllPoints = true;
        end
        
        if (nargin <4)
            new_fig_window = true;
        end
        
        idx = B.validBufferIndices;
        if (keepAllPoints)
            idx = B.getParetoInd();
        end
        
        colors = rand(B.nPatches, 3);
        if (new_fig_window)
            figure;
            %hold on;
        end
        
        % Group all points in individual cells for faster rendering
        a = cell(B.nPatches);
        ids  = zeros(B.nPatches,1);
        for i=idx'
            patchID = B.Buff(i).bestPatch;
            F = B.Buff(i).minF;
            if (patchID > 0)
                ids(patchID) =1;
                a{patchID} =  [a{patchID}; F];
            end
        end        

        % Display result
        for i = find(ids)'
            F = a{i};
            if(obj.rd == 2)
                Z = repmat(z, size(F, 1));
                scatter3(F(:,1),F(:,2), Z(:,1), 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));   
            else
                disp('Cannot visualize stacked pareto with rd > 2')              
            end
        end         
    end
    
    
    %% Display the current Pareto set (design space) as stored in the buffer
    %  Performs MDS using PCA for D > 3
    function showDiscoveredParetoSet(obj, keepAllPoints)
        b = obj.buffer;
        idx = b.validBufferIndices;
        if (keepAllPoints)
            idx = b.getParetoInd();
        end
        
        % Perform multidimensional scaling using PCA
        D = zeros(size(idx,1),obj.rD);
        for q = 1:size(idx,1)
            minD = b.Buff(idx(q)).minD;
            D(q,:) = minD;
        end
        
        if obj.rD > 3
            [~, reconstructed] = pca(D);
            paretoSet = reconstructed(:,1:3);  
        else
            paretoSet = D;
        end
        
        % Display the buffer points
        figure;
        hold on; grid on;
        scatter3(paretoSet(:,1),paretoSet(:,2), paretoSet(:,3),'r', 'filled');                
        axis equal; 
    end
     
    %% Display the current Pareto front (objective space) as stored in the buffer
    function showDiscoveredParetoFront(obj, keepAllPoints, newFigure)
        b = obj.buffer;
        ind = b.validBufferIndices;
        if (keepAllPoints)
            ind = b.getParetoInd();
        end
        
        % Put the patch points into a vector first for faster rendering
        paretoFront = [];
        for i = ind'
            if( b.emptyCells(i) == 0)
                bStruct = b.Buff(i);
                paretoFront = [paretoFront; bStruct.minF];
            end
        end
        
        % Display the buffer points
        if newFigure
            figure;
            hold on;            
        end
        
        if(obj.rd == 2)
            scatter(paretoFront(:,1),paretoFront(:,2),'r', 'filled');                
        else
            scatter3(paretoFront(:,1),paretoFront(:,2), paretoFront(:,3),'r', 'filled');                
        end
        obj.plotGroundTruth();
        axis equal; 
       
    end
    
    function plotGroundTruth(obj)
        id = obj.functionID;
        if id < 6
            obj.plotZDTGroundTruth(id);
        elseif id < 9
            obj.plotDTLZGroundTruth(id);
        else
            return
        end
    end
    
    
    %% Helper functions for plotting benchmark ground truths
    function plotZDTGroundTruth(obj, id)
        x = linspace(0,1,100);
        if id == 1
            ground_truth = (1 - sqrt(x))/10; %ground truth ZDT 1
        elseif id == 2
            ground_truth = (1 - x.^2)/10; %ground truth ZDT 2
        elseif id == 3
            % Ground truth for the ZDT 3 problem is disconnected on these
            % intervals
            x = linspace(0,.083001,10);
            x2 = linspace(.1822, .2577,10);
            x3 = linspace(.4093, .4538,10);
            x4 = linspace(.6183, .6525,10);
            x5 = linspace(.8233, .85183,10);

            ground_truth =  (1 - sqrt(x) - x.*(sin(10*pi*x)))/1; %ground truth ZDT 3
            ground_truth2 =  (1 - sqrt(x2) - x2.*(sin(10*pi*x2)))/1; %ground truth ZDT 3
            ground_truth3 =  (1 - sqrt(x3) - x3.*(sin(10*pi*x3)))/1; %ground truth ZDT 3
            ground_truth4 =  (1 - sqrt(x4) - x4.*(sin(10*pi*x4)))/1; %ground truth ZDT 3
            ground_truth5 =  (1 - sqrt(x5) - x5.*(sin(10*pi*x5)))/1; %ground truth ZDT 3

            plot(x,ground_truth, 'black','LineWidth',1);
            plot(x2,ground_truth2, 'k','LineWidth',1);
            plot(x3,ground_truth3, 'k','LineWidth',1);
            plot(x4,ground_truth4, 'k','LineWidth',1);
            plot(x5,ground_truth5, 'k','LineWidth',1);    
            return;
        elseif id == 4
            ground_truth = (1 - sqrt(x))/(obj.rD*30); %ground truth ZDT 4
        elseif id == 5
            x = linspace(.28, 1, 100);
            ground_truth = (1 - x.^2)/10; %ground truth ZDT 
        else
            return
        end
        plot(x,ground_truth, 'black','LineWidth',1);
    end

    function plotDTLZGroundTruth(obj, id)
        hold on;
        view(10, 60);
        if id == 6
            [x,y,z] = sphere;
            x(x < 0) = 0;
            y(y < 0) = 0;
            z(z < 0) = 0;

            surf(x/2.25, y/2.25, z/2.25);   
        elseif id == 7
            x = linspace(0,1,30);
            prs = permn(x,3);
            simplex =prs(abs(sum(prs,2)-1) < .01,:);
            simplex = simplex/(240*(obj.rD-2));
            a = simplex(:,1);
            b = simplex(:,2);
            c = simplex(:,3);
            tri = delaunay(a, b);
            trimesh(tri, a, b, c); 
            axis([0 0.05 0 0.05 0 0.05]);
        elseif id == 8
            [x,y,z] = sphere;
            x(x < 0) = 0;
            y(y < 0) = 0;
            z(z < 0) = 0;
            scalingFactor = 220*(obj.rD-2);
            surf(x/scalingFactor, y/scalingFactor, z/scalingFactor);  
        end
    end
    end
end
