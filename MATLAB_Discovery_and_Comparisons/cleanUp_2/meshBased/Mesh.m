classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vertices %nx3 list of vertices [x,y,z] each row
        faces %nx3 list of vertIDs (all tri only)
        
        numVerts
        numFaces
    end
    
    methods
        function obj = Mesh(vertices,faces)
            %MESH Construct an instance of this class
            obj.vertices = vertices;
            obj.faces = faces;
            
            obj.numVerts = length(obj.vertices);
            obj.numFaces = length(obj.faces);
        end
        
        function vert = getVertex(obj,vertID)
            %get vertID'th vertex from list
            vert = obj.vertices(vertID, :);
        end
        
        function face = getFace(obj,faceID)
            %get vertID'th vertex from list
            face = obj.faces(faceID, :);
        end
        
        function faceVerts = getFaceVerts(obj, faceIDs)
            % given a face (n vertex IDs), give back nx3 matrix of vertices
            faceVerts = [];
            for i=1:3
                faceVerts = [faceVerts; obj.getVertex(faceIDs(i))];
            end
        end
            
        function avgVertPos = getAvgVertPosition(obj)
            % get APPROX center of mass by averaging all vert positions
            % only approx because hi-def areas may skew vertex average
            avgVertPos = [0,0,0];
            for i=1:obj.numVerts
                avgVertPos = avgVertPos + obj.getVertex(i);
            end
            avgVertPos = avgVertPos / obj.numVerts;
        end
        
        function vol = getVolume(obj)
            % get mesh volume by using tetrahedral signed volume
            % turn each face into a tet wrt a common ref point
            
            ref = obj.getAvgVertPosition(); % use pt close to mesh for numerical stability; actual value or inside/outside doesn't matter
            
            vol = 0;
            for i=1:obj.numFaces
                faceVerts = getFaceVerts(obj, obj.getFace(i));
                v1 = faceVerts(1,:) - ref;
                v2 = faceVerts(2,:) - ref;
                v3 = faceVerts(3,:) - ref;
                
                vol = vol + obj.tetSignedVolume(v1, v2, v3);
            end
        end
        
        function vol = tetSignedVolume(~, v1, v2, v3)
            % get signed volume of tet defined by 3 vectors
            % tet is 1/6 of the parallelopiped spanned by vecs
            vol = dot(v1, cross(v2, v3)) / 6.0;
        end
      
    end
end

