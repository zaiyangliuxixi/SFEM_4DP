classdef SubRegion
    properties
        numTriangles;   Index_Triangles;
        numRectangles;  Index_Rectangles;
        numNodalPoints; Index_NodalPoints;
        Density; lambda; mu; lambda_vis; mu_vis; strain_mea;
        Partial_J_lambda; Partial_J_mu;
        color;
    end
    methods
        function obj = SubRegion (index, index_rects, index_npoints)
            arguments
                index;
                index_rects = [];
                index_npoints = [];
            end
            obj.Index_Triangles = index;
            obj.Index_Rectangles = index_rects;
            obj.Index_NodalPoints = index_npoints;
            obj.numTriangles = length(index);
            obj.numRectangles = length(index_rects);
            obj.numNodalPoints = length(index_npoints);
        end
    end
end
