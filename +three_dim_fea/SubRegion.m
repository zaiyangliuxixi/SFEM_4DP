classdef SubRegion
    properties
        numTetrahedrons;   Index_Tetrahedrons;
        numNodalPoints; Index_NodalPoints;
        Density; lambda; mu; lambda_vis; mu_vis;swelling_ratio;swelling_time;strain_mea;
        Partial_J_lambda; Partial_J_mu;
        color;
    end
    methods
        function obj = SubRegion (index, index_npoints)
            arguments
                index;
                index_npoints = [];
            end
            obj.Index_Tetrahedrons = index;
            obj.Index_NodalPoints = index_npoints;
            obj.numTetrahedrons = length(index);
            obj.numNodalPoints = length(index_npoints);
        end
    end
end
