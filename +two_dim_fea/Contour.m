classdef Contour
    properties
        numEdges;
        Edges;
        Area;
    end
    methods
        function obj = Contour(nedges, edges, area)
            obj.numEdges = nedges;
            obj.Edges = edges;
            obj.Area = area;
        end
    end
end
