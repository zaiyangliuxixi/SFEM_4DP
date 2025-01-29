classdef Surface
    properties
        numTriangles;
        Triangles;
        Volume;
    end
    methods
        function obj = Surface(ntris, tris, volume)
            obj.numTriangles = ntris;
            obj.Triangles = tris;
            obj.Volume = volume;
        end
    end
end
