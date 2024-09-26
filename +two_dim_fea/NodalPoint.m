classdef NodalPoint
    properties
        Coordinates;
        Displacement;
        Velocity
    end
    methods
        function obj = NodalPoint(p)
            obj.Coordinates = p;
        end
    end
end
