function [ points, triangles, rectangles ] = rectangular_object( m, n, width, height )
% nodal points and triangles of a rectangular body
    points = [];
    for j=1:n
        y = height*(j-1)/(n-1);
        for i=1:m
            x = width*(i-1)/(m-1);
            points = [points, [x;y]];
        end
    end

    triangles = [];
    for j=1:n-1
        for i=1:m-1
            triangles = [triangles; (j-1)*m+i, (j-1)*m+i+1, j*m+i];
        end
        for i=2:m
            triangles = [triangles; j*m+i, j*m+i-1, (j-1)*m+i];
        end
    end
    
    rectangles = [];
    for j=1:n-1
        for i=1:m-1
            rectangles = [rectangles; (j-1)*m+i, (j-1)*m+i+1, j*m+i+1, j*m+i];
        end
    end
end
