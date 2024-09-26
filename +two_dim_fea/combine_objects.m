function [ points, triangles, rectangles ] = combine_objects(points_a, triangles_a, rectangles_a, points_b, triangles_b, rectangles_b, eps)
    arguments
        points_a; triangles_a; rectangles_a; points_b; triangles_b; rectangles_b;
        eps = 1e-8;
    end
    npoints_a = size(points_a, 2);
    npoints_b = size(points_b, 2);
    index = zeros(1,npoints_b);
    for j = 1:npoints_b
        for i = 1:npoints_a
            d = points_a(:,i) - points_b(:,j);
            dist = sqrt(d.*d);
            if dist <= eps
                index(j) = i;
                break;
            end
        end
    end
    p = npoints_a + 1;
    for j = 1:npoints_b
        if index(j) == 0
            index(j) = p;
            p = p + 1;
        end
    end

    triangles = [ triangles_a; index(triangles_b) ];
    rectangles = [ rectangles_a; index(rectangles_b) ];
  
    points = points_a;
    for j = 1:npoints_b
        if (index(j) > npoints_a)
            points = [ points, points_b(:,j) ];
        end
    end
    
end
