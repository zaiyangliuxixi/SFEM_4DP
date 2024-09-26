function [ points, triangles ] = ring_object( m, n, router, rinner  )
% nodal points and triangles of a ring body
% the origin is the center of the ring
    points = [];
    dradius = (router - rinner) / (n-1);
    dtheta = 2*pi/m;
    for i=1:m
        theta = (i-1)*dtheta - pi/2;
        for j=1:n
            radius = rinner + (j-1)*dradius;
            points = [ points, radius*[ cos(theta); sin(theta) ] ];
        end
    end
    triangles = [];
    for i=1:m-1
        k = n*(i-1) + 1; kd = k+n;
        for j=1:n-1
            triangles = [ triangles; k, k+1, kd+1 ];
            triangles = [ triangles; k, kd+1,  kd ];
            k = k + 1; kd = kd +1;
        end
    end
    k = n*(m-1) + 1; kd = 1;
    for j=1:n-1
        triangles = [ triangles; k, k+1. kd+1 ];
        triangles = [ triangles; k, kd+1,  kd ];
        k = k + 1; kd = kd + 1;
     end
end
