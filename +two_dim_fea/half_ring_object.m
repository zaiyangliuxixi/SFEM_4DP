function [ points, triangles ] = half_ring_object( m, n, router, rinner, thetastart )
% nodal points and triangles of a half-ring body
    points = [];
    dradius = (router - rinner) / (n-1);
    dtheta = pi/(m-1);
    for i=1:m
        theta = (i-1)*dtheta + thetastart;
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
end
