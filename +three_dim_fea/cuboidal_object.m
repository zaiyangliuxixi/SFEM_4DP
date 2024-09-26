function [ points, tetrahedra, cuboids ] = cuboidal_object( l, m, n, a, b, c )
% nodal points and tetrahedra of a cuboidal body
% 直方体の節点と四面体のデータ
    points = [];
    for k=1:n
        z = c*(k-1)/(n-1);
        for j=1:m
            y = b*(j-1)/(m-1);
            for i=1:l
                x = a*(i-1)/(l-1);
                points = [points, [x;y;z]];
            end
        end
    end
    
    tetrahedra = [];
    for k=1:n-1
        for j=1:m-1
            for i=1:l-1
                down = 1 + l*m*(k-1) + l*(j-1) + (i-1);
                up = down + l*m;
                tetras = dividing_cuboid( down, down+1, down+l+1, down+l, ...
                                            up,   up+1,   up+l+1,   up+l);
                tetrahedra = [tetrahedra; tetras];
            end
        end
    end
    
    cuboids = [];
    for k=1:n-1
        for j=1:m-1
            for i=1:l-1
                down = 1 + l*m*(k-1) + l*(j-1) + (i-1);
                up = down + l*m;
                cuboid = [ down, down+1, down+l+1, down+l, ...
                             up,   up+1,   up+l+1,   up+l ];
                cuboids = [cuboids; cuboid];
            end
        end   
    end
end

function tetras = dividing_cuboid (i, j, k, l, m, n, r, s)
    tetras = [
        m, j, i, l;
        m, j, l, s;
        m, j, s, n;
        s, k, r, n;
        s, k, n, j;
        s, k, j, l;
    ];
end
