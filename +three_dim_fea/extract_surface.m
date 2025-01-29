function [ surface_triangles, volume ] = extract_surface (p, surfs, order, triangles, points)
%
% triangles and signed volume of the p-th surface
%
    index = find(surfs == p);
    triangles = triangles(index,:);
    order = order(index,:);
    [ ~, index ] = sort(order);
    surface_triangles = triangles(index,:);
    
    n = size(surface_triangles,1);
    volume = 0;
    for r=1:n
        xi = points(:, surface_triangles(r,1));
        xj = points(:, surface_triangles(r,2));
        xk = points(:, surface_triangles(r,3));
        volume = volume + det([xi,xj,xk]);
    end
    volume = -volume/6;
    % outer surface : positive value of signed volume
    % inner surface : negative value of signed volume
end
