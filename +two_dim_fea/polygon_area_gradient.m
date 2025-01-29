function [ area_x, area_y ] = polygon_area_gradient (coordinates)
    n = size(coordinates,2);
    x = coordinates(1,:)'; y = coordinates(2,:)';
    area_x = (1/2) * ( - [ y(n-1); y(n); y(1:n-2) ] + [ y(3:n); y(1); y(2) ] );
    area_y = (1/2) * (   [ x(n-1); x(n); x(1:n-2) ] - [ x(3:n); x(1); x(2) ] );
end
