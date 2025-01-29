function area = polygon_area (coordinates)
    n = size(coordinates,2);
    x = coordinates(1,:)'; y = coordinates(2,:)';
    dx = x - [ x(2:end); x(1) ];    % x_{k} - x_{k+1}
    dy = y + [ y(2:end); y(1) ];    % y_{k} + y_{k+1}
    area = (1/2)*sum(dx .* dy);
end
