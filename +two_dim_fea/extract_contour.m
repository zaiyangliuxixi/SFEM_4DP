function [ contour_edges, area ] = extract_contour (p, countours, order, edges, points)
%
% edges and signed area of the p-th contour
%
    index = find(countours == p);
    edges = edges(index,:);
    order = order(index,:);
    [ ~, index ] = sort(order);
    contour_edges = edges(index,:);
    
    n = size(contour_edges,1);
    area = 0;
    for k=1:n
        xs = points(:, contour_edges(k,1));
        xe = points(:, contour_edges(k,2));
        area = area + det([xs,xe]);
    end
    area = area/2;
    % outer contour : positive value of signed area
    % inner contour : negative value of signed area
end
