function [ ncontours, contours, order ] = edges_to_contours (edges)
    nedges = size(edges, 1);
    contours = zeros(nedges,1);  %contours是列向量
    ncontours = 0;
    order = zeros(nedges,1);
        
    while (1)
        p = find(contours == 0);
        if isempty(p)
            break;
        end
        k = p(1);
        ncontours = ncontours + 1;
        contours(k) = ncontours;
        d = 1;
        order(k) = d;
        while (1)
            nexts = find( edges(:,1) == edges(k,2) );
            if isempty(nexts)
                break;
            end
            next = nexts(1);
            if order(next) ~= 0
                break;
            end
            contours(next) = ncontours;
            d = d + 1;
            order(next) = d;
            k = next;
        end
    end     
end
