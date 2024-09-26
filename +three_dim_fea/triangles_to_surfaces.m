%
% classfy edges via connectivity
% when triangle i and j are connected, surfs(i) == surfs(j)
%
function [ nsurfs, surfs, order ] = triangles_to_surfaces (triangles)
    ntris = size(triangles, 1);
    surfs = zeros(ntris,1);
    nsurfs = 0;
    order = zeros(ntris,1);
    
    while (1)
        p = find(surfs == 0);
        if isempty(p)
            break;
        end
        k = p(1);
        nsurfs = nsurfs + 1;
        surfs(k) = nsurfs;
        d = 0;
        while (1)
            candidates = intersect( find(surfs == nsurfs), find(order == 0) );
            if isempty(candidates)
                break;
            end
            %fprintf('%d : ', nsurfs); candidates'
            for k = candidates'
                for l = k+1:ntris
                    common = intersect(triangles(k,:), triangles(l,:));
                    if ~isempty(common)
                        surfs(l) = nsurfs;
                    end
                end
            end
            for k = candidates'
                d = d + 1;
                order(k) = d;
            end
            %order(candidates)'
        end
    end
end
