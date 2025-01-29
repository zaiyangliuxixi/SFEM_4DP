function [ points, triangles, chamber1, chamber2, chamber3, chamber4 ] = ring_object_four_rooms(m, n, router, rinner, p, q)
    [ points, triangles ] = ring_object( m, n, router, rinner );
    npoints = size(points, 2);
    ntriangles = size(triangles, 1);
    
    bottom_center = points(:,1+0*n);
    bottom_right =  points(:,1+p*n);
    bottom_left =   points(:,1+(m-p)*n);
    %bottom_left, bottom_center, bottom_right
    
    top_center = points(:,1+(m/2)*n);
    top_right =  points(:,1+(m/2-p)*n);
    top_left =   points(:,1+(m/2+p)*n);
    %top_left, top_center, top_right
    
    left_center = points(:,1+(m/4*3)*n);
    left_bottom = points(:,1+(m/4*3+p)*n);
    left_top =    points(:,1+(m/4*3-p)*n);
    %left_bottom, left_center, left_top
    
    right_center = points(:,1+(m/4)*n);
    right_bottom = points(:,1+(m/4-p)*n);
    right_top   =  points(:,1+(m/4+p)*n);
    %right_bottom, right_center, right_top
    
    [ center_pt, center_tri ] = rectangular_object(2*p+1,2*p+1, top_right(1)-bottom_left(1), right_top(2)-left_bottom(2));
    center_pt = center_pt - center_pt(:,2*p^2+2*p+1); % { (2p+1)^2 +1 }/2 = 2p^2 + 2*p + 1
    
    [ left_pt, left_tri ] = rectangular_object(q,2*p+1, center_pt(1,1)-left_bottom(1), right_top(2)-left_bottom(2));
    left_pt = left_pt + left_bottom;
    for k=-p+1:p-1
        l = -k*q + (1+p*q);
        left_pt(:,l) = points(:,1+(m/4*3+k)*n);
    end
    
    [ right_pt, right_tri ] = rectangular_object(q,2*p+1, right_bottom(1)-center_pt(1,2*p+1), right_top(2)-left_bottom(2));
    right_pt = right_pt + center_pt(:,2*p+1);
    for k=-p+1:p-1
        l = k*q + (p+1)*q;
        right_pt(:,l) = points(:,1+(m/4+k)*n);
    end
    
    [ top_pt, top_tri ] = rectangular_object(2*p+1,q, top_right(1)-bottom_left(1), top_right(2)-center_pt(2,(2*p+1)*(2*p)+1));
    top_pt = top_pt + center_pt(:,(2*p+1)*(2*p)+1);
    %top_pt(:,q*(2*p+1)-1) = top_center;
    for k=-p+1:p-1
        l = -k + (2*p+1)*q - p;
        top_pt(:,l) = points(:,1+(m/2+k)*n);
    end
    
    [ bottom_pt, bottom_tri ] = rectangular_object(2*p+1,q, top_right(1)-bottom_left(1), center_pt(2,2*p+1)-bottom_left(2));
    bottom_pt = bottom_pt + bottom_left;
    %bottom_pt(:,2*p) = bottom_center;
    for k=-p+1:p-1
        l = k+(p+1);
        bottom_pt(:,l) = points(:,1+(k+m*(k<0))*n);
    end
    
    pt = center_pt; tri = center_tri;
    [pt, tri] = combine_objects(pt, tri, [], bottom_pt, bottom_tri, []);
    [pt, tri] = combine_objects(pt, tri, [], left_pt, left_tri, []);
    [pt, tri] = combine_objects(pt, tri, [], right_pt, right_tri, []);
    [pt, tri] = combine_objects(pt, tri, [], top_pt, top_tri, []);

    [points, triangles] = combine_objects(points, triangles, [], pt, tri, []);

    square_bottom_left = m*n+1;
    square_bottom_right = m*n+(2*p+1);
    square_top_left = m*n+(2*p+1)*(2*p)+1;
    square_top_right = m*n+(2*p+1)*(2*p+1);
    
    square_last = m*n+(2*p+1)*(2*p+1);
    nrect = (2*p+1)*(q-2);
    
    bottom_start = square_last+0*nrect+1;
    bottom_left = bottom_start:(2*p+1):bottom_start+(2*p+1)*((q-2)-1);
    bottom_right = bottom_left + (2*p);
    
    left_start = square_last+1*nrect+1;
    left_bottom = left_start:1:left_start+(q-2)-1;
    left_top = left_bottom + (2*p)*(2*p+1);
    
    right_start = square_last+2*nrect+1;
    right_bottom = right_start:1:right_start+(q-2)-1;
    right_top = right_bottom + (2*p)*(2*p+1);
    
    top_start = square_last+3*nrect+1;
    top_left = top_start:(2*p+1):top_start+(2*p+1)*((q-2)-1);
    top_right = top_left + (2*p);

    %chamber1 = [ 137, 162:1:164, 37:4:61, 173:(-3):167 ];
    %chamber2 = [ 135, 165:3:171, 69:4:93, 153:1:155 ];
    %chamber3 = [ 129, 149:(-1):147, 101:4:125, 138:3:144 ];
    %chamber4 = [ 131, 146:(-3):140, 5:4:29, 158:(-1):156 ];
    chamber1 = [ square_top_right, right_top, (1+(m/4+p)*n):n:1+(m/2-p)*n, flip(top_right) ];
    chamber2 = [ square_top_left, top_left, 1+(m/2+p)*n:n:1+(m/4*3-p)*n, left_top ];
    chamber3 = [ square_bottom_left, flip(left_bottom), 1+(m/4*3+p)*n:n:1+(m/4*4-p)*n, bottom_left ];
    chamber4 = [ square_bottom_right, flip(bottom_right), 1+p*n:n:1+(m/4-p)*n, flip(right_bottom) ];
end
