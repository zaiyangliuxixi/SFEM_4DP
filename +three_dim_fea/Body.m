classdef Body
    properties
        numNodalPoints; NodalPoints;
        numTetrahedrons; Tetrahedrons;
        BoundaryFaces;
        numSubRegions; SubRegions;
        numSurfaces; Surfaces;
        Density; lambda; mu; lambda_vis; mu_vis;
        strain_potential_energy;
        gravitational_potential_energy;
        J_lambda; J_mu;
        Stiffness_Matrix;
        Damping_Matrix;
        Inertia_Matrix;
        Gravitational_Vector;
    end
    methods
        function obj = Body(npoints, points, ntetras, tetras)
            obj.numNodalPoints = npoints;
            for k=1:npoints
                pt(k) = three_dim_fea.NodalPoint(points(:,k));
            end
            obj.NodalPoints = pt;
            %
            if ntetras > 0
                obj.numTetrahedrons = ntetras;
                for p=1:ntetras
                    i = tetras(p,1); j = tetras(p,2); k = tetras(p,3); l = tetras(p,4);
                    te(p) = three_dim_fea.Tetrahedron(i, j, k, l, points(:,i), points(:,j), points(:,k), points(:,l));
                end
                obj.Tetrahedrons = te;
            end
            
            obj.J_lambda = zeros(3*obj.numNodalPoints, 3*obj.numNodalPoints);
            obj.J_mu     = zeros(3*obj.numNodalPoints, 3*obj.numNodalPoints);
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                loc = [ 3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j, 3*k-2, 3*k-1, 3*k, 3*l-2, 3*l-1, 3*l ];
                obj.J_lambda(loc,loc) = obj.J_lambda(loc,loc) + tetra.Partial_J_lambda;
                obj.J_mu(loc,loc)     = obj.J_mu(loc,loc)     + tetra.Partial_J_mu;
            end
            
            faces = [];
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                faces = [ faces; tetra.Faces ];
            end
            faces_sorted = sort(faces, 2);
            n = size(faces_sorted, 1);
            onlyone = 1:n;
            for p=1:n-1
                [lia, loc] = ismember(faces_sorted(p,:), faces_sorted(p+1:n,:), 'rows');
                if lia
                    onlyone(p) = 0;
                    onlyone(p+loc) = 0;
                end
            end
            onlyone = setdiff(onlyone, 0);
            obj.BoundaryFaces = faces(onlyone, :);
            
            [ nsufs, surfs, order ] = three_dim_fea.triangles_to_surfaces (obj.BoundaryFaces);
            obj.numSurfaces = nsufs;
            for p=1:nsufs
                [surf_p, area_p] = three_dim_fea.extract_surface(p, surfs, order, obj.BoundaryFaces, points );
                sf(p) = three_dim_fea.Surface(size(surf_p,1), surf_p, area_p);
            end
            obj.Surfaces = sf;
        end
        
        function obj = mechanical_parameters(obj, rho, l, m)
            obj.Density = rho;
            obj.lambda = l;
            obj.mu     = m;
            for p=1:obj.numTetrahedrons
                obj.Tetrahedrons(p).Density = rho;
                obj.Tetrahedrons(p).lambda  = l;
                obj.Tetrahedrons(p).mu      = m;
            end
        end
        function obj = viscous_parameters(obj, lv, mv)
            obj.lambda_vis = lv;
            obj.mu_vis     = mv;
            for p=1:obj.numTetrahedrons
                obj.Tetrahedrons(p).lambda_vis  = lv;
                obj.Tetrahedrons(p).mu_vis      = mv;
            end
        end
  
        
        function energy = total_strain_potential_energy(obj, disps)
            energy = 0;
            for p=1:obj.numTetrahedrons
               tetra = obj.Tetrahedrons(p);
               energy = energy + tetra.partial_strain_potential_energy(disps);
            end
        end
        
        function energy = total_strain_potential_energy_Green_strain(obj, disps)
            energy = 0;
            for p=1:obj.numTetrahedrons
               tetra = obj.Tetrahedrons(p);
               energy = energy + tetra.partial_strain_potential_energy_Green_strain(disps);
            end
        end
        
        function energy = total_gravitational_potential_energy(obj, disps, grav)
            energy = 0;
            for p=1:obj.numTetrahedrons
               tetra = obj.Tetrahedrons(p);
               energy = energy + tetra.partial_gravitational_potential_energy(disps, grav);
            end
        end
        
        function output = draw( obj, disps, color )
            arguments
                obj;
                disps = zeros(3, obj.numNodalPoints);
                color = [0.9, 0.9, 0.9];
            end
            
            points = zeros(3, obj.numNodalPoints);
            np = obj.NodalPoints;
            for p=1:obj.numNodalPoints
                points(:,p) = np(p).Coordinates + disps(:,p);
            end
            
            h = patch('Faces', obj.BoundaryFaces, 'Vertices', points', 'FaceColor', color);
        end
        
       % lzy 2023 10 24 
        function output = draw_individual_app( obj, app_ax, disps, color )
            arguments
                obj;
                app_ax matlab.ui.control.UIAxes; % 必需参数，传递 UIAxes 对象
                disps = zeros(2, obj.numNodalPoints);
                color = [0.9, 0.9, 0.9];
            end
            faces = [1, 2, 3; 1, 2, 4; 1, 3, 4; 2, 3, 4];
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3);l=vs(4);
                np = obj.NodalPoints;
                xi = np(i).Coordinates + disps(:,i);
                xj = np(j).Coordinates + disps(:,j);
                xk = np(k).Coordinates + disps(:,k);
                xl = np(l).Coordinates + disps(:,l);
                points = [ xi, xj, xk, xl ];
                if isempty(obj.Tetrahedrons(p).color)
                    patch(app_ax,'Faces', faces, 'Vertices', points', 'FaceColor', color,'FaceAlpha', 0.5);
                else
                    patch(app_ax,'Faces', faces, 'Vertices', points', 'FaceColor', obj.Tetrahedrons(p).color,'FaceAlpha', 0.5);
                end
            end
        end

          %lzy  2024/8/18
       function output = draw_init_app( obj, app_ax, disps, color)
            arguments
                obj;
                app_ax matlab.ui.control.UIAxes; % 必需参数，传递 UIAxes 对象
                disps = zeros(3, obj.numNodalPoints);
                color = [0.9, 0.9, 0.9];
            end
            faces = [1, 2, 3; 1, 2, 4; 1, 3, 4; 2, 3, 4];
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3);l=vs(4);
                np = obj.NodalPoints;
                xi = np(i).Coordinates + disps(:,i);
                xj = np(j).Coordinates + disps(:,j);
                xk = np(k).Coordinates + disps(:,k);
                xl = np(l).Coordinates + disps(:,l);
                points = [ xi, xj, xk, xl ];
                if isempty(obj.Tetrahedrons(p).color)
                    patch(app_ax,'Faces', faces, 'Vertices', points', 'FaceColor', color,'FaceAlpha', 0.3);
                else
                    patch(app_ax,'Faces', faces, 'Vertices', points', 'FaceColor', obj.Tetrahedrons(p).color,'FaceAlpha', 0.3);
                end
                % calculate the center of tetrahedron
                center_tetra = (xi+xj+xk+xl)./4;
                % draw label
               % htext_in = text(app_ax,center_tetra(1), center_tetra(2), center_tetra(3), num2str(p), ...
               %     'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
               % set(htext_in, 'Clipping', 'off'); % 确保文本不会被遮挡
            end
            %Remove all duplicate values to form a vector
            surface_points_num = unique(obj.Surfaces.Triangles(:), 'stable');
            for i = 1:length(surface_points_num)
                   j = surface_points_num(i);
                   coor_sp = obj.NodalPoints(j).Coordinates;
                   % htext_ex= text(app_ax, coor_sp(1), coor_sp(2), coor_sp(3), num2str(j), ...
                   %    'HorizontalAlignment', 'center', 'FontSize', 8,'Color', 'red');
                   %set(htext_ex, 'Clipping', 'off'); % 确保文本不会被遮挡
            end        
        end
        

        function output = draw_subregion_app( obj, app_ax, disps, color)
            arguments
                obj;
                app_ax matlab.ui.control.UIAxes; % 必需参数，传递 UIAxes 对象
                disps = zeros(3, obj.numNodalPoints);
                color = [0.9, 0.9, 0.9];
            end
            faces = [1, 2, 3; 1, 2, 4; 1, 3, 4; 2, 3, 4];
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3);l=vs(4);
                np = obj.NodalPoints;
                xi = np(i).Coordinates + disps(:,i);
                xj = np(j).Coordinates + disps(:,j);
                xk = np(k).Coordinates + disps(:,k);
                xl = np(l).Coordinates + disps(:,l);
                points = [ xi, xj, xk, xl ];
                if isempty(obj.Tetrahedrons(p).color)
                    patch(app_ax,'Faces', faces, 'Vertices', points', 'FaceColor', color,'FaceAlpha', 0.3);
                else
                    patch(app_ax,'Faces', faces, 'Vertices', points', 'FaceColor', obj.Tetrahedrons(p).color,'FaceAlpha', 0.3);
                end
            end
        end
        
        function obj = calculate_stiffness_matrix(obj)
            obj.Stiffness_Matrix = zeros(3*obj.numNodalPoints, 3*obj.numNodalPoints);
            %
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [obj.Tetrahedrons(p), K_p] = tetra.partial_stiffness_matrix;
                loc = [ 3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j, 3*k-2, 3*k-1, 3*k, 3*l-2, 3*l-1, 3*l ];
                obj.Stiffness_Matrix(loc,loc) = obj.Stiffness_Matrix(loc,loc) + K_p;
            end
            %
        end
        
        function obj = calculate_damping_matrix(obj)
            obj.Damping_Matrix = zeros(3*obj.numNodalPoints, 3*obj.numNodalPoints);
            %
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [obj.Tetrahedrons(p), B_p] = tetra.partial_damping_matrix;
                loc = [ 3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j, 3*k-2, 3*k-1, 3*k, 3*l-2, 3*l-1, 3*l ];
                obj.Damping_Matrix(loc,loc) = obj.Damping_Matrix(loc,loc) + B_p;
            end
            %
        end
        
        function obj = calculate_inertia_matrix(obj)
            obj.Inertia_Matrix = zeros(3*obj.numNodalPoints, 3*obj.numNodalPoints);
            %
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [obj.Tetrahedrons(p), M_p] = tetra.partial_inertia_matrix;
                loc = [ 3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j, 3*k-2, 3*k-1, 3*k, 3*l-2, 3*l-1, 3*l ];
                obj.Inertia_Matrix(loc,loc) = obj.Inertia_Matrix(loc,loc) + M_p;
            end
            %
        end
        
        function obj = calculate_gravitational_vector(obj, g)
            obj.Gravitational_Vector = zeros(3*obj.numNodalPoints, 1);
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [obj.Tetrahedrons(p), grav_p] = tetra.partial_gravitational_vector(g);
                loc = [ 3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j, 3*k-2, 3*k-1, 3*k, 3*l-2, 3*l-1, 3*l ];
                obj.Gravitational_Vector(loc) = obj.Gravitational_Vector(loc) + grav_p;
            end
        end
        
        function A = constraint_matrix(obj, index)
            index = index(index <= obj.numNodalPoints);
            A = zeros(3*obj.numNodalPoints, 3*length(index));
            for k=1:length(index)
                A(3*index(k)-2, 3*k-2) = 1;
                A(3*index(k)-1, 3*k-1) = 1;
                A(3*index(k),   3*k)   = 1;
            end
        end
        
        function [index, index_npoints] = extract_index_for_subregion(obj, index)
            for r=1:obj.numSubRegions
                common = intersect(index, obj.SubRegions(r).Index_Tetrahedrons);
                if common
                    fprintf("%d ", common);
                    fprintf('already in another subregion\n');
                    index = setdiff(index, common);
                end
            end
            index_npoints = [];
            if index
                for p=index
                    tetra = obj.Tetrahedrons(p);
                    vs = tetra.Vertices;
                    index_npoints = union(index_npoints, vs);
                end
            end
        end
        
        function obj = define_subregion(obj, index)
            arguments
                obj; index;
            end
            if isempty(obj.numSubRegions)
                obj.numSubRegions = 0;
            end
            [index, index_npoints] = obj.extract_index_for_subregion(index);
            obj.numSubRegions = obj.numSubRegions + 1;
            obj.SubRegions = [ obj.SubRegions, three_dim_fea.SubRegion(index, index_npoints) ];
            obj = obj.subregion_partial_connection_matrices;
        end
        
        function obj = subregion_partial_connection_matrices(obj)
            r = obj.numSubRegions;
            npoints = obj.SubRegions(r).numNodalPoints;
            pJl = zeros(3*npoints, 3*npoints);
            pJm = zeros(3*npoints, 3*npoints);
            
            list = obj.SubRegions(r).Index_NodalPoints;
            for p=obj.SubRegions(r).Index_Tetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l=vs(4);
                ip = find(list==i); jp = find(list==j); kp = find(list==k);lp = find(list==l);
                loc = [ 3*ip-2, 3*ip-1, 3*ip, 3*jp-2, 3*jp-1, 3*jp, 3*kp-2, 3*kp-1,3*kp, 3*lp-2, 3*lp-1,3*lp];
                pJl(loc,loc) = pJl(loc,loc) + tetra.Partial_J_lambda;
                pJm(loc,loc) = pJm(loc,loc) + tetra.Partial_J_mu;
            end
            
            obj.SubRegions(r).Partial_J_lambda = pJl;
            obj.SubRegions(r).Partial_J_mu     = pJm;
        end
        
        function obj = subregion_mechanical_parameters(obj, rho, l, m, k)
            arguments
                obj; rho; l; m;
                k = obj.numSubRegions;
            end
            obj.SubRegions(k).Density = rho;
            obj.SubRegions(k).lambda  = l;
            obj.SubRegions(k).mu      = m;
            for p = obj.SubRegions(k).Index_Tetrahedrons
                obj.Tetrahedrons(p).Density = rho;
                obj.Tetrahedrons(p).lambda  = l;
                obj.Tetrahedrons(p).mu      = m;
            end
        end
        
        function obj = subregion_viscous_parameters(obj, lv, mv, k)
            arguments
                obj; lv; mv;
                k = obj.numSubRegions;
            end
            obj.SubRegions(k).lambda_vis  = lv;
            obj.SubRegions(k).mu_vis      = mv;
            for p = obj.SubRegions(k).Index_Tetrahedrons
                obj.Tetrahedrons(p).lambda_vis  = lv;
                obj.Tetrahedrons(p).mu_vis      = mv;
            end
        end
        


         function obj = subregion_strain_mea(obj,strain_mea,k)
            arguments
                obj; strain_mea;
                k = obj.numSubRegions;
            end
            obj.SubRegions(k).strain_mea  = strain_mea;
            for p = obj.SubRegions(k).Index_Tetrahedrons
                obj.Tetrahedrons(p).strain_mea  = strain_mea;
            end
         end
         
         function obj = subregion_color(obj, cl, k)
            arguments
                obj; cl;
                k = obj.numSubRegions;
            end
            obj.SubRegions(k).color = cl;
            if obj.SubRegions(k).Index_Tetrahedrons
                for p = obj.SubRegions(k).Index_Tetrahedrons
                    obj.Tetrahedrons(p).color = cl;
                end
            end
        end
        
        function forces = nodal_forces_Cauchy_strain_original(obj, disps)
            forces = zeros(3*obj.numNodalPoints,1);
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [fi, fj, fk, fl] = tetra.nodal_forces_Cauchy_strain(disps(:,i), disps(:,j), disps(:,k), disps(:,l));
                forces(3*i-2:3*i) = forces(3*i-2:3*i) + fi;
                forces(3*j-2:3*j) = forces(3*j-2:3*j) + fj;
                forces(3*k-2:3*k) = forces(3*k-2:3*k) + fk;
                forces(3*l-2:3*l) = forces(3*l-2:3*l) + fl;
            end
        end
        
        function forces = nodal_forces_Cauchy_strain(obj, disps)
            un = reshape(disps, [3*obj.numNodalPoints, 1]);
            forces = - obj.Stiffness_Matrix*un;
        end
        
        function forces = nodal_forces_Green_strain_original(obj, disps)
            forces = zeros(3*obj.numNodalPoints,1);
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [fi, fj, fk, fl] = tetra.nodal_forces_Green_strain(disps(:,i), disps(:,j), disps(:,k), disps(:,l));
                forces(3*i-2:3*i) = forces(3*i-2:3*i) + fi;
                forces(3*j-2:3*j) = forces(3*j-2:3*j) + fj;
                forces(3*k-2:3*k) = forces(3*k-2:3*k) + fk;
                forces(3*l-2:3*l) = forces(3*l-2:3*l) + fl;
            end
        end
        
        function forces = nodal_forces_Green_strain(obj, disps)
            forces = zeros(3*obj.numNodalPoints,1);
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [fi, fj, fk, fl] = tetra.nodal_forces_Green_strain(disps(:,i), disps(:,j), disps(:,k), disps(:,l));
                forces(3*i-2:3*i) = forces(3*i-2:3*i) + fi;
                forces(3*j-2:3*j) = forces(3*j-2:3*j) + fj;
                forces(3*k-2:3*k) = forces(3*k-2:3*k) + fk;
                forces(3*l-2:3*l) = forces(3*l-2:3*l) + fl;
            end
        end
        
         % LZY 2023/10/27  visela forces
        function forces = nodal_total_forces_Green_strain(obj, disps, vels)
            forces = zeros(3*obj.numNodalPoints,1);
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [fi, fj, fk, fl] = tetra.nodal_total_forces_Green_strain(disps(:,i), disps(:,j), disps(:,k), disps(:,l),vels(:,i), vels(:,j), vels(:,k), vels(:,l));
                forces(3*i-2:3*i) = forces(3*i-2:3*i) + fi;
                forces(3*j-2:3*j) = forces(3*j-2:3*j) + fj;
                forces(3*k-2:3*k) = forces(3*k-2:3*k) + fk;
                forces(3*l-2:3*l) = forces(3*l-2:3*l) + fl;
            end
        end
        
        % LZY 2024 7 30
        function forces = nodal_forces_external_stimuli(obj,disps,t)
            forces = zeros(3*obj.numNodalPoints,1);
            for p=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [fi, fj, fk, fl] = tetra.nodal_forces_external_stimuli(disps(:,i), disps(:,j), disps(:,k), disps(:,l),t);
                forces(3*i-2:3*i) = forces(3*i-2:3*i) + fi;
                forces(3*j-2:3*j) = forces(3*j-2:3*j) + fj;
                forces(3*k-2:3*k) = forces(3*k-2:3*k) + fk;
                forces(3*l-2:3*l) = forces(3*l-2:3*l) + fl;
            end
        end
        
        function position = positional_vectors(obj, disps)
            arguments
                obj;
                disps = zeros(3,obj.numNodalPoints);
            end
            position = [];
            for i=1:obj.numNodalPoints
                position = [ position, obj.NodalPoints(i).Coordinates + disps(:,i) ];
            end
        end
        
        function volume = volume_under_deformation(obj, disps)
            arguments
                obj;
                disps = zeros(3,obj.numNodalPoints);
            end
            position = obj.positional_vectors(disps);
            volume = 0;
            for i=1:obj.numTetrahedrons
                tetra = obj.Tetrahedrons(i);
                volume = volume + tetra.volume_under_deformation(position);
            end
        end
    end
end
