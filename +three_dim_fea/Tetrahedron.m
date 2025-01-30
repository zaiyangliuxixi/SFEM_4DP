classdef Tetrahedron
    properties
        Vertices;
        Faces;
        Volume;
        Density; lambda; mu; lambda_vis; mu_vis;
        color;
        vector_a; vector_b; vector_c;
        vector_a2;vector_b2;vector_c2;
        vector_a0;vector_b0;vector_c0;vector_T2;vector_T0;
        vector_a2_hat;vector_b2_hat;vector_c2_hat;
        vector_a0_hat;vector_b0_hat;vector_c0_hat;vector_T2_hat;vector_T0_hat;
        u_x; u_y; u_z; v_x; v_y; v_z; w_x; w_y; w_z;
        strain_mea;
        Cauchy_strain;
        Green_strain;
        Partial_J_lambda; Partial_J_mu;
        Partial_Stiffness_Matrix;
        Partial_Damping_Matrix;
        Partial_Inertia_Matrix;
        Partial_Gravitational_Vector;
        Laa; Lab; Lac; La2; La0;
        Lba; Lbb; Lbc; Lb2; Lb0; 
        Lca; Lcb; Lcc; Lc2; Lc0;
        L2a; L2b; L2c; L22; L20;
        L0a; L0b; L0c; L02; L00;
        Maa; Mab; Mac; Ma2; Ma0; 
        Mba; Mbb; Mbc; Mb2; Mb0; 
        Mca; Mcb; Mcc; Mc2; Mc0;
        M2a; M2b; M2c; M22; M20;
        M0a; M0b; M0c; M02; M00; 
        La; Lb; Lc; L2; L0;
        Ma; Mb; Mc;
        M2a_4DP; M2b_4DP; M2c_4DP;
        M0a_4DP; M0b_4DP; M0c_4DP;
        P_mat_origin;
    end
    methods
        function obj = Tetrahedron(i, j, k, l, pi, pj, pk, pl)
            obj.Vertices = [ i, j, k, l ];
            obj.Faces = [ i,j,k; j,l,k; l,i,k; i,l,j ];
            obj.Volume = det( [ pj-pi, pk-pi, pl-pi ] )/6;
            
            % above equation 2.4.6
            ajkl = (pj(2)*pk(3)-pk(2)*pj(3)) + (pk(2)*pl(3)-pl(2)*pk(3)) + (pl(2)*pj(3)-pj(2)*pl(3));
            akli = (pk(2)*pl(3)-pl(2)*pk(3)) + (pl(2)*pi(3)-pi(2)*pl(3)) + (pi(2)*pk(3)-pk(2)*pi(3));
            alij = (pl(2)*pi(3)-pi(2)*pl(3)) + (pi(2)*pj(3)-pj(2)*pi(3)) + (pj(2)*pl(3)-pl(2)*pj(3));
            aijk = (pi(2)*pj(3)-pj(2)*pi(3)) + (pj(2)*pk(3)-pk(2)*pj(3)) + (pk(2)*pi(3)-pi(2)*pk(3));
            
            bjkl = (pj(3)*pk(1)-pk(3)*pj(1)) + (pk(3)*pl(1)-pl(3)*pk(1)) + (pl(3)*pj(1)-pj(3)*pl(1));
            bkli = (pk(3)*pl(1)-pl(3)*pk(1)) + (pl(3)*pi(1)-pi(3)*pl(1)) + (pi(3)*pk(1)-pk(3)*pi(1));
            blij = (pl(3)*pi(1)-pi(3)*pl(1)) + (pi(3)*pj(1)-pj(3)*pi(1)) + (pj(3)*pl(1)-pl(3)*pj(1));
            bijk = (pi(3)*pj(1)-pj(3)*pi(1)) + (pj(3)*pk(1)-pk(3)*pj(1)) + (pk(3)*pi(1)-pi(3)*pk(1));
            
            cjkl = (pj(1)*pk(2)-pk(1)*pj(2)) + (pk(1)*pl(2)-pl(1)*pk(2)) + (pl(1)*pj(2)-pj(1)*pl(2));
            ckli = (pk(1)*pl(2)-pl(1)*pk(2)) + (pl(1)*pi(2)-pi(1)*pl(2)) + (pi(1)*pk(2)-pk(1)*pi(2));
            clij = (pl(1)*pi(2)-pi(1)*pl(2)) + (pi(1)*pj(2)-pj(1)*pi(2)) + (pj(1)*pl(2)-pl(1)*pj(2));
            cijk = (pi(1)*pj(2)-pj(1)*pi(2)) + (pj(1)*pk(2)-pk(1)*pj(2)) + (pk(1)*pi(2)-pi(1)*pk(2));
            

            
            %lzy 2023/10/27
            obj.vector_a = ( 1/(6*obj.Volume))*[ -ajkl; akli; -alij; aijk ];
            obj.vector_b = ( 1/(6*obj.Volume))*[ -bjkl; bkli; -blij; bijk ];
            obj.vector_c = ( 1/(6*obj.Volume))*[ -cjkl; ckli; -clij; cijk ];
            
            obj.vector_a2 = ( 1/(6*obj.Volume))^2 * [ -ajkl; akli; -alij; aijk ].^2;
            obj.vector_b2 = ( 1/(6*obj.Volume))^2 * [ -bjkl; bkli; -blij; bijk ].^2;
            obj.vector_c2 = ( 1/(6*obj.Volume))^2 * [ -cjkl; ckli; -clij; cijk ].^2;
            
            obj.vector_a0 = 2*( 1/(6*obj.Volume))^2*[-ajkl*akli; ajkl*alij; -ajkl*aijk;...
                                                  -akli*alij; akli*aijk; -alij*aijk];
            obj.vector_b0 = 2*( 1/(6*obj.Volume))^2*[-bjkl*bkli; bjkl*blij; -bjkl*bijk;...
                                                  -bkli*blij; bkli*bijk; -blij*bijk];
            obj.vector_c0 = 2*( 1/(6*obj.Volume))^2*[-cjkl*ckli; cjkl*clij; -cjkl*cijk;...
                                                  -ckli*clij; ckli*cijk; -clij*cijk];

            obj.vector_T2 = obj.vector_a2 + obj.vector_b2 +obj.vector_c2;
            obj.vector_T0 = obj.vector_a0 + obj.vector_b0 +obj.vector_c0;
            
            obj.vector_a2_hat = ( 1/(6*obj.Volume))^2 * [ bjkl*cjkl; bkli*ckli; blij*clij; bijk*cijk];
            obj.vector_b2_hat = ( 1/(6*obj.Volume))^2 * [ ajkl*cjkl; akli*ckli; alij*clij; aijk*cijk];
            obj.vector_c2_hat = ( 1/(6*obj.Volume))^2 * [ ajkl*bjkl; akli*bkli; alij*blij; aijk*bijk];

            obj.vector_a0_hat = ( 1/(6*obj.Volume))^2 * [-bjkl*ckli-cjkl*bkli; bjkl*clij+cjkl*blij; -bjkl*cijk-cjkl*bijk;...
                                                  -bkli*clij-ckli*blij; bkli*cijk+ckli*bijk; -blij*cijk-clij*bijk];
            obj.vector_b0_hat = ( 1/(6*obj.Volume))^2 * [-cjkl*akli-ajkl*ckli; cjkl*alij+ajkl*clij; -cjkl*aijk-ajkl*cijk;...
                                                  -ckli*alij-akli*clij; ckli*aijk+akli*cijk; -clij*aijk-alij*cijk];
            obj.vector_c0_hat = ( 1/(6*obj.Volume))^2 * [-ajkl*bkli-bjkl*akli; ajkl*blij+bjkl*alij; -ajkl*bijk-bjkl*aijk;...
                                                  -akli*blij-bkli*alij; akli*bijk+bkli*aijk; -alij*bijk-blij*aijk];

            a2 =  obj.vector_a2; b2 = obj.vector_b2; c2 = obj.vector_c2; 
            a0 =  obj.vector_a0; b0 = obj.vector_b0; c0 = obj.vector_c0;
            T2 =  obj.vector_T2; T0 = obj.vector_T0;
            a2_hat =  obj.vector_a2_hat; b2_hat = obj.vector_b2_hat; c2_hat = obj.vector_c2_hat; 
            a0_hat =  obj.vector_a0_hat; b0_hat = obj.vector_b0_hat; c0_hat = obj.vector_c0_hat;

              
            a = obj.vector_a;
            b = obj.vector_b;
            c = obj.vector_c;
            
            H_lambda = [ a*a', a*b', a*c';
                  b*a', b*b', b*c';
                  c*a', c*b', c*c' ].*obj.Volume;
            H_mu = [ 2*a*a'+b*b'+c*c', b*a', c*a';
                  a*b', 2*b*b'+c*c'+a*a', c*b';
                  a*c', b*c', 2*c*c'+a*a'+b*b' ].*obj.Volume;
            
            enum = [ 1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12 ];
            obj.Partial_J_lambda = H_lambda(enum, enum);
            obj.Partial_J_mu     = H_mu(enum, enum);
            %LZY  2023 10 27
            obj.Laa = obj.Volume.*a*transpose(a);
            obj.Lab = obj.Volume.*a*transpose(b);
            obj.Lac = obj.Volume.*a*transpose(c);
            obj.La2 = 0.5*obj.Volume.*a*transpose(T2);
            obj.La0 = 0.5*obj.Volume.*a*transpose(T0);

            obj.Lba = transpose(obj.Lab);
            obj.Lbb = obj.Volume.*b*transpose(b);
            obj.Lbc = obj.Volume.*b*transpose(c);
            obj.Lb2 = 0.5*obj.Volume.*b*transpose(T2);
            obj.Lb0 = 0.5*obj.Volume.*b*transpose(T0);

            obj.Lca = transpose(obj.Lac);
            obj.Lcb = transpose(obj.Lbc);
            obj.Lcc = obj.Volume.*c*transpose(c);
            obj.Lc2 = 0.5*obj.Volume.*c*transpose(T2);
            obj.Lc0 = 0.5*obj.Volume.*c*transpose(T0);

            obj.L2a = transpose(obj.La2);
            obj.L2b = transpose(obj.Lb2);
            obj.L2c = transpose(obj.Lc2);
            obj.L22 = 0.25*obj.Volume.*T2*transpose(T2);
            obj.L20 = 0.25*obj.Volume.*T2*transpose(T0);

            obj.L0a = transpose(obj.La0);
            obj.L0b = transpose(obj.Lb0);
            obj.L0c = transpose(obj.Lc0);
            obj.L02 = transpose(obj.L20);
            obj.L00 = 0.25*obj.Volume.*T0*transpose(T0);

            obj.La = obj.Volume.*a;   % row vector
            obj.Lb = obj.Volume.*b;
            obj.Lc = obj.Volume.*c;
            obj.L2 = 0.5* obj.Volume.*T2;
            obj.L0 = 0.5* obj.Volume.*T0;

            obj.Maa = 2*obj.Laa + obj.Lbb+ obj.Lcc;
            obj.Mab = obj.Lba;
            obj.Mac = obj.Lca;
            obj.Ma2 = obj.Volume.*(a*transpose(a2)+ b*transpose(c2_hat) + c*transpose(b2_hat));
            obj.Ma0 = obj.Volume.*(a*transpose(a0)+ b*transpose(c0_hat) + c*transpose(b0_hat));

            obj.Mba = transpose(obj.Mab);
            obj.Mbb = 2*obj.Lbb +obj.Lcc+ obj.Laa;
            obj.Mbc = obj.Lcb; 
            obj.Mb2 = obj.Volume.*(b*transpose(b2)+ c*transpose(a2_hat)+ a*transpose(c2_hat));  
            obj.Mb0 = obj.Volume.*(b*transpose(b0)+  c*transpose(a0_hat)+ a*transpose(c0_hat));
            
            obj.Mca = transpose(obj.Mac);
            obj.Mcb = transpose(obj.Mbc);
            obj.Mcc = 2*obj.Lcc +obj.Laa+ obj.Lbb;
            obj.Mc2 = obj.Volume.*(c*transpose(c2)+ a*transpose(b2_hat) + b*transpose(a2_hat));
            obj.Mc0 = obj.Volume.*(c*transpose(c0)+ a*transpose(b0_hat) + b*transpose(a0_hat));
            
            obj.M2a = transpose(obj.Ma2);
            obj.M2b = transpose(obj.Mb2);
            obj.M2c = transpose(obj.Mc2);
            obj.M22 = obj.Volume.*(0.5.*(a2*transpose(a2)+ b2*transpose(b2) + c2*transpose(c2))...
                + (a2_hat*transpose(a2_hat)+ b2_hat*transpose(b2_hat) + c2_hat*transpose(c2_hat)));
            obj.M20 = obj.Volume.*(0.5.*(a2*transpose(a0)+ b2*transpose(b0) + c2*transpose(c0))...
                + (a2_hat*transpose(a0_hat)+ b2_hat*transpose(b0_hat) + c2_hat*transpose(c0_hat)));

            obj.M0a = transpose(obj.Ma0);
            obj.M0b = transpose(obj.Mb0);
            obj.M0c = transpose(obj.Mc0);
            obj.M02 = transpose(obj.M20);
            obj.M00 = obj.Volume.*(0.5.*(a0*transpose(a0)+ b0*transpose(b0) + c0*transpose(c0))...
                + (a0_hat*transpose(a0_hat)+ b0_hat*transpose(b0_hat) + c0_hat*transpose(c0_hat)));

            obj.Ma = 2*obj.La;
            obj.Mb = 2*obj.Lb;
            obj.Mc = 2*obj.Lc;
            obj.M2a_4DP = obj.Volume.*a2;
            obj.M2b_4DP = obj.Volume.*b2;
            obj.M2c_4DP = obj.Volume.*c2;
            obj.M0a_4DP = obj.Volume.*a0;
            obj.M0b_4DP = obj.Volume.*b0;
            obj.M0c_4DP = obj.Volume.*c0; 

            P_mat = zeros(12,12);
            P_mat(1,1) = 1;P_mat(2,5) = 1;P_mat(3,9) = 1;
            P_mat(4,2) = 1;P_mat(5,6) = 1;P_mat(6,10) = 1;
            P_mat(7,3) = 1;P_mat(8,7) = 1;P_mat(9,11) = 1;
            P_mat(10,4) = 1;P_mat(11,8) = 1;P_mat(12,12) = 1;
            obj.P_mat_origin = P_mat;
        end
                        % 2025 1 28
        function obj = update_mass(obj, pi, pj, pk, pl)  
            obj.Volume = det( [ pj-pi, pk-pi, pl-pi ] )/6;
        end
        
        function obj = mechanical_parameters(obj, rho, lambda, mu)
            obj.Density = rho;
            obj.lambda = lambda; obj.mu = mu;
        end
        
        function obj = partial_derivaties(obj, ui, uj, uk, ul)
            gamma_u =  [ ui(1); uj(1); uk(1); ul(1) ];
            gamma_v =  [ ui(2); uj(2); uk(2); ul(2) ];
            gamma_w =  [ ui(3); uj(3); uk(3); ul(3) ];
            obj.u_x = obj.vector_a' * gamma_u;
            obj.u_y = obj.vector_b' * gamma_u;
            obj.u_z = obj.vector_c' * gamma_u;
            obj.v_x = obj.vector_a' * gamma_v;
            obj.v_y = obj.vector_b' * gamma_v;
            obj.v_z = obj.vector_c' * gamma_v;
            obj.w_x = obj.vector_a' * gamma_w;
            obj.w_y = obj.vector_b' * gamma_w;
            obj.w_z = obj.vector_c' * gamma_w;
        end
        
        function obj = calculate_Cauchy_strain(obj, ui, uj, uk, ul)
            obj = obj.partial_derivaties (ui, uj, uk, ul);
            obj.Cauchy_strain = [
                obj.u_x; obj.v_y;  obj.w_z;
                obj.v_z + obj.w_y;
                obj.w_x + obj.u_z;
                obj.u_y + obj.v_x ];
        end
        
        function energy = partial_strain_potential_energy(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3)); ul = disps(:,vs(4));
            obj = obj.calculate_Cauchy_strain(ui, uj, uk, ul);
            energy = (1/2) * (obj.Volume) * ...
                ( obj.lambda * ( obj.Cauchy_strain(1) + obj.Cauchy_strain(2) + obj.Cauchy_strain(3) )^2 + ...
                  obj.mu * ( 2*obj.Cauchy_strain(1)^2 + 2*obj.Cauchy_strain(2)^2 + 2*obj.Cauchy_strain(3)^2 + ...
                               obj.Cauchy_strain(4)^2 +   obj.Cauchy_strain(5)^2 +   obj.Cauchy_strain(6)^2 ) );
        end

        function obj = calculate_Green_strain(obj, ui, uj, uk, ul)
            obj = obj.partial_derivaties (ui, uj, uk, ul);
            obj.Green_strain = [ ...
                obj.u_x + (1/2)*(obj.u_x^2 + obj.v_x^2 + obj.w_x^2); ...
                obj.v_y + (1/2)*(obj.u_y^2 + obj.v_y^2 + obj.w_y^2); ...
                obj.w_z + (1/2)*(obj.u_z^2 + obj.v_z^2 + obj.w_z^2); ...
                obj.v_z + obj.w_y + (obj.u_y*obj.u_z + obj.v_y*obj.v_z + obj.w_y*obj.w_z); ...
                obj.w_x + obj.u_z + (obj.u_z*obj.u_x + obj.v_z*obj.v_x + obj.w_z*obj.w_x); ...
                obj.u_y + obj.v_x + (obj.u_x*obj.u_y + obj.v_x*obj.v_y + obj.w_x*obj.w_y) ];
        end
        
        function energy = partial_strain_potential_energy_Green_strain(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3)); ul = disps(:,vs(4));
            obj = obj.calculate_Green_strain(ui, uj, uk, ul);
            energy = (1/2) * (obj.Volume) * ...
                ( obj.lambda * ( obj.Green_strain(1) + obj.Green_strain(2) + obj.Green_strain(3) )^2 + ...
                  obj.mu * ( 2*obj.Green_strain(1)^2 + 2*obj.Green_strain(2)^2 + 2*obj.Green_strain(3)^2 + ...
                               obj.Green_strain(4)^2 +   obj.Green_strain(5)^2 +   obj.Green_strain(6)^2 ) );
        end
        
        function energy = partial_gravitational_potential_energy(obj, disps, grav)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3)); ul = disps(:,vs(4));
            energy = - (obj.Volume) * grav' *(ui+uj+uk+ul)/4;
        end
        
        function obj = calculate_partial_stiffness_matrix(obj)            
            obj.Partial_Stiffness_Matrix = ...
                obj.lambda * obj.Partial_J_lambda + obj.mu * obj.Partial_J_mu;
        end
        
        function [obj, K_p] = partial_stiffness_matrix(obj)
            if isempty( obj.Partial_Stiffness_Matrix )
                obj = obj.calculate_partial_stiffness_matrix;
            end
            K_p = obj.Partial_Stiffness_Matrix;
        end
        
        function obj = calculate_partial_damping_matrix(obj)
            obj.Partial_Damping_Matrix = ...
                obj.lambda_vis * obj.Partial_J_lambda + obj.mu_vis * obj.Partial_J_mu;
        end
        
        function [obj, B_p] = partial_damping_matrix(obj)
            if isempty( obj.Partial_Damping_Matrix )
                obj = obj.calculate_partial_damping_matrix;
            end
            B_p = obj.Partial_Damping_Matrix;
        end
        
        function obj = calculate_partial_inertia_matrix(obj)
            mass = obj.Density * obj.Volume;
            %I = eye(3);
            %obj.Partial_Inertia_Matrix = (mass/20) * ...
            %    [2*I, I, I, I; I, 2*I, I, I; I, I, 2*I, I; I, I, I, 2*I];
            obj.Partial_Inertia_Matrix = (mass/4) * eye(12);
        end
        
        function [obj, M_p] = partial_inertia_matrix(obj)
            %if isempty( obj.Partial_Inertia_Matrix )
            %    obj = obj.calculate_partial_inertia_matrix;
            %end
            obj = obj.calculate_partial_inertia_matrix;
            M_p = obj.Partial_Inertia_Matrix;
        end
        
        function obj = calculate_partial_gravitational_vector(obj, g)
            mass = obj.Density * obj.Volume;
            obj.Partial_Gravitational_Vector = mass/4 * [ g; g; g; g ]; 
        end
        
        function [obj, grav_p] = partial_gravitational_vector(obj, g)
            if isempty( obj.Partial_Gravitational_Vector )
                obj = obj.calculate_partial_gravitational_vector(g);
            end
            grav_p = obj.Partial_Gravitational_Vector;
        end
        
        function [fi, fj, fk, fl] = nodal_forces_Cauchy_strain(obj, ui, uj, uk, ul)
            a = obj.vector_a;
            b = obj.vector_b;
            c = obj.vector_c;
            obj = obj.calculate_Cauchy_strain (ui, uj, uk, ul);
            mat = obj.lambda*[ones(3,3), zeros(3,3); zeros(3,3), zeros(3,3)] + ...
                      obj.mu*diag([2, 2, 2, 1, 1, 1]);
            Up_e = mat*obj.Cauchy_strain*obj.Volume;
            Up_gammau = [ a, zeros(4,1), zeros(4,1), zeros(4,1), c, b ]*Up_e;
            Up_gammav = [ zeros(4,1), b, zeros(4,1), c, zeros(4,1), a ]*Up_e;
            Up_gammaw = [ zeros(4,1), zeros(4,1), c, b, a, zeros(4,1) ]*Up_e;
            fi = - [ Up_gammau(1); Up_gammav(1); Up_gammaw(1) ];
            fj = - [ Up_gammau(2); Up_gammav(2); Up_gammaw(2) ];
            fk = - [ Up_gammau(3); Up_gammav(3); Up_gammaw(3) ];
            fl = - [ Up_gammau(4); Up_gammav(4); Up_gammaw(4) ];
        end
        
        function [fi, fj, fk, fl] = nodal_forces_Green_strain(obj, ui, uj, uk, ul)
            a = obj.vector_a;
            b = obj.vector_b;
            c = obj.vector_c;
            obj = obj.calculate_Green_strain (ui, uj, uk, ul);
            mat = obj.lambda*[ones(3,3), zeros(3,3); zeros(3,3), zeros(3,3)] + ...
                      obj.mu*diag([2, 2, 2, 1, 1, 1]);
            Up_E = mat*obj.Green_strain*obj.Volume;
            ux = obj.u_x; uy = obj.u_y; uz = obj.u_z;
            vx = obj.v_x; vy = obj.v_y; vz = obj.v_z;
            wx = obj.w_x; wy = obj.w_y; wz = obj.w_z;
            e_gammau = [ 1+ux,  0,  0,  0,   uz,   uy ].*a + ...
                       [    0, uy,  0, uz,    0, 1+ux ].*b + ...
                       [    0,  0, uz, uy, 1+ux,    0 ].*c;
            e_gammav = [ vx,    0,  0,    0, vz, 1+vy ].*a + ...
                       [  0, 1+vy,  0,   vz,  0,   vx ].*b + ...
                       [  0,    0, vz, 1+vy, vx,    0 ].*c;
            e_gammaw = [ wx,  0,    0,    0, 1+wz, wy ].*a + ...
                       [  0, wy,    0, 1+wz,    0, wx ].*b + ...
                       [  0,  0, 1+wz,   wy,   wx,  0 ].*c;
            Up_gammau = e_gammau*Up_E;
            Up_gammav = e_gammav*Up_E;
            Up_gammaw = e_gammaw*Up_E;
            fi = - [ Up_gammau(1); Up_gammav(1); Up_gammaw(1) ];
            fj = - [ Up_gammau(2); Up_gammav(2); Up_gammaw(2) ];
            fk = - [ Up_gammau(3); Up_gammav(3); Up_gammaw(3) ];
            fl = - [ Up_gammau(4); Up_gammav(4); Up_gammaw(4) ];
        end
       
        %LZY 2023 10 27  total visela forces
        function [fi, fj, fk,fl] = nodal_total_forces_Green_strain(obj, dispi, dispj, dispk,displ, veli, velj, velk,vell)
            ui = dispi(1);uj = dispj(1);uk = dispk(1);ul = displ(1);
            vi = dispi(2);vj = dispj(2);vk = dispk(2);vl = displ(2);
            wi = dispi(3);wj = dispj(3);wk = dispk(3);wl = displ(3);
            ui_dot = veli(1);uj_dot = velj(1); uk_dot = velk(1);ul_dot = vell(1);
            vi_dot = veli(2);vj_dot = velj(2); vk_dot = velk(2);vl_dot = vell(2);
            wi_dot = veli(3);wj_dot = velj(3); wk_dot = velk(3);wl_dot = vell(3);  
            
            gamma_u = [ui;uj;uk;ul];
            gamma_v = [vi;vj;vk;vl];
            gamma_w = [wi;wj;wk;wl];
            gamma_u2 = [ui^2;uj^2;uk^2;ul^2];
            gamma_v2 = [vi^2;vj^2;vk^2;vl^2];
            gamma_w2 = [wi^2;wj^2;wk^2;wl^2];
            gamma_uu = [ui*uj;ui*uk;ui*ul;uj*uk;uj*ul;uk*ul];
            gamma_vv = [vi*vj;vi*vk;vi*vl;vj*vk;vj*vl;vk*vl];
            gamma_ww = [wi*wj;wi*wk;wi*wl;wj*wk;wj*wl;wk*wl];%equation 9.39
            gamma_2 = gamma_u2 + gamma_v2 + gamma_w2;
            gamma_0 = gamma_uu + gamma_vv + gamma_ww;%equation 9.46

            gamma_u_dot = [ui_dot;uj_dot;uk_dot;ul_dot];
            gamma_v_dot = [vi_dot;vj_dot;vk_dot;vl_dot];
            gamma_w_dot = [wi_dot;wj_dot;wk_dot;wl_dot];
            gamma_u2_dot = [2*ui*ui_dot;2*uj*uj_dot;2*uk*uk_dot;2*ul*ul_dot];
            gamma_v2_dot = [2*vi*vi_dot;2*vj*vj_dot;2*vk*vk_dot;2*vl*vl_dot];
            gamma_w2_dot = [2*wi*wi_dot;2*wj*wj_dot;2*wk*wk_dot;2*wl*wl_dot];
            gamma_2_dot = gamma_u2_dot + gamma_v2_dot + gamma_w2_dot;
            gamma_uu_dot = [ui_dot*uj + ui*uj_dot;ui_dot*uk + ui*uk_dot;ui_dot*ul + ui*ul_dot;uj_dot*uk + uj*uk_dot;uj_dot*ul + uj*ul_dot;uk_dot*ul + uk*ul_dot];
            gamma_vv_dot = [vi_dot*vj + vi*vj_dot;vi_dot*vk + vi*vk_dot;vi_dot*vl + vi*vl_dot;vj_dot*vk + vj*vk_dot;vj_dot*vl + vj*vl_dot;vk_dot*vl + vk*vl_dot];
            gamma_ww_dot = [wi_dot*wj + wi*wj_dot;wi_dot*wk + wi*wk_dot;wi_dot*wl + wi*wl_dot;wj_dot*wk + wj*wk_dot;wj_dot*wl + wj*wl_dot;wk_dot*wl + wk*wl_dot];
            gamma_0_dot = gamma_uu_dot + gamma_vv_dot + gamma_ww_dot;

            gamma2_gamma_u = diag([2*ui;2*uj;2*uk;2*ul]);  % equation 9.50下面
            gamma0_gamma_u = [uj uk ul 0 0 0;ui 0 0 uk ul 0;0 ui 0 uj 0 ul;0 0 ui 0 uj uk];
            gamma2_gamma_v = diag([2*vi;2*vj;2*vk;2*vl]);
            gamma0_gamma_v = [vj vk vl 0 0 0;vi 0 0 vk vl 0;0 vi 0 vj 0 vl;0 0 vi 0 vj vk];
            gamma2_gamma_w = diag([2*wi;2*wj;2*wk;2*wl]);
            gamma0_gamma_w = [wj wk wl 0 0 0;wi 0 0 wk wl 0;0 wi 0 wj 0 wl;0 0 wi 0 wj wk];

            gamma2_gamma_u_dot = diag([2*ui_dot;2*uj_dot;2*uk_dot;2*ul_dot]);
            gamma0_gamma_u_dot = [uj_dot uk_dot ul_dot 0 0 0;ui_dot 0 0 uk_dot ul_dot 0;0 ui_dot 0 uj_dot 0 ul_dot;0 0 ui_dot 0 uj_dot uk_dot];
            gamma2_gamma_v_dot = diag([2*vi_dot;2*vj_dot;2*vk_dot;2*vl_dot]);
            gamma0_gamma_v_dot = [vj_dot vk_dot vl_dot 0 0 0;vi_dot 0 0 vk_dot vl_dot 0;0 vi_dot 0 vj_dot 0 vl_dot;0 0 vi_dot 0 vj_dot vk_dot];
            gamma2_gamma_w_dot = diag([2*wi_dot;2*wj_dot;2*wk_dot;2*wl_dot]);
            gamma0_gamma_w_dot = [wj_dot wk_dot wl_dot 0 0 0;wi_dot 0 0 wk_dot wl_dot 0;0 wi_dot 0 wj_dot 0 wl_dot;0 0 wi_dot 0 wj_dot wk_dot];


            g_u_lam = obj.Laa*gamma_u + obj.Lab*gamma_v + obj.Lac*gamma_w + obj.La2*gamma_2 + obj.La0*gamma_0;
            g_v_lam = obj.Lba*gamma_u + obj.Lbb*gamma_v + obj.Lbc*gamma_w + obj.Lb2*gamma_2 + obj.Lb0*gamma_0;
            g_w_lam = obj.Lca*gamma_u + obj.Lcb*gamma_v + obj.Lcc*gamma_w + obj.Lc2*gamma_2 + obj.Lc0*gamma_0;
            
            g_2_lam = obj.L2a*gamma_u + obj.L2b*gamma_v + obj.L2c*gamma_w + obj.L22*gamma_2 + obj.L20*gamma_0;
            g_0_lam = obj.L0a*gamma_u + obj.L0b*gamma_v + obj.L0c*gamma_w + obj.L02*gamma_2 + obj.L00*gamma_0;
    
            g_u_mu = obj.Maa*gamma_u + obj.Mab*gamma_v + obj.Mac*gamma_w + obj.Ma2*gamma_2 + obj.Ma0*gamma_0;
            g_v_mu = obj.Mba*gamma_u + obj.Mbb*gamma_v + obj.Mbc*gamma_w + obj.Mb2*gamma_2 + obj.Mb0*gamma_0;
            g_w_mu = obj.Mca*gamma_u + obj.Mcb*gamma_v + obj.Mcc*gamma_w + obj.Mc2*gamma_2 + obj.Mc0*gamma_0;
            
            g_2_mu = obj.M2a*gamma_u + obj.M2b*gamma_v + obj.M2c*gamma_w + obj.M22*gamma_2 + obj.M20*gamma_0;
            g_0_mu = obj.M0a*gamma_u + obj.M0b*gamma_v + obj.M0c*gamma_w + obj.M02*gamma_2 + obj.M00*gamma_0;


            g_u_lam_dot = obj.Laa*gamma_u_dot + obj.Lab*gamma_v_dot + obj.Lac*gamma_w_dot + obj.La2*gamma_2_dot + obj.La0*gamma_0_dot;
            g_v_lam_dot = obj.Lba*gamma_u_dot + obj.Lbb*gamma_v_dot + obj.Lbc*gamma_w_dot + obj.Lb2*gamma_2_dot + obj.Lb0*gamma_0_dot;
            g_w_lam_dot = obj.Lca*gamma_u_dot + obj.Lcb*gamma_v_dot + obj.Lcc*gamma_w_dot + obj.Lc2*gamma_2_dot + obj.Lc0*gamma_0_dot;
            g_2_lam_dot = obj.L2a*gamma_u_dot + obj.L2b*gamma_v_dot + obj.L2c*gamma_w_dot + obj.L22*gamma_2_dot + obj.L20*gamma_0_dot;
            g_0_lam_dot = obj.L0a*gamma_u_dot + obj.L0b*gamma_v_dot + obj.L0c*gamma_w_dot + obj.L02*gamma_2_dot + obj.L00*gamma_0_dot;
    
            g_u_mu_dot = obj.Maa*gamma_u_dot + obj.Mab*gamma_v_dot + obj.Mac*gamma_w_dot + obj.Ma2*gamma_2_dot + obj.Ma0*gamma_0_dot;
            g_v_mu_dot = obj.Mba*gamma_u_dot + obj.Mbb*gamma_v_dot + obj.Mbc*gamma_w_dot + obj.Mb2*gamma_2_dot + obj.Mb0*gamma_0_dot;
            g_w_mu_dot = obj.Mca*gamma_u_dot + obj.Mcb*gamma_v_dot + obj.Mcc*gamma_w_dot + obj.Mc2*gamma_2_dot + obj.Mc0*gamma_0_dot;
            g_2_mu_dot = obj.M2a*gamma_u_dot + obj.M2b*gamma_v_dot + obj.M2c*gamma_w_dot + obj.M22*gamma_2_dot + obj.M20*gamma_0_dot;
            g_0_mu_dot = obj.M0a*gamma_u_dot + obj.M0b*gamma_v_dot + obj.M0c*gamma_w_dot + obj.M02*gamma_2_dot + obj.M00*gamma_0_dot;


            G_lam_gamma_u = g_u_lam + gamma2_gamma_u*g_2_lam + gamma0_gamma_u*g_0_lam;   %Prof Wang 12
            G_lam_gamma_v = g_v_lam + gamma2_gamma_v*g_2_lam + gamma0_gamma_v*g_0_lam;
            G_lam_gamma_w = g_w_lam + gamma2_gamma_w*g_2_lam + gamma0_gamma_w*g_0_lam;

            G_mu_gamma_u = g_u_mu + gamma2_gamma_u*g_2_mu + gamma0_gamma_u*g_0_mu;
            G_mu_gamma_v = g_v_mu + gamma2_gamma_v*g_2_mu + gamma0_gamma_v*g_0_mu;
            G_mu_gamma_w = g_w_mu + gamma2_gamma_w*g_2_mu + gamma0_gamma_w*g_0_mu;
            
            G_lam_gamma_u_dot = g_u_lam_dot + gamma2_gamma_u_dot*g_2_lam + gamma2_gamma_u*g_2_lam_dot + gamma0_gamma_u_dot*g_0_lam + gamma0_gamma_u*g_0_lam_dot;
            G_lam_gamma_v_dot = g_v_lam_dot + gamma2_gamma_v_dot*g_2_lam + gamma2_gamma_v*g_2_lam_dot + gamma0_gamma_v_dot*g_0_lam + gamma0_gamma_v*g_0_lam_dot;
            G_lam_gamma_w_dot = g_w_lam_dot + gamma2_gamma_w_dot*g_2_lam + gamma2_gamma_w*g_2_lam_dot + gamma0_gamma_w_dot*g_0_lam + gamma0_gamma_w*g_0_lam_dot;

            G_mu_gamma_u_dot = g_u_mu_dot + gamma2_gamma_u_dot*g_2_mu + gamma2_gamma_u*g_2_mu_dot + gamma0_gamma_u_dot*g_0_mu + gamma0_gamma_u*g_0_mu_dot;
            G_mu_gamma_v_dot = g_v_mu_dot + gamma2_gamma_v_dot*g_2_mu + gamma2_gamma_v*g_2_mu_dot + gamma0_gamma_v_dot*g_0_mu + gamma0_gamma_v*g_0_mu_dot;
            G_mu_gamma_w_dot = g_w_mu_dot + gamma2_gamma_w_dot*g_2_mu + gamma2_gamma_w*g_2_mu_dot + gamma0_gamma_w_dot*g_0_mu + gamma0_gamma_w*g_0_mu_dot;


            F_vis = obj.P_mat_origin*[obj.lambda_vis.*G_lam_gamma_u_dot + obj.mu_vis.*G_mu_gamma_u_dot;...
                     obj.lambda_vis.*G_lam_gamma_v_dot + obj.mu_vis.*G_mu_gamma_v_dot;...
                     obj.lambda_vis.*G_lam_gamma_w_dot + obj.mu_vis.*G_mu_gamma_w_dot];
        
            F_ela =  obj.P_mat_origin*[obj.lambda*G_lam_gamma_u + obj.mu*G_mu_gamma_u;...
                obj.lambda*G_lam_gamma_v + obj.mu*G_mu_gamma_v;...
                obj.lambda*G_lam_gamma_w + obj.mu*G_mu_gamma_w];

            fi = - [  F_vis(1)+ F_ela(1);  F_vis(2)+ F_ela(2);  F_vis(3)+ F_ela(3) ];
            fj = - [  F_vis(4)+ F_ela(4);  F_vis(5)+ F_ela(5);  F_vis(6)+ F_ela(6) ];
            fk = - [  F_vis(7)+ F_ela(7);  F_vis(8)+ F_ela(8);  F_vis(9)+ F_ela(9) ];
            fl = - [  F_vis(10)+ F_ela(10); F_vis(11)+ F_ela(11); F_vis(12)+ F_ela(12) ];
        end
       
        
         %LZY  2024 8 1       external stimuli forces
         function [fi, fj, fk,fl] = nodal_forces_external_stimuli(obj,dispi, dispj, dispk,displ,t)
            ui = dispi(1);uj = dispj(1);uk = dispk(1);ul = displ(1);
            vi = dispi(2);vj = dispj(2);vk = dispk(2);vl = displ(2);
            wi = dispi(3);wj = dispj(3);wk = dispk(3);wl = displ(3);
             
            mea = obj.strain_mea(t);
            alpha = mea + 0.5*mea.*mea;

            gamma2_gamma_u = diag([2*ui;2*uj;2*uk;2*ul]);  % equation 9.50下面
            gamma0_gamma_u = [uj uk ul 0 0 0;ui 0 0 uk ul 0;0 ui 0 uj 0 ul;0 0 ui 0 uj uk];
            gamma2_gamma_v = diag([2*vi;2*vj;2*vk;2*vl]);
            gamma0_gamma_v = [vj vk vl 0 0 0;vi 0 0 vk vl 0;0 vi 0 vj 0 vl;0 0 vi 0 vj vk];
            gamma2_gamma_w = diag([2*wi;2*wj;2*wk;2*wl]);
            gamma0_gamma_w = [wj wk wl 0 0 0;wi 0 0 wk wl 0;0 wi 0 wj 0 wl;0 0 wi 0 wj wk];

            g_u_lam_swell = (alpha(1) + alpha(2)+alpha(3)).*obj.La; %Prof Wang资料 12下面
            g_v_lam_swell = (alpha(1) + alpha(2)+alpha(3)).*obj.Lb;
            g_w_lam_swell = (alpha(1) + alpha(2)+alpha(3)).*obj.Lc;
            g_2_lam_swell = (alpha(1) + alpha(2)+alpha(3)).*obj.L2;
            g_0_lam_swell = (alpha(1) + alpha(2)+alpha(3)).*obj.L0;
            
            g_u_mu_swell = alpha(1).*obj.Ma;
            g_v_mu_swell = alpha(2).*obj.Mb;
            g_w_mu_swell = alpha(3).*obj.Mc;
            g_2_mu_swell = alpha(1).*obj.M2a_4DP + alpha(2).*obj.M2b_4DP+ alpha(3).*obj.M2c_4DP;
            g_0_mu_swell = alpha(1).*obj.M0a_4DP + alpha(2).*obj.M0b_4DP+ alpha(3).*obj.M0c_4DP;
            
            G_lam_gamma_u_swell = g_u_lam_swell + gamma2_gamma_u*g_2_lam_swell + gamma0_gamma_u*g_0_lam_swell;
            G_lam_gamma_v_swell = g_v_lam_swell + gamma2_gamma_v*g_2_lam_swell + gamma0_gamma_v*g_0_lam_swell;
            G_lam_gamma_w_swell = g_w_lam_swell + gamma2_gamma_w*g_2_lam_swell + gamma0_gamma_w*g_0_lam_swell;
            G_mu_gamma_u_swell = g_u_mu_swell + gamma2_gamma_u*g_2_mu_swell + gamma0_gamma_u*g_0_mu_swell;
            G_mu_gamma_v_swell = g_v_mu_swell + gamma2_gamma_v*g_2_mu_swell + gamma0_gamma_v*g_0_mu_swell;
            G_mu_gamma_w_swell = g_w_mu_swell + gamma2_gamma_w*g_2_mu_swell + gamma0_gamma_w*g_0_mu_swell;

            F_swell = obj.P_mat_origin*[obj.lambda*G_lam_gamma_u_swell + obj.mu*G_mu_gamma_u_swell;...
                         obj.lambda*G_lam_gamma_v_swell + obj.mu*G_mu_gamma_v_swell;...
                         obj.lambda*G_lam_gamma_w_swell + obj.mu*G_mu_gamma_w_swell];
                            
            fi =  [  F_swell(1);  F_swell(2);  F_swell(3) ];
            fj =  [  F_swell(4);  F_swell(5);  F_swell(6) ];
            fk =  [  F_swell(7);  F_swell(8);  F_swell(9) ];
            fl =  [  F_swell(10); F_swell(11); F_swell(12) ];
         end
         
        function volume = volume_under_deformation(obj, position)
            vs = obj.Vertices;
            pi = position(:,vs(1));
            pj = position(:,vs(2));
            pk = position(:,vs(3));
            pl = position(:,vs(4));
            volume = (1/6)*det([ pj-pi, pk-pi, pl-pi ]);
        end
    end
end
