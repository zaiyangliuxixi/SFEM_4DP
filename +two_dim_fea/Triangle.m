classdef Triangle
    properties
        Vertices;
        Area;
        Thickness;
        Density; lambda; mu; lambda_vis; mu_vis;
        strain_mea; x_i;x_j;x_k;y_i;y_j;y_k;
        color;
        vector_a; vector_b;vector_a2;vector_b2;vector_c2;vector_a0;vector_b0;vector_c0;vector_T2;vector_T0;
        u_x; u_y; v_x; v_y;
        Cauchy_strain;
        Green_strain;
        Partial_J_lambda; Partial_J_mu;
        Partial_Stiffness_Matrix;
        Partial_Damping_Matrix;
        Partial_Inertia_Matrix;
        Partial_Gravitational_Vector;
        Laa; Lab; La2; La0; 
        Lba; Lbb; Lb2; Lb0; 
        L2a; L2b; L22; L20; 
        L0a; L0b; L02; L00; 
        Maa; Mab; Ma2; Ma0; 
        Mba; Mbb; Mb2; Mb0; 
        M2a; M2b; M22; M20; 
        M0a; M0b; M02; M00;
        La; Lb; L2; L0; Ma; Mb; M2a_swell;M2b_swell; M0a_swell; M0b_swell;
        P_mat;
    end
    methods
        function obj = Triangle(i, j, k, pi, pj, pk, h)  
            obj.Vertices = [ i, j, k ];  % vetex的复数，顶点指数1 2 3
            obj.Area = det( [ pj-pi, pk-pi ] )/2; %det( [ pj-pi, pk-pi ] ) 表示平行四边形的面积
            obj.Thickness = h;
            obj.vector_a = ( 1/(2*obj.Area))*[ pj(2)-pk(2); pk(2)-pi(2); pi(2)-pj(2) ]; % 2.3.11   
            obj.vector_b = (-1/(2*obj.Area))*[ pj(1)-pk(1); pk(1)-pi(1); pi(1)-pj(1) ];
            kexi_i = pi(1);  obj.x_i = kexi_i;
            kexi_j = pj(1);  obj.x_j = kexi_j;
            kexi_k = pk(1);  obj.x_k = kexi_k;
    
            yita_i = pi(2);  obj.y_i = yita_i;
            yita_j = pj(2);  obj.y_j = yita_j;
            yita_k = pk(2);  obj.y_k = yita_k;
            %LZY   9.12上方
            obj.vector_a2 = 1/(4*obj.Area^2).*[(yita_j - yita_k)^2;(yita_k - yita_i)^2;(yita_i - yita_j)^2];
            obj.vector_b2 = 1/(4*obj.Area^2).*[(kexi_j - kexi_k)^2;(kexi_k - kexi_i)^2;(kexi_i - kexi_j)^2];
            obj.vector_c2 = -1/(4*obj.Area^2).*[(yita_j - yita_k)*(kexi_j - kexi_k);(yita_k - yita_i)*(kexi_k - kexi_i);(yita_i - yita_j)*(kexi_i - kexi_j)];
            obj.vector_a0 = 2/(4*obj.Area^2).*[(yita_j - yita_k)*(yita_k - yita_i);(yita_j - yita_k)*(yita_i - yita_j);(yita_k - yita_i)*(yita_i - yita_j)];
            obj.vector_b0 = 2/(4*obj.Area^2).*[(kexi_j - kexi_k)*(kexi_k - kexi_i);(kexi_j - kexi_k)*(kexi_i - kexi_j);(kexi_k - kexi_i)*(kexi_i - kexi_j)];
            obj.vector_c0 = -1/(4*obj.Area^2).*[(yita_j - yita_k)*(kexi_k - kexi_i) + (kexi_j - kexi_k)*(yita_k - yita_i);...
                                                (yita_j - yita_k)*(kexi_i - kexi_j) + (kexi_j - kexi_k)*(yita_i - yita_j);...
                                                (yita_k - yita_i)*(kexi_i - kexi_j) + (kexi_k - kexi_i)*(yita_i - yita_j)];
            obj.vector_T2 = obj.vector_a2 + obj.vector_b2;
            obj.vector_T0 = obj.vector_a0 + obj.vector_b0;
            
            a = obj.vector_a;  a2 =  obj.vector_a2; a0 = obj.vector_a0;
            b = obj.vector_b;  b2 =  obj.vector_b2; b0 = obj.vector_b0;
            c2 =  obj.vector_c2; c0 = obj.vector_c0;
            T2 =  obj.vector_T2; T0 = obj.vector_T0;
          
            %LZY   2023/10/26
            obj.Laa = h*obj.Area.*a*transpose(a);  
            obj.Lab = h*obj.Area.*a*transpose(b);  
            obj.La2 = 0.5*h*obj.Area.*a*transpose(T2);  
            obj.La0 = 0.5*h*obj.Area.*a*transpose(T0);  

            obj.Lba = transpose(obj.Lab);
            obj.Lbb = h*obj.Area.*b*transpose(b);
            obj.Lb2 = 0.5*h*obj.Area.*b*transpose(T2);
            obj.Lb0 = 0.5*h*obj.Area.*b*transpose(T0);

            obj.L2a = transpose(obj.La2);
            obj.L2b = transpose(obj.Lb2);
            obj.L22 = 0.25*h*obj.Area.*T2*transpose(T2);
            obj.L20 = 0.25*h*obj.Area.*T2*transpose(T0);

            obj.L0a = transpose(obj.La0);
            obj.L0b = transpose(obj.Lb0);
            obj.L02 = transpose(obj.L20);
            obj.L00 = 0.25*h*obj.Area.*T0*transpose(T0);

            obj.La = h*obj.Area.*a';
            obj.Lb = h*obj.Area.*b';
            obj.L2 = 0.5*h*obj.Area.*T2';
            obj.L0 = 0.5*h*obj.Area.*T0';

            obj.Maa = 2*obj.Laa + obj.Lbb;
            obj.Mab = obj.Lba;
            obj.Ma2 = h*obj.Area.*(a*transpose(a2) + b*transpose(c2));
            obj.Ma0 = h*obj.Area.*(a*transpose(a0) + b*transpose(c0));

            obj.Mba = transpose(obj.Mab);
            obj.Mbb = 2*obj.Lbb + obj.Laa;
            obj.Mb2 = h*obj.Area.*(b*transpose(b2) + a*transpose(c2));
            obj.Mb0 = h*obj.Area.*(b*transpose(b0) + a*transpose(c0));

            obj.M2a = transpose(obj.Ma2);
            obj.M2b = transpose(obj.Mb2);
            obj.M22 = h*obj.Area.*(0.5.*(a2*transpose(a2) + b2*transpose(b2)) + c2*transpose(c2));
            obj.M20 = h*obj.Area.*(0.5.*(a2*transpose(a0) + b2*transpose(b0)) + c2*transpose(c0));

            obj.M0a = transpose(obj.Ma0);
            obj.M0b = transpose(obj.Mb0);
            obj.M02 = transpose(obj.M20);
            obj.M00 = h*obj.Area.*(0.5.*(a0*transpose(a0)+ b0*transpose(b0)) + c0*transpose(c0));

            obj.Ma = 2*obj.La;
            obj.Mb = 2*obj.Lb;

            obj.M2a_swell = h*obj.Area.*a2';
            obj.M2b_swell = h*obj.Area.*b2';
            obj.M0a_swell = h*obj.Area.*a0';
            obj.M0b_swell = h*obj.Area.*b0';
            
            vol = obj.Area * obj.Thickness;
            luu = vol*a*a'; luv = vol*a*b';
            lvu = vol*b*a'; lvv = vol*b*b';
            muu = 2*luu + lvv; muv = lvu;
            mvu = luv; mvv = 2*lvv + luu;
            
            l = [ luu, luv; lvu, lvv ];
            m = [ muu, muv; mvu, mvv ];
    
            obj.Partial_J_lambda = l([1,4,2,5,3,6], [1,4,2,5,3,6]);
            obj.Partial_J_mu     = m([1,4,2,5,3,6], [1,4,2,5,3,6]); %2.3.14和2.3.15
            obj.P_mat =    [ 1 0 0 0 0 0;...
                             0 0 0 1 0 0;...
                             0 1 0 0 0 0;...
                             0 0 0 0 1 0;...
                             0 0 1 0 0 0;...
                             0 0 0 0 0 1];
        end

                % 2025 1 28
        function obj = update_mass(obj, pi, pj, pk)  
            obj.Area = det( [ pj-pi, pk-pi ] )/2; %det( [ pj-pi, pk-pi ] ) 表示平行四边形的面积
        end
        
        
        function obj = mechanical_parameters(obj, rho, l, m)
            obj.Density = rho;
            obj.lambda = l; obj.mu = m;
        end
        
        function obj = partial_derivaties(obj, ui, uj, uk)
            gamma_u =  [ ui(1); uj(1); uk(1) ];
            gamma_v =  [ ui(2); uj(2); uk(2) ];
            obj.u_x = obj.vector_a' * gamma_u;
            obj.u_y = obj.vector_b' * gamma_u;
            obj.v_x = obj.vector_a' * gamma_v;
            obj.v_y = obj.vector_b' * gamma_v;
        end
        %2.3.12
        function obj = calculate_Cauchy_strain(obj, ui, uj, uk)
            obj = obj.partial_derivaties (ui, uj, uk);
            obj.Cauchy_strain = [ obj.u_x; obj.v_y;  obj.u_y + obj.v_x ];
        end
        
        function energy = partial_strain_potential_energy(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3));
            obj = obj.calculate_Cauchy_strain(ui, uj, uk);
            energy = (1/2) * (obj.Area * obj.Thickness) * ...
                ( obj.lambda * ( obj.Cauchy_strain(1) + obj.Cauchy_strain(2) )^2 + ...
                  obj.mu * ( 2*obj.Cauchy_strain(1)^2 + 2*obj.Cauchy_strain(2)^2 + obj.Cauchy_strain(3)^2 ) );
        end

        function obj = calculate_Green_strain(obj, ui, uj, uk)
            obj = obj.partial_derivaties (ui, uj, uk);
            obj.Green_strain = [ ...
                obj.u_x + (1/2)*(obj.u_x^2 + obj.v_x^2); ...
                obj.v_y + (1/2)*(obj.u_y^2 + obj.v_y^2); ...
                obj.u_y + obj.v_x + (obj.u_x*obj.u_y + obj.v_x*obj.v_y) ];
        end
        
        function energy = partial_strain_potential_energy_Green_strain(obj, disps)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3));
            obj = obj.calculate_Green_strain(ui, uj, uk);
            energy = (1/2) * (obj.Area * obj.Thickness) * ...
                ( obj.lambda * ( obj.Green_strain(1) + obj.Green_strain(2) )^2 + ...
                  obj.mu * ( 2*obj.Green_strain(1)^2 + 2*obj.Green_strain(2)^2 + obj.Green_strain(3)^2 ) );
        end

        function energy = partial_strain_potential_energy_4D(obj,  disps, t)
            mea = obj.strain_mea(t);
            stain_4D = zeros(3,1);
            stain_4D(1) = mea(1)*mea(1) + 2*mea(1);
            stain_4D(2) = mea(2)*mea(2) + 2*mea(2);
            stain_4D(3) = 0;
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3));
            obj = obj.calculate_Green_strain(ui, uj, uk);

            stain_diff = obj.Green_strain - stain_4D;
            energy = (1/2) * (obj.Area * obj.Thickness) * ...
                ( obj.lambda * ( stain_diff(1) + stain_diff(2) )^2 + ...
                  obj.mu * ( 2*stain_diff(1)^2 + 2*stain_diff(2)^2 + stain_diff(3)^2 ) );
        end
        
        function energy = partial_gravitational_potential_energy(obj, disps, grav)
            vs = obj.Vertices;
            ui = disps(:,vs(1)); uj = disps(:,vs(2)); uk = disps(:,vs(3));
            energy = - (obj.Density*obj.Area*obj.Thickness) * grav' *(ui+uj+uk)/3;
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
            mass = obj.Density * obj.Area * obj.Thickness;
            %I = eye(2);
            %obj.Partial_Inertia_Matrix = (mass/12) * ...
            %    [2*I, I, I; I, 2*I, I; I, I, 2*I];
            obj.Partial_Inertia_Matrix = (mass/3) * eye(6); %lumped mass
        end
        
        function [obj, M_p] = partial_inertia_matrix(obj)
            if isempty( obj.Partial_Inertia_Matrix )
                obj = obj.calculate_partial_inertia_matrix;
            end
            M_p = obj.Partial_Inertia_Matrix;
        end
        
        function obj = calculate_partial_gravitational_vector(obj, g)
            mass = obj.Density * obj.Area * obj.Thickness;
            obj.Partial_Gravitational_Vector = mass/3 * [ g; g; g ]; 
        end
        
        function [obj, grav_p] = partial_gravitational_vector(obj, g)
            if isempty( obj.Partial_Gravitational_Vector )
                obj = obj.calculate_partial_gravitational_vector(g);
            end
            grav_p = obj.Partial_Gravitational_Vector;
        end
        
        function [fi, fj, fk] = nodal_forces_Cauchy_strain(obj, ui, uj, uk)
            a = obj.vector_a;
            b = obj.vector_b;
            obj = obj.calculate_Cauchy_strain (ui, uj, uk);
            mat = obj.lambda*[1,1,0; 1,1,0; 0,0,0] + obj.mu*[2,0,0; 0,2,0; 0,0,1];
            Up_e = mat*obj.Cauchy_strain*obj.Area*obj.Thickness;
            Up_gammau = [ a, zeros(3,1), b ]*Up_e;
            Up_gammav = [ zeros(3,1), b, a ]*Up_e;
            fi = - [ Up_gammau(1); Up_gammav(1) ];
            fj = - [ Up_gammau(2); Up_gammav(2) ];
            fk = - [ Up_gammau(3); Up_gammav(3) ];
        end
        
        function [fi, fj, fk] = nodal_forces_Green_strain(obj, ui, uj, uk)
            a = obj.vector_a;
            b = obj.vector_b;
            obj = obj.calculate_Green_strain (ui, uj, uk);
            mat = obj.lambda*[1,1,0; 1,1,0; 0,0,0] + obj.mu*[2,0,0; 0,2,0; 0,0,1];
            Up_E = mat*obj.Green_strain*obj.Area*obj.Thickness;
            ux = obj.u_x; uy = obj.u_y; vx = obj.v_x; vy = obj.v_y;
            Up_gammau = [ (1+ux)*a, uy*b, (1+ux)*b+uy*a ]*Up_E;
            Up_gammav = [ vx*a, (1+vy)*b, (1+vy)*a+vx*b ]*Up_E;
            fi = - [ Up_gammau(1); Up_gammav(1) ];
            fj = - [ Up_gammau(2); Up_gammav(2) ];
            fk = - [ Up_gammau(3); Up_gammav(3) ];
        end
         
       %LZY   2024/3/11
       function [fi, fj, fk] = nodal_total_forces_Green_strain(obj, dispi, dispj, dispk, veli, velj, velk)
            ui = dispi(1);uj = dispj(1);uk = dispk(1);
            vi = dispi(2);vj = dispj(2);vk = dispk(2);
            ui_dot = veli(1);uj_dot = velj(1); uk_dot = velk(1);
            vi_dot = veli(2);vj_dot = velj(2); vk_dot = velk(2);
            
            gamma_u = [ui;uj;uk];
            gamma_v = [vi;vj;vk];
            gamma_u2 =[ui*ui;uj*uj;uk*uk];
            gamma_v2 = [vi*vi;vj*vj;vk*vk];
            gamma_uu = [ui*uj;ui*uk;uj*uk];
            gamma_vv = [vi*vj;vi*vk;vj*vk];
            gamma_2 = gamma_u2 + gamma_v2;
            gamma_0 = gamma_uu + gamma_vv;

            gamma_u_dot = [ui_dot;uj_dot;uk_dot];
            gamma_v_dot = [vi_dot;vj_dot;vk_dot];
            gamma_u2_dot = [2*ui*ui_dot;2*uj*uj_dot;2*uk*uk_dot];
            gamma_v2_dot = [2*vi*vi_dot;2*vj*vj_dot;2*vk*vk_dot];
            gamma_2_dot = gamma_u2_dot + gamma_v2_dot;
            gamma_uu_dot =[ui_dot*uj + ui*uj_dot;ui_dot*uk + ui*uk_dot;uj_dot*uk + uj*uk_dot];
            gamma_vv_dot =[vi_dot*vj + vi*vj_dot;vi_dot*vk + vi*vk_dot;vj_dot*vk + vj*vk_dot];
            gamma_0_dot = gamma_uu_dot + gamma_vv_dot;

            gamma2_gamma_u = diag([2*ui;2*uj;2*uk]);
            gamma0_gamma_u =[uj uk 0;ui 0 uk;0 ui uj];
            gamma2_gamma_v = diag([2*vi;2*vj;2*vk]);
            gamma0_gamma_v = [vj vk 0;vi 0 vk;0 vi vj];

            gamma2_gamma_u_dot = diag([2*ui_dot;2*uj_dot;2*uk_dot]);
            gamma0_gamma_u_dot = [uj_dot uk_dot 0;ui_dot 0 uk_dot;0 ui_dot uj_dot];
            gamma2_gamma_v_dot = diag([2*vi_dot;2*vj_dot;2*vk_dot]);
            gamma0_gamma_v_dot = [vj_dot vk_dot 0;vi_dot 0 vk_dot;0 vi_dot vj_dot];

       
            g_u_lam = obj.Laa*gamma_u + obj.Lab*gamma_v + obj.La2*gamma_2 + obj.La0*gamma_0;
            g_v_lam = obj.Lba*gamma_u + obj.Lbb*gamma_v + obj.Lb2*gamma_2 + obj.Lb0*gamma_0;
            g_2_lam = obj.L2a*gamma_u + obj.L2b*gamma_v + obj.L22*gamma_2 + obj.L20*gamma_0;
            g_0_lam = obj.L0a*gamma_u + obj.L0b*gamma_v + obj.L02*gamma_2 + obj.L00*gamma_0;
            g_u_mu =  obj.Maa*gamma_u + obj.Mab*gamma_v + obj.Ma2*gamma_2 + obj.Ma0*gamma_0;
            g_v_mu =  obj.Mba*gamma_u + obj.Mbb*gamma_v + obj.Mb2*gamma_2 + obj.Mb0*gamma_0;
            g_2_mu =  obj.M2a*gamma_u + obj.M2b*gamma_v + obj.M22*gamma_2 + obj.M20*gamma_0;
            g_0_mu =  obj.M0a*gamma_u + obj.M0b*gamma_v + obj.M02*gamma_2 + obj.M00*gamma_0;

            g_u_lam_dot = obj.Laa*gamma_u_dot + obj.Lab*gamma_v_dot + obj.La2*gamma_2_dot + obj.La0*gamma_0_dot;
            g_v_lam_dot = obj.Lba*gamma_u_dot + obj.Lbb*gamma_v_dot + obj.Lb2*gamma_2_dot + obj.Lb0*gamma_0_dot;
            g_2_lam_dot = obj.L2a*gamma_u_dot + obj.L2b*gamma_v_dot + obj.L22*gamma_2_dot + obj.L20*gamma_0_dot;
            g_0_lam_dot = obj.L0a*gamma_u_dot + obj.L0b*gamma_v_dot + obj.L02*gamma_2_dot + obj.L00*gamma_0_dot;

            g_u_mu_dot = obj.Maa*gamma_u_dot + obj.Mab*gamma_v_dot + obj.Ma2*gamma_2_dot + obj.Ma0*gamma_0_dot;
            g_v_mu_dot = obj.Mba*gamma_u_dot + obj.Mbb*gamma_v_dot + obj.Mb2*gamma_2_dot + obj.Mb0*gamma_0_dot;
            g_2_mu_dot = obj.M2a*gamma_u_dot + obj.M2b*gamma_v_dot + obj.M22*gamma_2_dot + obj.M20*gamma_0_dot;
            g_0_mu_dot = obj.M0a*gamma_u_dot + obj.M0b*gamma_v_dot + obj.M02*gamma_2_dot + obj.M00*gamma_0_dot;

            G_lam_gamma_u = g_u_lam + gamma2_gamma_u*g_2_lam + gamma0_gamma_u*g_0_lam;
            G_lam_gamma_v = g_v_lam + gamma2_gamma_v*g_2_lam + gamma0_gamma_v*g_0_lam;

            G_mu_gamma_u = g_u_mu + gamma2_gamma_u*g_2_mu + gamma0_gamma_u*g_0_mu;
            G_mu_gamma_v = g_v_mu + gamma2_gamma_v*g_2_mu + gamma0_gamma_v*g_0_mu;
    
            G_lam_gamma_u_dot = g_u_lam_dot + gamma2_gamma_u_dot*g_2_lam + gamma2_gamma_u*g_2_lam_dot + gamma0_gamma_u_dot*g_0_lam + gamma0_gamma_u*g_0_lam_dot;
            G_lam_gamma_v_dot = g_v_lam_dot + gamma2_gamma_v_dot*g_2_lam + gamma2_gamma_v*g_2_lam_dot + gamma0_gamma_v_dot*g_0_lam + gamma0_gamma_v*g_0_lam_dot;

            G_mu_gamma_u_dot = g_u_mu_dot + gamma2_gamma_u_dot*g_2_mu + gamma2_gamma_u*g_2_mu_dot + gamma0_gamma_u_dot*g_0_mu + gamma0_gamma_u*g_0_mu_dot;
            G_mu_gamma_v_dot = g_v_mu_dot + gamma2_gamma_v_dot*g_2_mu + gamma2_gamma_v*g_2_mu_dot + gamma0_gamma_v_dot*g_0_mu + gamma0_gamma_v*g_0_mu_dot;
            
          
            F_ela = obj.P_mat*([obj.lambda.*G_lam_gamma_u + obj.mu.*G_mu_gamma_u;...
                         obj.lambda.*G_lam_gamma_v + obj.mu.*G_mu_gamma_v]);
            
            F_vis =obj.P_mat*([obj.lambda_vis.*G_lam_gamma_u_dot + obj.mu_vis.*G_mu_gamma_u_dot;...
                             obj.lambda_vis.*G_lam_gamma_v_dot + obj.mu_vis.*G_mu_gamma_v_dot]);
                         
            fi = - [  F_ela(1)+F_vis(1);  F_ela(2)+F_vis(2) ];
            fj = - [  F_ela(3)+F_vis(3);  F_ela(4)+F_vis(4) ];
            fk = - [  F_ela(5)+F_vis(5);  F_ela(6)+F_vis(6) ];
       end


       % LZY 2024 7 26  nodal_forces external stimuli  best 
        function [fi, fj, fk] = nodal_forces_external_stimuli(obj, dispi, dispj, dispk, t)
         
            ui = dispi(1); uj = dispj(1); uk = dispk(1);
            vi = dispi(2); vj = dispj(2); vk = dispk(2);
        
            mea = obj.strain_mea(t);

            alpha = mea + 0.5*mea.*mea;
           
            alpha_sum = sum(alpha);
            g_u_lam_swell = alpha_sum * obj.La';
            g_v_lam_swell = alpha_sum * obj.Lb';
            g_2_lam_swell = alpha_sum * obj.L2';
            g_0_lam_swell = alpha_sum * obj.L0';
        
            g_u_mu_swell =  alpha(1) * obj.Ma';
            g_v_mu_swell =  alpha(2) * obj.Mb';
            g_2_mu_swell =  alpha(1) * obj.M2a_swell' + alpha(2) * obj.M2b_swell';
            g_0_mu_swell =  alpha(1) * obj.M0a_swell' + alpha(2) * obj.M0b_swell';
        
            gamma2_gamma_u = diag([2*ui, 2*uj, 2*uk]);
            gamma2_gamma_v = diag([2*vi, 2*vj, 2*vk]);
            gamma0_gamma_u = [uj, uk, 0; ui, 0, uk; 0, ui, uj];
            gamma0_gamma_v = [vj, vk, 0; vi, 0, vk; 0, vi, vj];
        
            G_lam_gamma_u_swell = g_u_lam_swell + gamma2_gamma_u * g_2_lam_swell + gamma0_gamma_u * g_0_lam_swell;
            G_lam_gamma_v_swell = g_v_lam_swell + gamma2_gamma_v * g_2_lam_swell + gamma0_gamma_v * g_0_lam_swell;
            G_mu_gamma_u_swell = g_u_mu_swell + gamma2_gamma_u * g_2_mu_swell + gamma0_gamma_u * g_0_mu_swell;
            G_mu_gamma_v_swell = g_v_mu_swell + gamma2_gamma_v * g_2_mu_swell + gamma0_gamma_v * g_0_mu_swell;
        
            F_swell = obj.P_mat * [obj.lambda * G_lam_gamma_u_swell + obj.mu * G_mu_gamma_u_swell; ...
                                   obj.lambda * G_lam_gamma_v_swell + obj.mu * G_mu_gamma_v_swell];
        
            fi = F_swell(1:2);
            fj = F_swell(3:4);
            fk = F_swell(5:6);
        end

    end
end
