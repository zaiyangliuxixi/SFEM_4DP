function [] = setSubregion(app)
     value = app.TypeDropDown.Value;   
     if strcmp(value, '2D')
         setSubregion2D(app);
     elseif strcmp(value, '3D')
         setSubregion3D(app);
     end  
end
function [] = setSubregion2D(app)
    NumSubregions = length(app.SubregionProperties);
    % delete all subregion before apply
    app.object4D.SubRegions =[];
    app.object4D.numSubRegions = [];
    for p = 1:app.object4D.numTriangles
        app.object4D.Triangles(p).color =[];
    end
    epsilon = 1e-10; 

    try
        for i = 1:NumSubregions
            subtext = app.SubregionProperties(i).ElementIndices;
            parts = strsplit(subtext, ' ');
            % 初始化一个空数组来存放结果
            sub = [];
            % 遍历每个部分
            for j = 1:length(parts)
                % 对每个部分使用 eval 来解析范围表达式
                range = eval(parts{j});
                % 将解析后的结果添加到结果数组中
                sub = [sub, range];
            end         
             %check Young's modulus
            Young =app.SubregionProperties(i).YoungsModulus;
            if abs(app.SubregionProperties(i).YoungsModulus) < epsilon
                errordlg([' Young''s Modulus of Subregion ' num2str(i) ' is 0.'], 'Error','modal');
                return;
            end
            %check density
             density = app.SubregionProperties(i).Density;
            if abs(app.SubregionProperties(i).Density) < epsilon
                errordlg([' Density of Subregion ' num2str(i) ' is 0.'], 'Error','modal');
                return;
            end
            c = app.SubregionProperties(i).ViscoCoefficient; 
            nu = app.SubregionProperties(i).PoissonsRatio; 
            [ lambda_ela, mu_ela ] = Lame_constants( Young, nu );
            [ lambda_vis, mu_vis ] = Lame_constants( c, nu );          
    
            strain_mea_x = parseMeasuredStrain(app.SubregionProperties(i).MeasuredStrainX);
            strain_mea_y = parseMeasuredStrain(app.SubregionProperties(i).MeasuredStrainY);
            strain_mea = @(t) [strain_mea_x(t); strain_mea_y(t)];
            
            app.object4D= app.object4D.define_subregion(sub);
    
            app.object4D = app.object4D.subregion_mechanical_parameters(density, lambda_ela, mu_ela);
            app.object4D = app.object4D.subregion_viscous_parameters( lambda_vis, mu_vis);
            app.object4D = app.object4D.subregion_strain_mea(strain_mea);        
            
            color = app.SubregionProperties(i).Color;
            app.object4D= app.object4D.subregion_color(color./255 );
            app.issetSubregion = true;
        end
    catch ME   
        app.issetSubregion =false;
        errordlg(ME.message, 'Error', 'modal');
        return;
    end
    app_ax = app.UIAxes;
    hold(app_ax, 'on');
    disps_init = zeros([2,app.object4D.numNodalPoints]);
    app.object4D.draw_subregion_app(app_ax, disps_init);  
    hold(app_ax, 'off'); 
end


function [] = setSubregion3D(app)
    NumSubregions = length(app.SubregionProperties);
    app.object4D.SubRegions =[];
    app.object4D.numSubRegions = [];
    for p = 1:app.object4D.numTetrahedrons
        app.object4D.Tetrahedrons(p).color =[];
    end
    epsilon = 1e-10;  % 设定一个很小的阈值
    try
        for i = 1:NumSubregions
            subtext = app.SubregionProperties(i).ElementIndices;
            parts = strsplit(subtext, ' ');
            % 初始化一个空数组来存放结果
            sub = [];
            % 遍历每个部分
            for j = 1:length(parts)
                % 对每个部分使用 eval 来解析范围表达式
                range = eval(parts{j});
                % 将解析后的结果添加到结果数组中
                sub = [sub, range];
            end        
            Young = app.SubregionProperties(i).YoungsModulus;
            %check Young's modulus
            if abs(app.SubregionProperties(i).YoungsModulus) < epsilon
                errordlg([' Young''s Modulus of Subregion ' num2str(i) ' is 0.'], 'Error','modal');
                return;
            end
            
            %check density
             density = app.SubregionProperties(i).Density;
            if abs(app.SubregionProperties(i).Density) < epsilon
                errordlg([' Density of Subregion ' num2str(i) ' is 0.'], 'Error','modal');
                return;
            end
            
            c = app.SubregionProperties(i).ViscoCoefficient; 
            nu = app.SubregionProperties(i).PoissonsRatio; 
            [ lambda_ela, mu_ela ] = Lame_constants( Young, nu );
            [ lambda_vis, mu_vis ] = Lame_constants( c, nu );          
    
    
            strain_mea_x = parseMeasuredStrain(app.SubregionProperties(i).MeasuredStrainX);
            strain_mea_y = parseMeasuredStrain(app.SubregionProperties(i).MeasuredStrainY);
            strain_mea_z = parseMeasuredStrain(app.SubregionProperties(i).MeasuredStrainZ);
            strain_mea = @(t) [strain_mea_x(t); strain_mea_y(t);strain_mea_z(t)];
            
            app.object4D= app.object4D.define_subregion(sub);
    
            app.object4D = app.object4D.subregion_mechanical_parameters(density, lambda_ela, mu_ela);
            app.object4D = app.object4D.subregion_viscous_parameters( lambda_vis, mu_vis);
            app.object4D = app.object4D.subregion_strain_mea(strain_mea);        
            
            color = app.SubregionProperties(i).Color;
            app.object4D= app.object4D.subregion_color(color./255 );
            app.issetSubregion = true;
        end
    catch ME   
        errordlg(ME.message, 'Error', 'modal');
        app.issetSubregion =false;
        return;
    end

    app_ax = app.UIAxes;
    hold(app_ax, 'on');
    disps_init = zeros([3,app.object4D.numNodalPoints]);
    app.object4D.draw_subregion_app(app_ax, disps_init);  

    drawnow;
    hold(app_ax, 'off'); 
end


function strain_mea = parseMeasuredStrain(inputStr)
    % default
    defaultValues = {'$u^{msr}_x(t)$', '$v^{msr}_y(t)$', '$w^{msr}_z(t)$'};
    try
        if any(strcmp(inputStr, defaultValues))
            error('Please enter a valid MATLAB expression for measured strain.');
        end
        strain_mea = eval(['@(t)', inputStr]);
        %testValue = strain_mea(1); 
    catch ME
        errordlg(ME.message, 'Error', 'modal');
        return;  
    end
end


function [ lambda, mu ] = Lame_constants ( E, nu )
    % calculating Lame's constants
    lambda = nu*E/(1+nu)/(1-2*nu);
    mu = E/2/(1+nu);
end
