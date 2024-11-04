function elastic = meshObject(app) 
     elastic = [];
     try
         value = app.TypeDropDown.Value;   
         if strcmp(value, '2D')
             elastic= mesh2DObject(app);
         elseif strcmp(value, '3D')
             elastic= mesh3DObject(app); 
         end  
         initPlot(elastic, app);
     catch ME
          msgbox(['Error occurred: ', ME.message], 'Error', 'error','modal');
     end
end
 
function  elastic = mesh2DObject(app)
    scale = app.SizeScale;

    shape = app.ShapeDropDown2D.Value;
    switch shape
        case "Rectangular"
            x_length = app.x_2D; y_length = app.y_2D; thickness = app.z_2D;   
            x_length =x_length*scale; y_length = y_length*scale; thickness = thickness*scale; 
            m = app.numNode_x_2D; n = app.numNode_y_2D; 
            [points, triangles] = two_dim_fea.rectangular_object(m, n, x_length, y_length); 
        case "Doughnut"
             router = app.x_2D; rinner = app.y_2D; thickness = app.z_2D;   
            if router > rinner
                router =router*scale; rinner = rinner*scale; thickness = thickness*scale; 
                m = app.numNode_x_2D; n = app.numNode_y_2D; 
                [points, triangles] = two_dim_fea.ring_object(m, n, router, rinner);    
            else
                 error("Outer R should be bigger than inner R.")
            end
        case "Half Doughnut"
             router = app.x_2D; rinner = app.y_2D; thickness = app.z_2D;  thetastart = app.theta_2D/180*pi; 
            if router > rinner
                router =router*scale; rinner = rinner*scale; thickness = thickness*scale; 
                m = app.numNode_x_2D; n = app.numNode_y_2D; 
                [points, triangles] = two_dim_fea.half_ring_object( m, n, router, rinner, thetastart);    
            else
                 error("Outer R should be bigger than inner R.")
            end
    end

    npoints = size(points,2);
    ntriangles = size(triangles,1);
    elastic = two_dim_fea.Body(npoints, points, ntriangles, triangles, thickness);
end


function  elastic = mesh3DObject(app)
    scale = app.SizeScale;
    x_length = app.x_3D; y_length = app.y_3D; z_length = app.z_3D;   
    x_length =x_length*scale; y_length = y_length*scale; z_length = z_length*scale; 
    l = app.numNode_x_3D; m = app.numNode_y_3D; n = app.numNode_z_3D;       
    [ points, tetrahedron ] = three_dim_fea.cuboidal_object( l, m, n, x_length, y_length, z_length );
    npoints = size(points,2);
    ntetrahedron = size(tetrahedron,1);
    elastic = three_dim_fea.Body(npoints, points, ntetrahedron, tetrahedron);
end

function [] = initPlot(elastic, app)
    cla(app.UIAxes, 'reset');   % clear app.UIAxes
    hold(app.UIAxes, 'on');
    npoints = elastic.numNodalPoints;
    init_coor = [];
    for i = 1:npoints
        init_coor = [init_coor, elastic.NodalPoints(i).Coordinates];
    end
    coor_row = size(init_coor,1);
    disps_init = zeros([coor_row,elastic.numNodalPoints]);
    elastic.draw_init_app(app.UIAxes, disps_init); 
    sizeScale = app.SizeScale;

    overall_max = max(init_coor,[],2);  
    overall_min = min(init_coor,[],2);
    diff_max_min = abs(overall_max - overall_min);
    
    axis(app.UIAxes, 'equal');     % must be set before the XYZ ticklabel changed 
    xlimRangeMin = overall_min(1) - 0.1*diff_max_min(1);
    xlimRangeMax = overall_max(1) + 0.1*diff_max_min(1);

    ylimRangeMin = overall_min(2) - 0.1*diff_max_min(2);
    ylimRangeMax = overall_max(2) + 0.1*diff_max_min(2);  

    xlim(app.UIAxes, [ xlimRangeMin, xlimRangeMax]);
    ylim(app.UIAxes, [ ylimRangeMin, ylimRangeMax]);
    view(app.UIAxes,2);
    if coor_row ==3
         zlimRangeMin = overall_min(3) - 0.1*diff_max_min(3);
         zlimRangeMax = overall_max(3) + 0.1*diff_max_min(3);  
         zlim(app.UIAxes, [ zlimRangeMin, zlimRangeMax]);
         view(app.UIAxes,[30,10])
     end
   
    % set tickLabel
    % xticks = get(app.UIAxes, 'XTick');
    % yticks = get(app.UIAxes, 'YTick');
    % zticks = get(app.UIAxes, 'ZTick');
    % 
    % %epsilon = 1e-10;
    % xticklabels = arrayfun(@(x) num2str(x*(1/sizeScale)), xticks, 'UniformOutput', false);
    % yticklabels = arrayfun(@(y) num2str(y*(1/sizeScale)), yticks, 'UniformOutput', false);
    % zticklabels = arrayfun(@(z) num2str(z*(1/sizeScale)), zticks, 'UniformOutput', false);
    % 
    % % xticklabels(abs(xticks) < epsilon) = {'0'};
    % % yticklabels(abs(yticks) < epsilon) = {'0'};
    % % zticklabels(abs(zticks) < epsilon) = {'0'};
    % 
    % set(app.UIAxes, 'XTickLabel', xticklabels);
    % set(app.UIAxes, 'YTickLabel', yticklabels);
    % set(app.UIAxes, 'ZTickLabel', zticklabels);


    xticks = get(app.UIAxes, 'XTick');
    yticks = get(app.UIAxes, 'YTick');
    zticks = get(app.UIAxes, 'ZTick');
    
    %epsilon = 1e-10;
    
    xticklabelsArray = arrayfun(@(x) num2str(x *(1/sizeScale)), xticks, 'UniformOutput', false);
    yticklabelsArray = arrayfun(@(y) num2str(y *(1/sizeScale)), yticks, 'UniformOutput', false);
    zticklabelsArray = arrayfun(@(z) num2str(z *(1/sizeScale)), zticks, 'UniformOutput', false);
    
    % 
    % xticklabelsArray(abs(xticks) < epsilon) = {'0'};
    % yticklabelsArray(abs(yticks) < epsilon) = {'0'};
    % zticklabelsArray(abs(zticks) < epsilon) = {'0'};
    % 

    xticklabels(app.UIAxes, xticklabelsArray);
    yticklabels(app.UIAxes, yticklabelsArray);
    zticklabels(app.UIAxes, zticklabelsArray);

    app.UIAxes.LineWidth = 1.5;
    xlabel(app.UIAxes,"$x$ [m]",'Interpreter','latex','FontName','Times New Roman','FontSize',20);
    ylabel(app.UIAxes,"$y$ [m]",'Interpreter','latex','FontName','Times New Roman','FontSize',20);
    zlabel(app.UIAxes,"$z$ [m]",'Interpreter','latex','FontName','Times New Roman','FontSize',20);


    hold(app.UIAxes, 'off'); % Turn off the hold state to avoid affecting subsequent plotting.     
end