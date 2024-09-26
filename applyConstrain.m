function [A,b0,b1,constraintMarkers] = applyConstrain(app)
   try
         value = app.TypeDropDown.Value;   
         if strcmp(value, '2D')
             [A,b0,b1,constraintMarkers] = applyConstrain2D(app);
         elseif strcmp(value, '3D')
             [A,b0,b1,constraintMarkers] = applyConstrain3D(app);
         end  
   catch ME
         msgbox(['Error occurred: ', ME.message], 'Error', 'error','modal');
   end
end

function [A,b0,b1,constraintMarkers] = applyConstrain2D(app)
    app_ax = app.UIAxes;
    x_con = app.xCon2D;
    y_con = app.yCon2D;
    xy_con = app.xyCon2D;
    elastic = app.object4D;
    hold(app_ax, 'on');
    cn = length(x_con) + length(y_con) + 2*length(xy_con);
    numMaker = length(x_con) + length(y_con) + length(xy_con);
    constraintMarkers = gobjects(1,numMaker);
    if cn == 0
        A_trans = 0;
        A = A_trans.';
        b0 = 0;
        b1 = 0;  
        constraintMarkers = gobjects(1); 
    else
        A_trans = zeros(cn,2*elastic.numNodalPoints);
        row = 1;   nMarker = 1;
        if ~isempty(x_con)
            for i = 1:length(x_con)
                A_trans(row, (x_con(i)*2) - 1) = 1;
                row = row + 1;
                % draw constrain in x direction
                coor_con = elastic.NodalPoints(x_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot(app_ax, coor_con(1), coor_con(2),'rx', 'MarkerSize', 10, 'LineWidth', 2);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end
        if ~isempty(y_con)
            for i = 1:length(y_con)
                A_trans(row, y_con(i)*2) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(y_con(i)).Coordinates;  
                constraintMarkers(nMarker)  = plot(app_ax, coor_con(1), coor_con(2),'gx', 'MarkerSize', 10, 'LineWidth', 2);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end
        
        % xy
        if ~isempty(xy_con)
            for i = 1:length(xy_con)
                %x
                A_trans(row, (xy_con(i)*2) - 1) = 1;
                row = row + 1;
                %y
                A_trans(row, xy_con(i)*2) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(xy_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot(app_ax, coor_con(1), coor_con(2),'bx', 'MarkerSize', 10, 'LineWidth', 2);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end
        A = A_trans.';
        b0 = zeros(cn,1); b1 = zeros(cn,1);
    end
    hold(app_ax, 'off');
end

function [A,b0,b1,constraintMarkers] = applyConstrain3D(app)
    app_ax = app.UIAxes;
    x_con = app.xCon3D;
    y_con = app.yCon3D;
    z_con = app.zCon3D;
    xy_con = app.xyCon3D;
    xz_con = app.xzCon3D;
    yz_con = app.yzCon3D;
    xyz_con = app.xyzCon3D;
    
    elastic = app.object4D;
    hold(app_ax, 'on');
    cn = length(x_con) + length(y_con) + length(z_con)+ 2*length(xy_con)+2*length(xz_con)+2*length(yz_con)+3*length(xyz_con);
    
    allpoints = [x_con,y_con,z_con,xy_con,xz_con,yz_con,xyz_con];
    numMaker = length(allpoints);
    constraintMarkers = gobjects(1,numMaker);
    if cn == 0
        A_trans = 0;
        A = A_trans.';
        b0 = 0;
        b1 = 0;  
        constraintMarkers = gobjects(1);  
    else
        A_trans = zeros(cn,3*elastic.numNodalPoints);
        row = 1;   nMarker = 1;
        if ~isempty(x_con)
            for i = 1:length(x_con)
                A_trans(row, x_con(i)*3 - 2) = 1;
                row = row + 1;
                % draw constrain in x direction
                coor_con = elastic.NodalPoints(x_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot3(app_ax, coor_con(1), coor_con(2), coor_con(3),...
                                            'x', 'MarkerSize', 10, 'LineWidth', 2, ...
                                            'Color', [0.255, 0.412, 0.882]);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end
        if ~isempty(y_con)
            for i = 1:length(y_con)
                A_trans(row, y_con(i)*3-1) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(y_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot3(app_ax, coor_con(1), coor_con(2), coor_con(3),...
                                            'x', 'MarkerSize', 10, 'LineWidth', 2, ...
                                            'Color', [0.133, 0.545, 0.133]);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end      
       % z
        if ~isempty(z_con)
            for i = 1:length(z_con)
                A_trans(row, z_con(i)*3) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(z_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot3(app_ax, coor_con(1), coor_con(2), coor_con(3),...
                                            'x', 'MarkerSize', 10, 'LineWidth', 2, ...
                                            'Color', [0.855, 0.647, 0.125]);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end
 
        % xy 
        if ~isempty(xy_con)
            for i = 1:length(xy_con)
                A_trans(row, xy_con(i)*3 - 2) = 1;
                row = row + 1;
                A_trans(row, xy_con(i)*3- 1) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(xy_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot3(app_ax, coor_con(1), coor_con(2), coor_con(3),...
                                            'x', 'MarkerSize', 10, 'LineWidth', 2, ...
                                            'Color', [0.576, 0.439, 0.859]);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end

         %xz
        if ~isempty(xz_con)
            for i = 1:length(xz_con)
                A_trans(row, xz_con(i)*3 - 2) = 1;
                row = row + 1;
                A_trans(row, xz_con(i)*3) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(xz_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot3(app_ax, coor_con(1), coor_con(2), coor_con(3),...
                                            'x', 'MarkerSize', 10, 'LineWidth', 2, ...
                                            'Color', [0.439, 0.502, 0.565]);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end

         %yz
        if ~isempty(yz_con)
            for i = 1:length(yz_con)
                A_trans(row, yz_con(i)*3 - 1) = 1;
                row = row + 1;
                A_trans(row, yz_con(i)*3) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(yz_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot3(app_ax, coor_con(1), coor_con(2), coor_con(3),...
                                            'x', 'MarkerSize', 10, 'LineWidth', 2, ...
                                            'Color', [1.0, 0.549, 0.0]);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end

          %xyz
        if ~isempty(xyz_con)
            for i = 1:length(xyz_con)
                A_trans(row, xyz_con(i)*3 - 2) = 1;
                row = row + 1;
                A_trans(row, xyz_con(i)*3 - 1) = 1;
                row = row + 1;
                A_trans(row, xyz_con(i)*3) = 1;
                row = row + 1;
                coor_con = elastic.NodalPoints(xyz_con(i)).Coordinates;  
                constraintMarkers(nMarker) = plot3(app_ax, coor_con(1), coor_con(2), coor_con(3),...
                                            'x', 'MarkerSize', 10, 'LineWidth', 2, ...
                                            'Color', [0.863, 0.078, 0.235]); %[0.576, 0.439, 0.859]);
                uistack(constraintMarkers(nMarker), 'top');
                nMarker = nMarker+1;
            end
        end
        A = A_trans.';
        b0 = zeros(cn,1); b1 = zeros(cn,1);
    end
    hold(app_ax, 'off');
end
