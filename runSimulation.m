function [time, q] = runSimulation(app)  
    time =[];
    q = [];
    sizeScale = app.SizeScale;  
    try
        %initial coordinate of point 1
        init_coor_p1 = app.object4D.NodalPoints(1).Coordinates;
        coor_row = size(init_coor_p1,1);
        timeUnit = app.TimeUnitDropDown.Value;
        if timeUnit == "s"
             g_scale = 1;
        elseif timeUnit == "min"
             g_scale = 60*60;
        elseif timeUnit == "h"
             g_scale = 60*60*60*60;
        end

        if coor_row ==2
            g = [0;-9.81]*g_scale*sizeScale;
        elseif coor_row ==3
            g = [0;0;-9.81]*g_scale*sizeScale;
        end
        app.object4D = app.object4D.calculate_gravitational_vector(g);

        app.object4D = app.object4D.calculate_inertia_matrix;

        alpha = 2000;    
        interval = [0,app.SimulationTime];    %simulation time
    
        A = app.A;
        b0 = app.b0;
        b1 =app.b1;
        M = app.object4D.Inertia_Matrix;
        
        if  nnz(A) == 0
            coef = M;
        else
            coef = [ M, -A; -A', zeros(size(A,2),size(A,2))];
        end
        coef_inv = inv(coef);

    

        qinit = zeros(2*coor_row*app.object4D.numNodalPoints,1);
        f_hydrogel_swelling = @(t,q) dynamic_equation(t,q, app, A,b0,b1, alpha,coef_inv,coor_row);
    
        options = odeset('RelTol',1e-03,'AbsTol',1e-06,'Events', @(time, q) cancelEventFcn(time, q, app));                        
        [time, q] = ode45(f_hydrogel_swelling, interval, qinit,options);
    catch ME
        msgbox(['Error occurred: ', ME.message], 'Error', 'error','modal');
    end
end


function [value,isterminal,direction] = cancelEventFcn(~,~, app)   
    % 如果用户按下了取消按钮，事件值设为 0，触发停止
    if app.SimulationProgressDlg.CancelRequested
        value = 0; % 事件发生时 value = 0
    else
        value = 1; % 正常情况下 value 不为 0，不触发事件
    end
    isterminal = 1; % 事件发生时终止解算
    direction = 0;  % 无论递增或递减均触发事件
end


function dotq = dynamic_equation(t,q, app, A,b0,b1, alpha,coef_inv,coor_row)
   % 计算进度百分比
    app.newSimProgress = t / app.SimulationTime;    
    % 更新进度条
    if  app.newSimProgress > max(app.oldSimProgress)
        app.SimulationProgressDlg.Value = app.newSimProgress;
        app.SimulationProgressDlg.Message = sprintf('Simulation Progress: %.0f%%', app.newSimProgress * 100);
    end
    app.oldSimProgress = [app.oldSimProgress,app.newSimProgress];


    npoints = app.object4D.numNodalPoints;
    
    un = q(1:coor_row*npoints);
    vn = q(coor_row*npoints+1 : 2*coor_row*npoints); 
    dotun = vn;
    forces_visela = app.object4D.nodal_total_forces_Green_strain(reshape(un, [coor_row,npoints]),reshape(vn, [coor_row,npoints]));   
    forces_4D = app.object4D.nodal_forces_external_stimuli(reshape(un, [coor_row,npoints]),t);
    forces_gravity = app.object4D.Gravitational_Vector;

    if app.isGravity
          total_forces = forces_4D+ forces_visela+forces_gravity;
    else 
          total_forces = forces_4D+ forces_visela;
    end
    if  nnz(A) == 0
        vec = total_forces;
    else
        vec = [total_forces; 2*alpha*(A'*vn-b1)+(alpha^2)*(A'*un-(b0+b1*t))];
    end
 
    sol = coef_inv * vec;
    dotvn = sol(1:coor_row*npoints);
    dotq = [dotun; dotvn];
end