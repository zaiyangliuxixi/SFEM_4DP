function [] = generateVideo(app)    
     M = [];                                           % 
     T_end = app.SimulationTime;                       % simulation time
     time = app.result_time;                           % time array
     q = app.result_q;                                 % displacement
     npoints = app.object4D.numNodalPoints;
     sizeScale = app.SizeScale;        
     timeScale = app.TimeScale;                        % real time v.s simulation time
     timeUnit = app.TimeUnitDropDown.Value;
     videoScale = app.VideoScale;
     num_pic = 0;                                      %number of pictures
     % to calcute the maximal and mininal coordinates
     max_per_time = [];
     min_per_time = [];
     init_coor =[];
     for i = 1:npoints
         init_coor = [init_coor, app.object4D.NodalPoints(i).Coordinates];
     end
     % to obtain top points and bottom points
     npoints_m = app.numNode_x_2D;
     bottom_point_init_coor = zeros(2,npoints_m);
     top_point_init_coor    = zeros(2,npoints_m);
    
     for i = 1:npoints_m
         bottom_point_init_coor(:,i) = app.object4D.NodalPoints(1,i).Coordinates;
         top_point_init_coor(:,i) = app.object4D.NodalPoints(1,npoints-npoints_m+i).Coordinates;
     end
     coor_row = size(init_coor,1);
     t_values = linspace(0, T_end, round(288 / videoScale) + 1);
     for t = t_values
            %fprintf("time %f\n", t);
            index = nearest_index(time, t);
            disps = reshape(q(index, 1:npoints * coor_row), [coor_row, npoints]);        
            cur_coor = init_coor + disps;        
            % to get the max and min at current time
            max_current_time = max(cur_coor, [], 2);  % max of per row
            min_current_time = min(cur_coor, [], 2);  
            max_per_time = [max_per_time,max_current_time];
            min_per_time = [min_per_time,min_current_time];
     end
     overall_max = max(max_per_time,[],2);  
     overall_min = min(min_per_time,[],2);  
     diff_max_min = abs(overall_max - overall_min);

     axis(app.UIAxesVideo, 'equal');     %放在limit前面
     xlimRangeMin = overall_min(1) - 0.1*diff_max_min(1);
     xlimRangeMax = overall_max(1) + 0.1*diff_max_min(1);

     ylimRangeMin = overall_min(2) - 0.8*diff_max_min(2);
     ylimRangeMax = overall_max(2) + 0.8*diff_max_min(2);  

     xlim(app.UIAxesVideo, [ xlimRangeMin, xlimRangeMax]);
     ylim(app.UIAxesVideo, [ ylimRangeMin, ylimRangeMax]);
     if coor_row ==3
         zlimRangeMin = overall_min(3) - 0.5*diff_max_min(3);
         zlimRangeMax = overall_max(3) + 0.5*diff_max_min(3);  
         zlim(app.UIAxesVideo, [ zlimRangeMin, zlimRangeMax]);
     end
  
     % obtain view of app.UIAxes
     [az, el] = view(app.UIAxes);
     % same as  UIAxes 
     view(app.UIAxesVideo, [az, el]);
   
    % set tickLabel
    xticks = get(app.UIAxesVideo, 'XTick');
    yticks = get(app.UIAxesVideo, 'YTick');
    zticks = get(app.UIAxesVideo, 'ZTick');

    epsilon = 1e-10;
    xticklabels = arrayfun(@(x) num2str(x/sizeScale), xticks, 'UniformOutput', false);
    yticklabels = arrayfun(@(y) num2str(y/sizeScale), yticks, 'UniformOutput', false);
    zticklabels = arrayfun(@(z) num2str(z/sizeScale), zticks, 'UniformOutput', false);
    
    xticklabels(abs(xticks) < epsilon) = {'0'};
    yticklabels(abs(yticks) < epsilon) = {'0'};
    zticklabels(abs(zticks) < epsilon) = {'0'};
    
    set(app.UIAxesVideo, 'XTickLabel', xticklabels);
    set(app.UIAxesVideo, 'YTickLabel', yticklabels);
    set(app.UIAxesVideo, 'ZTickLabel', zticklabels);
  
    % Generating figures 
    for t = t_values      
         % Check for Cancel button press
        if app.VideoProgressDlg.CancelRequested
           return
        else
            progress = t / T_end;   
            app.VideoProgressDlg.Value = progress;
            app.VideoProgressDlg.Message = sprintf('Video generating: %.0f%%', progress * 100);
        end
        % clear fig every time to avoid overlap   
        cla(app.UIAxesVideo);
        index = nearest_index(time, t);
        disps = reshape(q(index, 1:npoints * coor_row), [coor_row, npoints]);
        % app.UIAxesVideo 
        hold(app.UIAxesVideo, 'on');
        %  draw figure at every time
        app.object4D.draw_individual_app(app.UIAxesVideo,disps); 
         
        % add bending angle change
        value = app.TypeDropDown.Value;   
        if strcmp(value, '2D') && app.isAngleOutput == 1
            top_point_init_start = top_point_init_coor(:,1);
            top_point_init_end = top_point_init_coor(:,end);
            final_point = top_point_init_end  + reshape(q(index, 2*npoints-1:2*npoints), [2, 1]);
            drawBendingAngle(top_point_init_start,top_point_init_end, final_point,app)
        elseif strcmp(value, '2D') && app.isAngleOutput == 2
            bottom_points_cur_coor = bottom_point_init_coor + reshape(q(index, 1:2*npoints_m), [2, npoints_m]);
            top_points_cur_coor = top_point_init_coor + reshape(q(index, 2*(npoints-npoints_m)+1:2*npoints), [2, npoints_m]);
            mid_points_cur_coor = (bottom_points_cur_coor + top_points_cur_coor)./2;
            %plot(app.UIAxesVideo, mid_points_cur_coor(1,:), mid_points_cur_coor(2,:),'bo');
            drawFittingArc(mid_points_cur_coor,app);
        else
        end

        real_time = t*timeScale;
        epsilon = 1e-8;
        if real_time <epsilon
           time_show = real_time;
           timeUnit_show = '';
        else
           [time_show, timeUnit_show] =updateTimeUnit(real_time, timeUnit);
        end 
        title(app.UIAxesVideo,['Time: ' num2str(time_show,"%3.2f") timeUnit_show]);     
         try
                %save current pic
                if  strlength(app.SelectedFolderPath) ~= 0
                    filename = fullfile(app.SelectedFolderPath, sprintf('deformation_%04d.png', num_pic));
                    exportgraphics(app.UIAxesVideo, filename, 'Resolution', 300);
                    % read pic from file
                    img = imread(filename);                    
                    % 调整图像大小与第一个帧相同
                    if isempty(M)
                        frame_size = size(img);  % 使用第一张图像的尺寸作为标准
                    else
                        img = imresize(img, [frame_size(1), frame_size(2)]);
                    end
                    frame = im2frame(img);                    
                    % 将帧加入到视频中
                    M = [M, frame];
                    num_pic = num_pic +1;
                end
          catch ME
             msgbox(['Error occurred: ', ME.message], 'Error', 'error','modal');
             return;
          end

    end
    % save video
    try
        if  strlength(app.SelectedFolderPath) ~= 0
            % 将 videoScale 放大10倍并取整
            videoScaleInteger = round(videoScale * 10);
            filename_video = fullfile(app.SelectedFolderPath, sprintf('video_%d.mp4', videoScaleInteger));
            v = VideoWriter(filename_video, 'MPEG-4');
            v.Quality = 100;  % 设置视频质量为最高
            open(v);
            writeVideo(v, M);
            close(v);
            uialert(app.UIFigure, 'Video generated successfully!', 'Success');
        end       
    catch ME
         msgbox(['Error occurred: ', ME.message], 'Error', 'error','modal');
         return;
    end
end


function [time_show, timeUnit_show] = updateTimeUnit(real_time, timeUnit)

    switch timeUnit
        case 's'
            if real_time >= 3600  % 大于等于 3600秒，转换为小时
                time_show = real_time / 3600;
                timeUnit_show = 'h';
            elseif real_time >= 60  % 大于等于 60秒，转换为分钟
                time_show = real_time / 60;
                timeUnit_show = 'min';
            else
                time_show = real_time;
                timeUnit_show = 's';
            end
        case 'min'
            if real_time >= 60  % 大于等于 60分钟，转换为小时
                time_show = real_time / 60;
                timeUnit_show = 'h';
            elseif real_time < 1  % 小于 1分钟，转换为秒
                time_show = real_time * 60;
                timeUnit_show = 's';
            else
                time_show = real_time;
                timeUnit_show = 'min';
            end
        case 'h'
            if real_time < 1/60  % 小于 1分钟，转换为秒
                time_show = real_time * 3600;
                timeUnit_show = 's';
            elseif real_time < 1  % 小于 1小时，转换为分钟
                time_show = real_time * 60;
                timeUnit_show = 'min';
            else
                time_show = real_time;
                timeUnit_show = 'h';
            end
    end
end


function index = nearest_index ( vec, value )
    sz = size(vec); n = sz(1);
    i = 1; j = n;
    
    if vec(i) >= value
        index = i;
        return;
    end
    if vec(j) <= value
        index = j;
        return;
    end
    
    while (i < j-1)
        k = floor((i+j)/2);
        if vec(k) <= value
            i = k;
        end
        if vec(k) >= value
            j = k;
        end
    end
    
    if value - vec(i) < vec(j) - value
        index = i;
    else
        index = j;
    end
end

function [] =drawFittingArc(mid_points, app)
    x = mid_points(1, :);
    y = mid_points(2, :);
      
    % objective function
    circleFit = @(params) sum((sqrt((x - params(1)).^2 + (y - params(2)).^2) - params(3)).^2);
    
    % initial guess [xc, yc, r]
    x0 = mean(x); 
    y0 = mean(y);
    r0 = mean(sqrt((x - x0).^2 + (y - y0).^2));
    initial_guess = [x0, y0, r0];
    
    % use fminsearch fit a circle
    optimal_params = fminsearch(circleFit, initial_guess);
    
    % get center and radius
    xc = optimal_params(1);
    yc = optimal_params(2);
    r = optimal_params(3);
    
    % calculated start angle and end angle
    theta_start = atan2(y(1) - yc, x(1) - xc);
    theta_end = atan2(y(end) - yc, x(end) - xc);
    
    % polar coordinate
    theta_arc = linspace(theta_start, theta_end, 100);
    x_arc = xc + r * cos(theta_arc);
    y_arc = yc + r * sin(theta_arc);
    
    % slopes
    m1 = -(x(1) - xc) / (y(1) - yc); % 起点切线斜率
    m2 = -(x(end) - xc) / (y(end) - yc); % 终点切线斜率
    
    % line equation 
    % L1: y = m1 * (x - x(1)) + y(1)
    % L2: y = m2 * (x - x(end)) + y(end)
    
    % obtain intersection poins
    A = [-m1, 1; -m2, 1]; % 
    b = [y(1) - m1 * x(1); y(end) - m2 * x(end)];
    intersection = A\b; % 
    
    % calculate bending angle
    angle_rad = abs(theta_end - theta_start);
    angle_deg = rad2deg(angle_rad);
    

    % get a tangent line from start point
    direction_x_start = intersection(1) - x(1);
    direction_y_start = intersection(2) - y(1);
    x_line_start = [x(1), x(1) + 2 * direction_x_start];
    y_line_start = [y(1), y(1) + 2 * direction_y_start];
  
    % 
    x_line_end = [intersection(1), x(end)];
    y_line_end = [intersection(2), y(end)];

    % get ceter line direction between two lines
    vector_start = [direction_x_start, direction_y_start]; 
    vector_end = [x(end) - intersection(1), y(end) - intersection(2)]; 
    center_line_vector = (vector_start / norm(vector_start) + vector_end / norm(vector_end)); 
    center_line_vector = center_line_vector / norm(center_line_vector); 
    
  
    offset = 0.5*sqrt(direction_x_start^2+direction_y_start^2);
    x_angle_label = intersection(1) + offset * center_line_vector(1);
    y_angle_label = intersection(2) + offset * center_line_vector(2);
    

    plot(app.UIAxesVideo,x_line_start, y_line_start, 'k--', 'LineWidth', 1.5);     
    plot(app.UIAxesVideo,x_arc, y_arc, 'k-', 'LineWidth', 2); 
    plot(app.UIAxesVideo,x_line_end, y_line_end, 'k--', 'LineWidth', 1.5);
    plot(app.UIAxesVideo,intersection(1), intersection(2), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    text(app.UIAxesVideo,x_angle_label, y_angle_label, sprintf('%.2f°', angle_deg), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 16, 'Color', 'black');  
end


function [] = drawBendingAngle(init_start, init_end,final_point,app)

 
    direction_vector = init_end - init_start;
    extended_point = init_end + direction_vector;
    
    % calculate bending angle
    vec1 = init_end - init_start;
    vec2 = final_point - init_start;
    angle_rad = acos(dot(vec1, vec2) / (norm(vec1) * norm(vec2)));
    angle_deg = rad2deg(angle_rad);
    
    %
    angle_bisector = (vec1 / norm(vec1) + vec2 / norm(vec2)) / 2;
    angle_bisector = angle_bisector / norm(angle_bisector); % nomalization
    offset_multiplier = 0.5; %
    offset_distance = offset_multiplier * norm(vec1);
    label_point = init_start + offset_distance * angle_bisector;
    
    % 
    plot(app.UIAxesVideo,[init_start(1), init_end(1)], [init_start(2), init_end(2)], '--k', 'LineWidth', 1.5);
    plot(app.UIAxesVideo,[init_end(1), extended_point(1)], [init_end(2), extended_point(2)], '--k', 'LineWidth', 1.5);
    plot(app.UIAxesVideo,[init_start(1), final_point(1)], [init_start(2), final_point(2)], '--k', 'LineWidth', 1.5);
    
    % 
    scatter(app.UIAxesVideo,[init_start(1), final_point(1)], [init_start(2), final_point(2)], 25, 'filled','k');
    text(app.UIAxesVideo,label_point(1), label_point(2), sprintf('%.2f°', angle_deg), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 16, 'Color', 'black');  

end