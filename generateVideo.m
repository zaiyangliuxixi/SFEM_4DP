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
     coor_row = size(init_coor,1);
     for t = 0:T_end/288*videoScale:T_end
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

     ylimRangeMin = overall_min(2) - 0.1*diff_max_min(2);
     ylimRangeMax = overall_max(2) + 0.1*diff_max_min(2);  

     xlim(app.UIAxesVideo, [ xlimRangeMin, xlimRangeMax]);
     ylim(app.UIAxesVideo, [ ylimRangeMin, ylimRangeMax]);
     if coor_row ==3
         zlimRangeMin = overall_min(3) - 0.1*diff_max_min(3);
         zlimRangeMax = overall_max(3) + 0.1*diff_max_min(3);  
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
    for t = 0:T_end/288*videoScale:T_end       
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
        app.object4D.draw_individual_app(app.UIAxesVideo,disps); 
         
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