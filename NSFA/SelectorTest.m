%% NSFA Selector

% load file
[file,folder] = uigetfile('*.*');
S = ephysIO(fullfile(folder,file));
Time = S.array(:,1);
Array = S.array(:,2:end);

figure;
LineList = plot(Time,Array);
set(LineList, 'ButtonDownFcn', {@myLineCallback, LineList});


function myLineCallback(LineH, EventData, LineList)
%disp(LineH);                    % The handle
%disp(get(LineH, 'YData'));      % The Y-data
disp(find(LineList == LineH));  % Index of the active line in the list
set(LineList, 'LineWidth', 0.5);
set(LineH,    'LineWidth', 2.5);
uistack(LineH, 'top');  % Set active line before all others
end

% What do we need the script to do?
% logic should be ...
    % plot trace
    % ginput to square off the region you want the peaks to be within
    % series of logical statements
        % does this wave have a peak in that square region?
        % if yes, merge into new recording
        % if no, discard
