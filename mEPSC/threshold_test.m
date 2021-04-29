% tidy the environment
clear
close all

% load ensemble averages
[file,directory] = uigetfile('*.phy');
cd(directory)
S = ephysIO(file);

% apply a further 1000 hz filter
YF = S.array;
for s = 2:size(S.array,2)
t = (0:size(YF(:,s),1)-1)';
t = ((0:size(YF(:,s),1)-1)')./20000;
YF(:,s) = filter1(YF(:,s), t, 0, 1000);
end

% plot filtered traces and 3SD of the mean overlays
for i = 2:size(S.array,2)
figure; plot(YF(:,i)); ylim([-10e-12 10e-12])
mean_trace = mean(YF(:,i)); stdev_trace = std(YF(:,i)); threshold_dev = 3.5;
yline(mean_trace);
yline(mean_trace + (threshold_dev * stdev_trace));
yline(mean_trace - (threshold_dev * stdev_trace));

    % does it cross the threshold
    if min(YF(1:2000,i)) < (mean_trace - (threshold_dev * stdev_trace))
        event(i-1) = 1;
    else
        event(i-1) = 0;
    end
    
dim = [.2 .5 .3 .3];
    if event(i-1) == 1
        str = 'Event detected';
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    else
    end
    events_detected = sum(event);
end

% happy with selection decision tree
answer_happy = questdlg('Happy with selection?', ...
'Happy?', ...
'Yarp','Narp','dft');

if answer_happy == 'Yarp'
    % if yes, concatenate array of event-positive only traces
    positive = [];
    for ii = 1:size(event,2)
        if event(ii) == 1
            temp = S.array(:,ii+1);
            positive = [positive temp];
        end 
    end 
    % if no, then bring up selection dialogue
elseif answer_happy == 'Narp'
    for ii = 1:size(event,2)
        prompt = {'1 = event, 0 = not event'};
        dlgtitle = num2str(ii);
        dims = [1 35];
        definput = {num2str(event(ii))};
        answer_selection(ii) = inputdlg(prompt,dlgtitle,dims,definput);        
    end
    answer_selection = str2double(answer_selection); % convert to double
    % now concatenate arrays
    positive = [];
        for ii = 1:size(answer_selection,2)
            if answer_selection(ii) == 1
             temp = S.array(:,ii+1);
                positive = [positive temp];
            end 
        end  
        
end

% display percent positive
percent_positive = (size(positive,2)/(size(S.array,2)-1))*100
num_pos = size(positive,2)
num_total = size(S.array,2)-1

% add time as the first and save as ensemble average
time = 5e-5*[1:size(positive,1)]';

% prompt user to name the output
prompt = {'Name the concatenated output'};
def = {'Genotype Gender Receptor'};
dlgtitle = 'output name';
dims = [1 50];
name = inputdlg(prompt,dlgtitle,dims,def);
% create the output in the structure
filename = append(name,'_ensemble_averages_positive.phy');
% save the ensemble average
ephysIO(char(filename),[time positive],'s','A') 
