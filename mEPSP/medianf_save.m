function medianf_save
% Median filter code to apply a user set number ( or 15 pole by default )
% of poled median filter and then resave as .phy recording

% load in file w/ ephysIO
[filename,path] = uigetfile('*.tdms'); % ui file name
disp(append('User selected ', filename, ' ...'))
disp(append('Loading ', filename, ' with epysIO.m ...'))
S=ephysIO(fullfile(path,filename)); % load file selected
t=S.array(:,1); % time
y=S.array(:,2)*0.01*1000; % array

% set number of poles (commented out user selection)

% prompt = {...
%     'Number of poles to filter',...
%     };
% defaults = {'15'};
% dlgtitle = 'Filter Settings';
% dims = [1 50];
% poles = inputdlg(prompt,dlgtitle,dims,defaults);
% poles_num = str2double(poles); % convert to number
% poles_str = append('_',sprintf('%.0f',poles_num)); % convert to string

% comment out if the above is there
poles_num = 3; % convert to number
poles_str = append('_',sprintf('%.0f',poles_num)); % convert to string

% set filetype to save to
filetype = '.phy'; % file type to save to

% apply median filter and plot
disp(append('Applying ',string(poles_num), ' pole median filter to ', filename, ' ...'))
yf=medianf(y,t,poles_num); % filtered array
disp(append('Plotting an overlay of the filtered and raw traces ...'))
plot(t,y,'r');hold on;plot(t,yf,'b'); legend('raw','filtered');hold off; xlabel('Time (s)'); ylabel('Vm (mV)')% plot traces

% save filtered file in current working directory
cd(path) % change directory to location of file
new_filename = append(filename(1:end-5),'_medianf',poles_str,filetype); % create filenames
disp(append('Saving ',new_filename, ' in current directory ...'))
ephysIO(new_filename,[t,yf],'s','mV') % save as new filename
end