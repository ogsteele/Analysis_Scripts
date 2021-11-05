function [splits, nf_splits] = TP_null_MedianF_function(Param,splits,x,poles)

%% exlcudedata, filter and restitch
% % warning off
warning('off','MATLAB:colon:nonIntegerIndex')
nf_splits = splits;

for i = 1:size(splits,2)
    array = splits(:,i);
    time = 0:5.0000e-05:Param.split_length;
    time = time(1:Param.sample_rate*Param.split_length)';

    % select the range to exclude
    % tf = excludedata(ydata,xdata,'range',[100, 150]);
    range_dp = [x(1), x(1)+(0.035*Param.sample_rate)];
    range_s = range_dp/Param.sample_rate;
    
    tf = excludedata(array,time,'range',range_s); % not excl.
    tx = ~excludedata(array,time,'range',range_s); % excl.
    
    % filter the region, bar exclusions
    yf = medianf(array(tf), time(tf) ,poles); % filtered array
    
    % now add back in the bit you didn't want to be filtered
    full = vertcat(yf(1:x(1)),array(tx),yf((x(1)):end));
    splits(:,i) = full;
end

% bug testing plots
% plot region with exclusion zone
% % plot(time(tf),array(tf))
% % plot the excluded region over the top (another color by default
% % hold on; plot(time(tx),array(tx))
% % title('Region excluded from median filtering')
% % figure; plot(full)
% % 
% % figure; plot(tf); hold on; plot(tx); ylim([-0.5 1.5]); hold on
% % legend('tf','tx')
% % xline(range_dp(1)); hold on; xline(range_dp(2)); hold off
% % legend('tf','tx')
% % title('tf / tx')
