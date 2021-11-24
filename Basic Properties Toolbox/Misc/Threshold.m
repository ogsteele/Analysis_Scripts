function [] = Threshold(filename2,sweep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CS = abfload(filename2);
CS10 = smooth(CS(:,sweep));
time = (1:103224)';
plot (time,CS10);
[x,y] = getpts;
    peak = round(x);
CS10win = CS10(1:peak);
timewin = (1:peak);
TDiff = smooth(diff(diff(diff(CS10win))));
TimeDiff = (1:(peak-3))';
plot(timewin,CS10win,TimeDiff,TDiff);
[x,y] = getpts;
    pre = round(x(1));
    post = round(x(2));
TDiff2 = TDiff(pre:post);
[a,b] = max(TDiff2);
CS10win2 = CS10(pre:post);
Threshold = CS10win2(b);
end

