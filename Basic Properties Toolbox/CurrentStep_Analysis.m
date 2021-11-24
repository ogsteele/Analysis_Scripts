function [Rheobase_pA,IR_MOhm,n_APs] = CurrentStep_Analysis(filename)

[I_Step,~,~,fn] = abfload(filename);
I_Step = squeeze(I_Step);
[~,name,~] = fileparts(fn);

% preallocate
n_APs = zeros(size(I_Step,2),1);

for i = 1:size(I_Step,2)
    max_val = max(I_Step(:,i));
    if max_val > 0 
       [events,index] = findpeaks(I_Step(:,i),'minPeakDistance',100,'minPeakHeight',0);
       n_APs(i) = size(events,1);
       if index > 14500
           n_APs(i) = 0;
       end
    else
       n_APs(i) = 0;
    end
end

pA = [-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]';
pA_text = [{'-10 pA'},{'0 pA'},{'+10 pA'},{'+20 pA'},{'+30 pA'},{'+40 pA'},{'+50 pA'},...
    {'+60 pA'},{'+70 pA'},{'+80 pA'},{'+90 pA'},{'+100 pA'},{'+110 pA'},{'+120 pA'},...
    {'+130 pA'},{'+140 pA'},{'+150 pA'},{'+160 pA'},{'+170 pA'},{'+180 pA'}]';
%T = table(pA,n_APs);

first = find(n_APs,1);
    if first = []
        Rheobase_pA = NaN;
    else
        Rheobase_pA = pA(first);
    end

x = linspace(0,size(I_Step,1)/20000,size(I_Step,1));
figure 
plot(x,I_Step(:,first))
hold on
plot(x,I_Step(:,1))
plot(x,I_Step(:,2))
title("Rheobase Plot of recording " + name + ".abf")
ylabel('Membrane potential (mV)');
xlabel('Time (s)');
xlim([0,3]);
box off
hold off
legend(cell2mat(pA_text(first)),cell2mat(pA_text(1)),cell2mat(pA_text(2)),...
    'AutoUpdate','off')

[x2] = getpts;
    a = round(x2(2))*20000;
    b = round(x2(1))*20000;
IR_MOhm = (((mean(I_Step(b:a,2) - mean(I_Step(b:a,1))))/1000)/(1e-11))/1000000;
    %output.IR = IR;
%for j = 1:size(I_Step,2)
    %figure
    %plot(I_Step(:,j))
%end
end