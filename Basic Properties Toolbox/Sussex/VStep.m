function out = VStep(clampfile, steps, holding, capacitance, N)
%% Decription
% Voltage step protocol analysis

% adapted from the analysis performed in Telezhkin et al., 2016

% @OGSteele, 2022


% example use
%       [file,path] = uigetfile('*.ma');
%       clampfile = fullfile(path,file);
%       steps = [-120,80];
%       holding = -65;
%       capacitance = 220; 
%       N = 10;
%       output = VStep(clampfile, steps, holding, capacitance, N);

% input
%       S = ephysIO output structure of first 'ClampX.ma' recording
%       in a current step sequence, assuming the the rest are in the same
%       directory structure
%       clampfile = full file path of S
%       steps = [min,max] Vm step vals
%       holding = holding potential in mV
%       capacitance = cell capacitance in pF
%       N = P/N fraction (10 by default)

% output
%       

%% Dev Mode
%path = 'D:\Oli Steele\Documents\GitHub\Analysis Scripts\Basic Properties Toolbox\Sussex\exampleData\VoltageStep\OGS_ActInact_Raw_000\000';
%clampfile = [path '\Clamp1.ma'];
%S = ephysIO(clampfile);
%steps = [-120,80];
%holding = -65;
%capacitance = 220; 
%N = 10; 
%output = VStep(S,clampfile,steps,holding);

%% TO DO 
% calculate actual ENa in our experiments
%% Clear the environment
close all

%% Organise data and plot

w = warning ('off','all');
warning(w)

S = ephysIO(clampfile);

% Split into time and waves
Time = S.array(:,1);
Waves = S.array(:,2:end);

% cd to path
splitPath = split(clampfile,filesep);
newPath = char(string(join(splitPath(1:end-3),"\")));
cd(newPath)


%% run analysis or not?
% gives the user the option to abort if the data looks horrendous

% plot the trace
fh = figure();
plot(Time,Waves*1000,'color','black')
box off; set(gcf,'color','white'); set(gca,'linewidth',2)
ylabel('Amplitude (A)'); xlabel('Time (s)')
title('Voltage Step Waveform')

% ask the user 
dlgTitle = 'Run Analysis';
dlgQuestion = 'Would you like to run the Current Step Analysis?';
run = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');

close(fh) % close figure

if run == "Yes"
    % run analysis
    
    %% Perform P/N subtraction
    % ask the user 
    dlgTitle    = 'P/N Subtraction';
    dlgQuestion = 'Is there an associated P/N subtraction recording?';
    subtracted = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
    
    if subtracted == "Yes" 
        [subfile, subpath] = uigetfile('*.*');
        cd(subpath)
        sub = ephysIO(subfile);
        subarray = sub.array(:,2:end);
        fsub = (filter1(subarray,Time, 0, 500))*N; % heavily filtered trace to be subtracted
        Waves = Waves-fsub; % P/N subtracted trace
    elseif subtracted == "No"
        % carry on with the analysis anyway, noting that the recording was
        % not leak subtracted
        subpath = "N/A";
    end
    
    cd ../.. % change back two directories
    
    % create wavelength and Hz variables for ease later
    wavelength = Time(end);
    Hz = size(S.array,1)/wavelength;
    
    
    %% Set cursors
    % NaA start and finish (variable)
    figure;
    plot(Time,Waves,'color','black')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (A)');
    title('Voltage Step Waveform')
    title('highlight the initial sodium spikes')
    [A,~] = ginput(2);
    xline(A(1),'--b','linewidth',1,'HandleVisibility','off')
    xline(A(2),'--b','HandleVisibility','off','linewidth',1)
    % K start and finish (consistent)
    k1 = 7990;
    k2 = 9980;
    % NaI start and finish (consistent)
    i1 = 10005;
    i2 = 10505;
    close(gcf) % close the figure
    
    %% Plot the trace   
    % plot the overall wave above the waveform
    trace_fig = figure;
    subplot(5,1,(1:3))
    plot(Time,Waves*1000,'color','black','HandleVisibility', 'off')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    ylabel('Amplitude (A)');
    title('Voltage Step Waveform')
    ax = gca; xax = ax.XAxis; set(xax,'visible','off')
    
    % plot the regions onto the waveform
    xline(A(1),'--b'); xline(A(2),'--b','HandleVisibility','off');
    xline(i1/Hz,'--g'); xline(i2/Hz,'--g','HandleVisibility','off');
    xline(k1/Hz,'--r'); xline(k2/Hz,'--r','HandleVisibility','off');
    legend('NaA','K','NaI','linewidth',1)
    
    %% plot the waveform 
    Vm = linspace(steps(1),steps(2),size(Waves,2))';
    Vm_waveform = zeros(size(Waves,1),size(Waves,2)) + holding; % sets the flat line at -65 mV (the holding potential) 
    for i = 1:size(Vm,1)
        Vm_waveform(round(0.05*Hz):round(0.25*Hz),i) = Vm(i); % sets the steps
        Vm_waveform(round(0.251*Hz):round(0.35*Hz),i) = 0;
    end
    subplot(5,1,(4:5)); plot(Time,Vm_waveform(:,:),'linewidth',1,'color','black')
    box off; set(gca,'linewidth',2); set(gcf,'color','white');
    xlabel('Time (s)'); ylabel('Membrane Potential (mV)'); ylim([(steps(1)-50),(steps(2)+50)])
    ax = gca; yax = ax.YAxis; set(yax,'TickDirection','out')
    sgtitle('Waveform')
    
    % plot the regions onto the waveform
    xline(A(1),'--b'); xline(A(2),'--b','HandleVisibility','off');
    xline(i1/Hz,'--g'); xline(i2/Hz,'--g','HandleVisibility','off');
    xline(k1/Hz,'--r'); xline(k2/Hz,'--r','HandleVisibility','off');
    
    %% Na Activation values
    NaA_fig = figure; 
    subplot(1,2,1)
    plot(...
        Time((round((A(1)*Hz)-100)):(round((A(2)*Hz)+100))),...
        Waves((round((A(1)*Hz)-100)):(round((A(2)*Hz)+100)),:),...
        'color','black','HandleVisibility','off')
    xline(A(1),'--b'); xline(A(2),'--b','HandleVisibility','off');
    NaApeakval = zeros(1,size(Waves,2));
    ind = zeros(1,size(Waves,2));
    for i = 1:size(Waves,2)
        [NaApeakval(i),ind(i)] = min(Waves(round(A(1)*Hz):round(A(2)*Hz),i));
        hold on; plot((ind(i)+A(1)*Hz)/Hz,NaApeakval(i),'ob'); hold off
    end
    ylabel('Amplitude (A)'); xlabel('Time (s)');
    set(gca,'linewidth',2); set(gcf,'color','white'); box off
    legend('Analysis Region','Peak Amplitude','linewidth',1,'location','southeast')
    
    subplot(1,2,2)
    plot(Vm,NaApeakval,'-ob')
    ylabel('Amplitude (A)'); xlabel('Membrane Potential (Vm)');
    set(gca,'linewidth',2); set(gcf,'color','white'); box off
    sgtitle('NaA')
    
    %% Potassium values
    K_fig = figure; 
    subplot(1,2,1)
    plot(Time(k1-100:k2+100),Waves(k1-100:k2+100,:),'color','black','HandleVisibility','off')
    xline(k1/Hz,'--r'); xline(k2/Hz,'--r','HandleVisibility','off');
    Kmeanval = zeros(1,size(Waves,2));
    for i = 1:size(Waves,2)
        Kmeanval(i) = mean(Waves(k1:k2,i));
        hold on; plot((((k2-k1)/2)+k1-100)/Hz,Kmeanval(i),'or'); hold off
    end
    ylabel('Amplitude (A)'); xlabel('Time (s)');
    set(gca,'linewidth',2); set(gcf,'color','white'); box off
    legend('Analysis Region','Mean Amplitude','linewidth',1,'location','southeast')
    
    subplot(1,2,2)
    plot(Vm,Kmeanval,'-or')
    ylabel('Amplitude (A)'); xlabel('Membrane Potential (Vm)');
    set(gca,'linewidth',2); set(gcf,'color','white'); box off
    sgtitle('K')
    
    %% Na Inactivation values
    NaI_fig = figure; 
    subplot(1,2,1)
    plot(Time(i1-100:i2+100),Waves(i1-100:i2+100,:),'color','black','HandleVisibility','off')
    xline(i1/Hz,'--g'); xline(i2/Hz,'--g','HandleVisibility','off');
    % preallocate for speed
    NaIpeakval = zeros(1,size(Waves,2));
    for i = 1:size(Waves,2)
        [NaIpeakval(i),ind(i)] = min(Waves(i1:i2,i));
        hold on; plot((ind(i)+i1)/Hz,NaIpeakval(i),'og'); hold off
    end
    ylabel('Amplitude (A)'); xlabel('Time (s)');
    set(gca,'linewidth',2); set(gcf,'color','white'); box off
    legend('Analysis Region','Peak Amplitude','linewidth',1,'location','southeast')
    
    subplot(1,2,2)
    plot(Vm,NaIpeakval,'-og')
    ylabel('Amplitude (A)'); xlabel('Membrane Potential (Vm)');
    set(gca,'linewidth',2); set(gcf,'color','white'); box off
    sgtitle('NaI')
    
    %% Na/K Current Density Plots
    % convert to pA and divide by pF to calculate the current density in
    % in pA/pF
    NaK_ID_fig = figure;
    NaID = (NaApeakval*1e12)/capacitance;
    KID = (Kmeanval*1e12)/capacitance;
    plot(Vm,NaID,'-ob')
    hold on; plot(Vm,KID,'-or')
    ax = gca;
    ax.XAxisLocation = 'origin'; xlabel('Membrane Potential (mV)')
    ax.YAxisLocation = 'origin'; ylabel('Current Density (pA/pF)')
    box off; set(gcf,'color','white'); set(gca,'linewidth',2)
    legend('Sodium Activation','Potassium','location','southwest','linewidth',1);
    sgtitle('Na/K Current Density')
    %% Sodium 'Availability' Plots
    % mean fractional conductance plots
    % "For activation and inactivation curves, conductance (G) was 
    % calculated by dividing current by the appropriate driving force 
    % (Vc − ENa), where Vc = command potential and ENa = +66.7 mV."
    %                                         ~ Telezhkin et al., 2016
    
    % calculate ENa here properly later
    ENa = 66.7; 
    % convert to conductance (G)
    GNaA =  zeros(1,size(Vm,1));
    GNaI =  zeros(1,size(Vm,1));
    for i = 1:size(Vm,1)
        GNaA(i) = NaApeakval(i) / (Vm(i) - ENa);
        GNaI(i) = NaIpeakval(i) / (Vm(i) - ENa);
    end
    % convert to fractional conductance (G/GMax)
    Frac_GNaA = GNaA./max(GNaA(1:25));
    Frac_GNaI = GNaI./max(GNaI(1:25));
    % plot
    NaFracG_fig = figure; plot(Vm,Frac_GNaA,'-ob')
    hold on; plot(Vm,Frac_GNaI,'-og')
    %hold on; plot(-68, 0.1, '^k') % temporarily remove the membrane
    %potential bit
    xlim([-120 0]); box off; set(gcf,'color','white')
    set(gca,'linewidth',2); xlabel('Membrane Potential (mV)');
    ylabel('Fractional Conductance (G/GMax)');
    legend('Sodium Activation','Sodium Inactivation','linewidth',1,'location','west');
    %legend('Sodium Activation','Sodium Inactivation', 'Membrane Potential' ,'linewidth',1,'location','west');
    sgtitle('Sodium Mean Fractional Conductance')

    
    %% Conditions
    % What Genotype was the animal?
    dlgTitle    = 'Genotpye';
    dlgQuestion = 'What Genotype was this animal?';
    genotype = questdlg(dlgQuestion,dlgTitle,'APOE3','APOE4','APOE3');

    % Was the recording exposed to ketamine?
    dlgTitle    = 'Ketamine';
    dlgQuestion = 'Was this recorded in the presence of ketamine?';
    ketamine = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
 
    % Was the recording exposed to memantine?
    dlgTitle    = 'Memantine';
    dlgQuestion = 'Was this recorded in the presence of meamntine?';
    memantine = questdlg(dlgQuestion,dlgTitle,'Yes','No','No');
    
    %% Create the output
    % Output Generation
    out.Genotype = genotype;
    out.Ketamine = ketamine;
    out.Memantine = memantine;
    out.rawFilepath = clampfile;
    out.subfilepath = subpath;
    out.leakSubtracted = subtracted;
    out.raw = S;
    out.sub = sub;
    out.Time = Time;
    out.Waves = Waves;
    out.Vm = Vm;
    out.Cap = capacitance;
    out.NaAPeak = NaApeakval;
    out.KMean = Kmeanval;
    out.NaIPeak = NaIpeakval;
    out.NaID = NaID;
    out.KID = KID;
    out.GNaA = GNaA;
    out.GNaI = GNaI;
    
    % create output directory
    dirname = [char(splitPath(end-2)),'_output'];
    mkdir(dirname); cd(dirname);
    
    % save output & figures
    save('out.mat','out')
    savefig(trace_fig,'trace_fig')
    savefig(NaA_fig, 'NaA_fig')
    savefig(K_fig,'K_fig')
    savefig(NaI_fig,'NaI_fig')
    savefig(NaK_ID_fig,'NaK_ID_fig')
    savefig(NaFracG_fig,'NaFracG_fig')
end