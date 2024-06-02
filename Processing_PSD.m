clear all
close all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))

matfiles = ["c7.20240417.111204.calibration.cvsa_lbrb_filt.mat","c7.20240417.111557.calibration.cvsa_lbrb_filt.mat"];
selectedFreq = [4:2:48];
%% INFO
wlength = 0.5; % seconds. Length of the external window
pshift = 0.25; % seconds. Shift of the internal windows
wshift = 0.0625; % seconds. Shift of the external window
samplerate = 512;   % curr_h.SampleRate (??)
mlength = 1; % seconds

%% CONCATENAZIONE
PSD = []; events = struct('POS', [], 'DUR', [], 'TYP', [], 'SampleRate', [], 'info',[]);
for i=1:length(matfiles)
    cdata = load(matfiles(i));
    cpsd = cdata.psd;
    cevents = cdata.event;
    events.POS = cat(1,events.POS,cevents.POS + size(PSD,1));
    events.DUR = cat(1,events.DUR,cevents.DUR);
    events.TYP = cat(1,events.TYP,cevents.TYP);
    events.SampleRate = cdata.samplerate;   %(??)
    PSD = cat(1,PSD,cpsd);
end
% Create Activity matrix
feedb_dur = events.DUR(events.TYP == 781);
feedb_pos = events.POS(events.TYP == 781);
fix_pos = events.POS(events.TYP == 786);
fix_dur = events.DUR(events.TYP == 786);
cue_pos = events.POS(events.TYP == 730 | events.TYP == 731);
cue_dur = events.DUR(events.TYP == 730 | events.TYP == 731);

nwindows = size(PSD,1);
nfreq = size(PSD,2);
nchannels = size(PSD,3);
trials = length(feedb_pos); %ntrials=sum(events.TYP==781)

% Siamo interessati al periodo dalla fixation alla fine del feedback per le
% due MI task (730,731 in cue)
Ck = zeros(nwindows,1);
Tk = zeros(nwindows,1);
TrialStart = NaN(trials,1);
TrialStop = NaN(trials,1);
FixStop = NaN(trials,1);
for trId=1:trials
    cstart = fix_pos(trId);
    cstop = feedb_pos(trId) + feedb_dur(trId) - 1;
    Ck(cstart:cstop) = events.TYP(cue_pos(trId)==events.POS);
    Tk(cstart:cstop) = trId;

    TrialStart(trId) = cstart;
    TrialStop(trId) = cstop;
    FixStop(trId) = cstart+fix_dur(trId)-1;
end

trial_dur = max(TrialStop-TrialStart);
Activity = NaN(trial_dur,nfreq,nchannels,trials);
tCk = zeros(trials,1);
for trId=1:trials
    cstart = fix_pos(trId);
    cstop = cstart + max(trial_dur) - 1;
    Activity(:,:,:,trId) = PSD(cstart:cstop,:,:);
    tCk(trId)=unique(nonzeros(Ck(cstart:cstop)));    
    %posso mettere al posto di cstop
    %la lunghezza di Tk==trId, così prende tutti i sample dei trial senza
    %sforare?
end

% Reference Matrix referred to the fixation period
MaxFix_dur = max(fix_dur); 
Reference = NaN(MaxFix_dur, nfreq, nchannels, trials);
for trId=1:trials
    cstart = fix_pos(trId);
    cstop = cstart+MaxFix_dur - 1;
    Reference(:,:,:,trId) = PSD(cstart:cstop,:,:);
end

% Compute ERD/ERS for each trial 
Baseline = repmat(mean(Reference), [size(Activity,1) 1 1]);
ERD = log(Activity./Baseline);

freqBand = [8 12];
bandIndices = find(selectedFreq >= freqBand(1) & selectedFreq <= freqBand(2));
f_alpha = selectedFreq(bandIndices);

%% IMAGESC??
% figure
% t = linspace(0,trial_dur*info.wshift,trial_dur);
% ChanSelected = [7 9 11];
% tasks = [730 731];
% taskLb = {'Both Hands', 'Both Feet'};
% for taskId=1:length(tasks)
%     for chId=1:length(ChanSelected)
%         subplot(2,3,(taskId-1)*length(ChanSelected)+chId)
%         cERD = mean(ERD(:,:,ChanSelected(chId),tCk==tasks(taskId)),4);
%         imagesc(t,selectedFreq,cERD');
%         set(gca,'YDir','normal');
%         colormap("hot");
%         colorbar;
%         title(['Channel ' channelLb{ChanSelected(chId)} ' | ' taskLb{taskId}]);
%         xlabel('Time [s]');
%         ylabel('Frequency [Hz]');
%     end
%     
% end

%% VISULIZATION (topoplot)
%% Media degli ERD su windows e trials ad una data frequenza sui 39 canali
load("240503\new_chanlocs64.mat")

% Definizione dei canali registrati (esempio)
recorded_channels = {'Fp1','Fp2','F3','Fz','F4','FC1','FC2','C3','Cz','C4','CP1','CP2','P3','Pz','P4','POz','O1','O2','F1','F2','FC3','FCz','FC4','C1','C2','CP3','CPz','CP4','P5','P1','P2','P6','PO5','PO3','PO4','PO6','PO7','PO8','Oz'};
% Trova gli indici dei canali registrati nella struct chanlocs
% ELETTRODI PO5 E PO6 DA ACQUISIZIONE NON SONO PRESENTI, SI DEVE TROVARE IL
% CORRISPETTIVO
recorded_indices = zeros(1,length(recorded_channels));
for i = 1:length(recorded_channels)
    recorded_indices(i) = find(strcmp({chanlocs.labels},recorded_channels{i}));
end

% Estrai i dati ERD per i canali registrati
tasks = [730, 731];
figure(1);
%N.B.: devo distinguere le due task
for taskId=1:length(tasks)
    for fId=1:length(bandIndices)
ERD_task = ERD(:,:,:,find(tCk==tasks(taskId)));
mean_ERD = squeeze(mean(ERD_task, 4)); % Media su dimensione 1 (windows) e 4 (trials)
% Visualizzazione con topoplot per una frequenza specifica
data_for_topoplot = mean_ERD(fId, :);
% Visualizza i dati con topoplot
subplot(2,3,(taskId-1)*length(f_alpha)+fId)
topoplot(data_for_topoplot, chanlocs(recorded_indices), 'electrodes', 'on');
colorbar;
title(sprintf('ERD at %d Hz', f_alpha(fId)));
    end
end

%% Differenza tra i due trial



%% Media ERD sulle finestre temporali dei due tipi di trial
% Estrai i dati ERD per i canali registrati per entrambi i task
% Media degli ERD su tutte le prove
mean_ERD_task1 = mean(ERD(:, :, :, find(tCk==730)), 4); % Media su dimensione 4 (trials)
mean_ERD_task2 = mean(ERD(:, :, :, find(tCk==731)), 4); % Media su dimensione 4 (trials)

% Differenza tra i due task
difference_ERD = mean_ERD_task1 - mean_ERD_task2;
num_windows = size(mean_ERD_task1, 1); % Numero di finestre temporali

% Creare una figura per visualizzare i risultati per ogni window e
% frequenza alfa
for fId=1:length(bandIndices)
    curr_f_index = bandIndices(fId);
    figure(fId);
    for window_index = 1:num_windows
    % Prepara i dati per topoplot per il task 1
        data_for_topoplot_task1 = squeeze(mean_ERD_task1(window_index, curr_f_index, :));

    % Prepara i dati per topoplot per il task 2
        data_for_topoplot_task2 = squeeze(mean_ERD_task2(window_index, curr_f_index, :));

    % Prepara i dati per topoplot della differenza
        data_for_topoplot_difference = squeeze(difference_ERD(window_index, curr_f_index, :));

    % Visualizza i dati con topoplot per il task 1
        subplot(3, num_windows, window_index);
        topoplot(data_for_topoplot_task1, chanlocs(recorded_indices), 'electrodes', 'on');
        colorbar;
        title(sprintf('Task 1 - Window %d', window_index));

    % Visualizza i dati con topoplot per il task 2
        subplot(3, num_windows, num_windows + window_index);
        topoplot(data_for_topoplot_task2, chanlocs(recorded_indices), 'electrodes', 'on');
        colorbar;
        title(sprintf('Task 2 - Window %d', window_index));

    % Visualizza i dati con topoplot della differenza
        subplot(3, num_windows, 2 * num_windows + window_index);
        topoplot(data_for_topoplot_difference, chanlocs(recorded_indices), 'electrodes', 'on');
        colorbar;
        title(sprintf('Difference - Window %d', window_index));
    end
% Imposta un titolo generale per la figura
sgtitle(sprintf('ERD Comparison at %d Hz', curr_f_index));

end

