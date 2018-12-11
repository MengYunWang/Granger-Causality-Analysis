%% Pre-processing EEG data
% The preprocessing consisted of two stages:
% Stage I: preICA and run ICA
%   preICA: Import data> Channel Location > Filter > resample the data>  
%           rereference the data > extract all epoches > reject bad trials
%   run ICA: runica > exclude artifacts with ADjust
%   Baseline correction

% Stage II: extract conditions

% If you do not konw which function you should use, 
% please run one of your data, and type EEG.history.
% created by M.-Y. Wang
% 05-09-2017

%% Stage1-preICA
clear all
clc
cd ('F:\EEG\face-random\Raw_data')

for subi = 101:121;
    EEG = pop_biosig([num2str(subi),'-random.bdf'], 'channels',1:64,'ref',47);
    EEG.setname = [num2str(subi),'_preICA'];
    EEG = pop_chanedit(EEG, 'lookup','D:\\Program Files\\MATLAB\\R2014a\\matlabtoolbox\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
    EEG =pop_eegfiltnew(EEG, 0.5,45,13518,0,[],0);
%     EEG = pop_eegfiltnew(EEG,1,30,6760,0,[],0);
%     EEG = pop_eegfiltnew(EEG, [],0.5,13518,1,[],0);
%     EEG = pop_eegfiltnew(EEG, [],1,6760,1,[],0);
%     EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:64] ,'computepower',1,'linefreqs',[50 100] ,'normSpectrum',0,...
%         'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',4);
    EEG = pop_resample( EEG, 250);
    EEG = pop_reref( EEG, []);
    EEG = pop_epoch( EEG, {'1' '2' '3' '4'}, [-0.7 2.5], 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-200 0]);
    EEG.setname = [num2str(subi),'_PreICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), '_PreICA.set'],'filepath','F:\EEG\face-random\Preprocessing\PreICA');
end
%% Check the data and disgard bad trials.


%% Stage1-RunICA
clear all
clc
cd ('F:\EEG\face-random\Preprocessing\PreICA')

for subi=101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_PreICA.set'],'filepath','F:\EEG\face-random\Preprocessing\PreICA\');
    EEG = pop_runica(EEG,'icatype','fastica','approach','symm');
    EEG.setname = [num2str(subi),'_runICA'];
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), '_runICA.set'],'filepath','F:\EEG\face-random\Preprocessing\ICA');
end
% Stage1-ICA exclude (refer)
clear all
clc
cd ('F:\EEG\face-random\Preprocessing\ICA')

for subi=101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_runICA.set'],'filepath','F:\EEG\face-random\Preprocessing\ICA\');
    EEG = interface_ADJ (EEG, [num2str(subi),'.txt']);
end
%% Stage2-select the conditions
clear all
clc
cd ('F:\EEG\face-random\Preprocessing\ICA')

data(1).name = '_Neutral.set';data(2).name = '_Happy.set';data(3).name = '_N2H.set';data(4).name = '_H2N.set';
data(1).file = 'Condition1_Neutral'; data(2).file = 'Condition2_Happy';data(3).file = 'Condition3_N2H';data(4).file = 'Condition4_H2N';
for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA3.set'],'filepath','F:\EEG\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',1,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(1).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(1).name],'filepath',['F:\EEG\face-random\Preprocessing\Conditions\',data(1).file]);
end

for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA3.set'],'filepath','F:\EEG\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',2,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(2).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(2).name],'filepath',['F:\EEG\face-random\Preprocessing\Conditions\',data(2).file]);
end

for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA3.set'],'filepath','F:\EEG\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',3,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(3).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(3).name],'filepath',['F:\EEG\face-random\Preprocessing\Conditions\',data(3).file]);
end

for subi = 101:121;
    EEG = pop_loadset('filename',[num2str(subi),'_pruned with ICA3.set'],'filepath','F:\EEG\face-random\Preprocessing\ICA');
    EEG = pop_selectevent( EEG, 'type',4,'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [num2str(subi),data(4).name];
    EEG.CSD = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
    EEG = pop_saveset( EEG, 'filename',[num2str(subi), data(4).name],'filepath',['F:\EEG\face-random\Preprocessing\Conditions\',data(4).file]);
end
