%% GC-Frequency domain

% 2018-01-23
% Meng-Yun Wang
%%
clear all
clc
%% -----------------------------------------------------------Condition1
data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
%---------------------------------------innitialize parameter
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
% EEG = pop_resample( EEG, 250);
timewin = 500; % Granger prediction parametersin ms
timewin_points = round (timewin/(1000/EEG.srate));% convert parameters to indices

min_freq = 0; % in Hz
max_freq = 40; 
frequencies = linspace(min_freq,max_freq,41);

order_time = 40; %40ms
order_points = round (order_time/(1000/EEG.srate)); % how many data pnts in 100ms

time2save = -400:50:2200;
time2saveidx = dsearchn (EEG.times',time2save');

% initialize
tf_granger1=zeros(2,length(frequencies),length (time2save),length(data1_name));

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    
    for timei = 1:length(time2save);
        tempdata = squeeze (eegdata(:,time2saveidx(timei)-floor(timewin_points/2):time2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
        
   %---------------------------------------Main function pwcausal
        [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);
        tf_granger1(1,:,timei,ii) = Fx2y;
        tf_granger1(2,:,timei,ii) = Fy2x;
    end
end

%% ---------------------------------------------------------Condition 2

data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
% initialize
tf_granger2=zeros(2,length(frequencies),length(time2save),length(data2_name));
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    
    for timei = 1:length(time2save);
        
       tempdata = squeeze (eegdata(:,time2saveidx(timei)-floor(timewin_points/2):time2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
                % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
     %---------------------------------------Main function pwcausal        
        [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);
        tf_granger2(1,:,timei,ii) = Fx2y;
        tf_granger2(2,:,timei,ii) = Fy2x;
    end
end 
%% ---------------------------------------------------------Condition 3

data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
% initialize
tf_granger3=zeros(2,length(frequencies),length(time2save),length(data3_name));
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    for timei=1:length(time2save);
        tempdata = squeeze (eegdata(:,time2saveidx(timei)-floor(timewin_points/2):time2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
     %---------------------------------------Main function pwcausal
        [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);
        tf_granger3(1,:,timei,ii) = Fx2y;
        tf_granger3(2,:,timei,ii) = Fy2x;
    end
end

%% ---------------------------------------------------------Condition 4

data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
% initialize
tf_granger4=zeros(2,length(frequencies),length(time2save),length(data4_name));
for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    for timei = 1:length(time2save);
        tempdata = squeeze (eegdata(:,time2saveidx(timei)-floor(timewin_points/2):time2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
    %---------------------------------------Main function pwcausal
    [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);    
    tf_granger4(1,:,timei,ii) = Fx2y;
    tf_granger4(2,:,timei,ii) = Fy2x;
    end
end
% %% Baseline Correction
% baseline_wind = [-500 -200];
% [~,baseline_idx1] = find (time2save == -500);
% [~,baseline_idx2] = find (time2save == -200);
% 
% mean_granger = (tf_granger1 + tf_granger2 + tf_granger3 + tf_granger4)/4;
% baseline_granger = mean (mean_granger(:,:,baseline_idx1:baseline_idx2,:),3);
% 
% tfbs_granger1 = 100 * (bsxfun(@rdivide, bsxfun(@minus, tf_granger1,baseline_granger),baseline_granger));
% tfbs_granger2 = 100 * (bsxfun(@rdivide, bsxfun(@minus, tf_granger2,baseline_granger),baseline_granger));
% tfbs_granger3 = 100 * (bsxfun(@rdivide, bsxfun(@minus, tf_granger3,baseline_granger),baseline_granger));
% tfbs_granger4 = 100 * (bsxfun(@rdivide, bsxfun(@minus, tf_granger4,baseline_granger),baseline_granger));
% 
% tfbs_diff1 = squeeze(mean(tfbs_granger1(1,:,:,:),4))-squeeze(mean(tfbs_granger1(2,:,:,:),4));
% tfbs_diff2 = squeeze(mean(tfbs_granger2(1,:,:,:),4))-squeeze(mean(tfbs_granger2(2,:,:,:),4));
% tfbs_diff3 = squeeze(mean(tfbs_granger3(1,:,:,:),4))-squeeze(mean(tfbs_granger3(2,:,:,:),4));
% tfbs_diff4 = squeeze(mean(tfbs_granger4(1,:,:,:),4))-squeeze(mean(tfbs_granger4(2,:,:,:),4));
% tfbs_Static= (tfbs_diff1 + tfbs_diff2)/2;
% tfbs_Dynamic = (tfbs_diff3 + tfbs_diff4)/2;

save TF_granger_40-500 tf_granger1 tf_granger2 tf_granger3 tf_granger4 frequencies time2save
%     tfbs_granger1 tfbs_granger2 tfbs_granger3 tfbs_granger4 ...
%     tfbs_diff1 tfbs_diff2 tfbs_diff3 tfbs_diff4 tfbs_Static tfbs_Dynamic... 
    