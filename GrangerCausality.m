
clear all
clc
%--------------------------------------------------------------------Condition 1
data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral');
timeperiod = [300 2100]; % Granger prediction parametersin ms
times_idx = round ((timeperiod+700)/(1000/EEG.srate)+1);% convert parameters to indices

Neutral_gc = zeros (226,2,2,length(data1_name));

for ii = 1:length(data1_name);
    
    data2go = zeros (2,size (EEG.CSD,2),size(EEG.CSD,3));
    data2go (1,:,:) = mean( EEG.CSD([1 33 34],:,:)); data2go (2,:,:) = mean( EEG.CSD([27:30,64],:,:));
    
%     eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    tempdata = data2go(:,times_idx(1):times_idx(2),:);
%     % detrend and zscore all data
%     for triali=1:size(tempdata)
%         tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
%         tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
%         % At this point with real data, you might want to check for stationarity
%         % and possibly discard or mark data epochs that are non-stationary.
%     end
    data2gc = permute (tempdata,[2,3,1]);
    
    [S,f]= sig2mTspect_nv(data2gc,EEG.srate,1);
    spectra = permute(S,[3 1 2]);Neutral_coh = S2coh(spectra);
    for ichan = 1: size (data2gc,3),
        Neutral_pw(:,ichan) = 2*spectra(:,ichan,ichan);%one-sided power
    end
    [H, Z] = wilson_sf(S,EEG.srate);
    Neutral_gc(:,:,:,ii) = hz2cgcAll(S,H,Z,f,EEG.srate);
end

figure, clf
set (gcf,'color','w')

plot (f, squeeze( mean (Neutral_gc(:,1,2,:),4)),'r')
hold on;
plot (f, squeeze(mean(Neutral_gc(:,2,1,:),4)),'b')
 set (gca,'xlim',[0 50])

%% -------------------------------------------------------------------Condition 2
data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
EEG = pop_loadset('filename',data2_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy');
timeperiod = [300 2100]; % Granger prediction parametersin ms
times_idx = round ((timeperiod+700)/(1000/EEG.srate)+1);% convert parameters to indices

Happy_gc = zeros (226,2,2,length(data2_name));

for ii = 1:length(data2_name);
    
    data2go = zeros (2,size (EEG.CSD,2),size(EEG.CSD,3));
    data2go (1,:,:) = mean( EEG.CSD([1 33 34],:,:)); data2go (2,:,:) = mean( EEG.CSD([27:30,64],:,:));
    
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    tempdata = eegdata(:,times_idx(1):times_idx(2),:);
    % detrend and zscore all data
    for triali=1:size(tempdata)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        % At this point with real data, you might want to check for stationarity
        % and possibly discard or mark data epochs that are non-stationary.
    end
    data2gc = permute (tempdata,[2,3,1]);
    
    [S,f]= sig2mTspect_nv(data2gc,EEG.srate,1);
    spectra = permute(S,[3 1 2]);Happy_coh = S2coh(spectra);
    for ichan = 1: size (data2gc,3),
        Happy_pw(:,ichan) = 2*spectra(:,ichan,ichan);%one-sided power
    end
    [H, Z] = wilson_sf(S,EEG.srate);
    Happy_gc(:,:,:,ii) = hz2cgcAll(S,H,Z,f,EEG.srate);
end

% 
figure, clf
set (gcf,'color','w')

plot (f, squeeze( mean (Happy_gc(:,1,2,:),4)),'r')
hold on;
plot (f, squeeze(mean (Happy_gc(:,2,1,:),4)),'b')
 set (gca,'xlim',[0 50])
 %% -------------------------------------------------------Condition 3
data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
EEG = pop_loadset('filename',data3_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H');
timeperiod = [300 2100]; % Granger prediction parametersin ms
times_idx = round ((timeperiod+700)/(1000/EEG.srate)+1);% convert parameters to indices

N2H_gc = zeros (226,2,2,length(data3_name));

for ii = 1:length(data3_name);
    
    data2go = zeros (2,size (EEG.CSD,2),size(EEG.CSD,3));
    data2go (1,:,:) = mean( EEG.CSD([1 33 34],:,:)); data2go (2,:,:) = mean( EEG.CSD([27:30,64],:,:));
    
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    tempdata = eegdata(:,times_idx(1):times_idx(2),:);
    % detrend and zscore all data
    for triali=1:size(tempdata)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        % At this point with real data, you might want to check for stationarity
        % and possibly discard or mark data epochs that are non-stationary.
    end
    data2gc = permute (tempdata,[2,3,1]);
    
    [S,f]= sig2mTspect_nv(data2gc,EEG.srate,1);
    spectra = permute(S,[3 1 2]);N2H_coh = S2coh(spectra);
    for ichan = 1: size (data2gc,3),
        N2H_pw(:,ichan) = 2*spectra(:,ichan,ichan);%one-sided power
    end
    [H, Z] = wilson_sf(S,EEG.srate);
    N2H_gc(:,:,:,ii) = hz2cgcAll(S,H,Z,f,EEG.srate);
end
% 
figure, clf
set (gcf,'color','w')

plot (f, squeeze( mean (N2H_gc(:,1,2,:),4)),'r')
hold on
plot (f, squeeze( mean (N2H_gc(:,2,1,:),4)),'b')
 set (gca,'xlim',[0 50])
 
 %% -------------------------------------------------------------Condition 4
data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
EEG = pop_loadset('filename',data4_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N');
timeperiod = [300 2100]; % Granger prediction parametersin ms
times_idx = round ((timeperiod+700)/(1000/EEG.srate)+1);% convert parameters to indices

H2N_gc = zeros (226,2,2,length(data4_name));

for ii = 1:length(data4_name);
    
    data2go = zeros (2,size (EEG.CSD,2),size(EEG.CSD,3));
    data2go (1,:,:) = mean( EEG.CSD([1 33 34],:,:)); data2go (2,:,:) = mean( EEG.CSD([27:30,64],:,:));
    
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    tempdata = eegdata(:,times_idx(1):times_idx(2),:);
    % detrend and zscore all data
    for triali=1:size(tempdata)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        % At this point with real data, you might want to check for stationarity
        % and possibly discard or mark data epochs that are non-stationary.
    end
    data2gc = permute (tempdata,[2,3,1]);
    
    [S,f]= sig2mTspect_nv(data2gc,EEG.srate,1);
    spectra = permute(S,[3 1 2]);H2N_coh = S2coh(spectra);
    for ichan = 1: size (data2gc,3),
        H2N_pw(:,ichan) = 2*spectra(:,ichan,ichan);%one-sided power
    end
    [H, Z] = wilson_sf(S,EEG.srate);
    H2N_gc(:,:,:,ii) = hz2cgcAll(S,H,Z,f,EEG.srate);
end
% 
figure, clf
set (gcf,'color','w')

plot (f, squeeze( mean (H2N_gc(:,1,2,:),4)),'r')
hold on
plot (f, squeeze( mean (H2N_gc(:,2,1,:),4)),'b')
 set (gca,'xlim',[0 50])
