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
timeperiod = [300 2100]; % Granger prediction parametersin ms
times_idx = round ((timeperiod+1000)/(1000/EEG.srate)+1);% convert parameters to indices
timewin_points = times_idx(2)-times_idx(1)+1;

min_freq = 0; % in Hz
max_freq = 40; 
frequencies = linspace(min_freq,max_freq,41);

order = 40;% ms
order_points   = round(order/(1000/EEG.srate));
% initialize
freq_granger1=zeros(2,length(frequencies),length(data1_name));

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
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
    
    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
    
    %---------------------------------------Main function pwcausal
    [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);    
    freq_granger1(1,:,ii) = Fx2y; % Occi to Front
    freq_granger1(2,:,ii) = Fy2x; % Front to Occi
    
%     % fit AR models
%     [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
%     [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
%     [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
%     
%     % code below is adapted from bsmart toolbox function pwcausal.m
%     % corrected covariance
%     eyx = E(2,2) - E(1,2)^2/E(1,1); 
%     exy = E(1,1) - E(2,1)^2/E(2,2);
%     N = size(E,1);
%     
%     for fi=1:length(frequencies)
%         
%         % transfer matrix (note the similarity to Fourier transform)
%         H = eye(N);
%         for m = 1:order_points
%             H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*frequencies(fi)/EEG.srate);
%         end
%         
%         Hi = inv(H);
%         S  = H\E*Hi'/EEG.srate;
%         
%         % granger prediction per frequency
%         freq_granger1(1,fi,ii) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/EEG.srate) ); % from channel 1 to channel 2
%         freq_granger1(2,fi,ii) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/EEG.srate) ); % from channel 2 to channel 1
%     end
end 

%% ---------------------------------------------------------Condition 2

data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
% initialize
freq_granger2=zeros(2,length(frequencies),length(data2_name));
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    tempdata = squeeze(eegdata(:,times_idx(1):times_idx(2),:));
    
    % detrend and zscore all data
    for triali=1:size(tempdata)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
    end
    
    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
    %---------------------------------------Main function pwcausal
    [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);    
    freq_granger2(1,:,ii) = Fx2y;
    freq_granger2(2,:,ii) = Fy2x;
    
end 
%% ---------------------------------------------------------Condition 3

data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
% initialize
freq_granger3=zeros(2,length(frequencies),length(data3_name));
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    tempdata = squeeze(eegdata(:,times_idx(1):times_idx(2),:));
    
    % detrend and zscore all data
    for triali=1:size(tempdata)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
    end
    
    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
    
    %---------------------------------------Main function pwcausal
    [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);    
    freq_granger3(1,:,ii) = Fx2y;
    freq_granger3(2,:,ii) = Fy2x;
    
end

%% ---------------------------------------------------------Condition 4

data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
% initialize
freq_granger4=zeros(2,length(frequencies),length(data4_name));
for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    % (note that the ERP-subtracted data are used)
    tempdata = squeeze(eegdata(:,times_idx(1):times_idx(2),:));
    
    % detrend and zscore all data
    for triali=1:size(tempdata)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
    end
    
    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
    
    %---------------------------------------Main function pwcausal
    [~,~,Fx2y,Fy2x] = pwcausal (tempdata,EEG.trials,timewin_points,order_points,EEG.srate,frequencies);    
    freq_granger4(1,:,ii) = Fx2y;
    freq_granger4(2,:,ii) = Fy2x;

end
%%
diff1 = squeeze(mean(freq_granger1(1,:,:),3))-squeeze(mean(freq_granger1(2,:,:),3));
diff2 = squeeze(mean(freq_granger2(1,:,:),3))-squeeze(mean(freq_granger2(2,:,:),3));
diff3 = squeeze(mean(freq_granger3(1,:,:),3))-squeeze(mean(freq_granger3(2,:,:),3));
diff4 = squeeze(mean(freq_granger4(1,:,:),3))-squeeze(mean(freq_granger4(2,:,:),3));
Static= (diff1 + diff2)/2;
Dynamic = (diff3 + diff4)/2;

save freq_ganger_300+2300 freq_granger1 freq_granger2 freq_granger3 freq_granger4 diff1 diff2 diff3 diff4 Static Dynamic frequencies 