%% Granger Causality
% 2018-01-23
% Meng-yun Wang
%% ---------------------------------------------------------STATIC FACE
clear all
clc
%% --------------------------------------------------Condition 1--Neutral
data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
%---------------------------------------innitialize parameter
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
% EEG = pop_resample( EEG, 100);
timewin = 400; % Granger prediction parametersin ms
timewin_points = round(timewin/(1000/EEG.srate));% convert parameters to indices

times2save = -400:200:2200; % temporal down-sample results (but not data!)in ms
times2saveidx = dsearchn(EEG.times' ,times2save');% convert requested times to indices

% initialize
bic1 = zeros(length(data1_name),length(times2save),40); % Bayes info criteria (hard-coded to order=40)
data_station1 = zeros (length(data1_name),length(times2save));

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
%     EEG = pop_resample( EEG, 100);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save);
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);% channels * trials
        data_station1 (ii,timei) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,(timewin_points)*EEG.trials);
  % %----------------------------------------Get the order
        % test bic1 for optimal model order at each time point
        for bici=1:size(bic1,3)
            % run model
            [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
            % compute Bayes Information Criteria
            bic1(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
    end
end
%
order = 40;% ms
order_points   = round(order/(1000/EEG.srate));
[x2y1,y2x1] = deal(zeros(length(data1_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3)); % remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save);
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
  % %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x1(ii,timei)=log(Ex/E(1,1));
        x2y1(ii,timei)=log(Ey/E(2,2));
    end
end

%% ---------------------------------------------------------Condition 2

data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
% initializeg
bic2 = zeros(length(data2_name),length(times2save),40); % Bayes info criteria (hard-coded to order=15)
data_station2 = zeros (length(data2_name),length(times2save));
[x2y2,y2x2] = deal(zeros(length(data2_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);
        data_station2 (ii,timei) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
        %----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic2,3)
            % run model
            [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
            % compute Bayes Information Criteria
            bic2(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
        %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x2(ii,timei)=log(Ex/E(1,1)); % Front to Occi
        x2y2(ii,timei)=log(Ey/E(2,2)); % Occi to Front
    end
end

%% ------------------------------------------------------------------Condition 3
data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
% initialize
bic3 = zeros(length(data3_name),length(times2save),40); % Bayes info criteria (hard-coded to order=15)
data_station3 = zeros (length(data3_name),length(times2save));
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
%---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);
        data_station3 (ii,timei) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
%----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic3,3)
            % run model
            [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
            % compute Bayes Information Criteria
            bic3(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
    end
end 

% initialize
[x2y3,y2x3] = deal(zeros(length(data3_name),length(times2save))); % the function deal assigns inputs to all outputs
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
%---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
%---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x3(ii,timei)=log(Ex/E(1,1));
        x2y3(ii,timei)=log(Ey/E(2,2));
    end
end
%% --------------------------------------------------------------------Condition 4
data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
% initialize
bic4 = zeros(length(data4_name),length(times2save),40); % Bayes info criteria (hard-coded to order=15)
data_station4 = zeros (length(data4_name),length(times2save));

for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
 %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,EEG.trials,timewin_points,[]);
        data_station4 (ii,timei) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
%----------------------------------------Get the order
        % test bic2 for optimal model order at each time point
        for bici=1:size(bic4,3)
            % run model
            [~,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
            % compute Bayes Information Criteria
            bic4(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
    end
end

% initialize
[x2y4,y2x4] = deal(zeros(length(data4_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
 %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save)
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
%----------------------------------------main procedure
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);

        % time-domain causal  estimate
        y2x4(ii,timei)=log(Ex/E(1,1));
        x2y4(ii,timei)=log(Ey/E(2,2));
    end
end

% %% Baseline Correction
% baseline_wind = [-500 -200];
% [~,baseline_idx1] = find (times2save == -500);
% [~,baseline_idx2] = find (times2save == -200);
% 
% mean_granger = ( x2y1 + y2x1 + x2y2 + y2x2 + x2y3 + y2x3 + x2y4 + y2x4)/8;
% baseline_granger = mean (mean_granger(:,baseline_idx1:baseline_idx2),2);
% 
% timebs_granger1(1,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, x2y1,baseline_granger),baseline_granger));
% timebs_granger1(2,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, y2x1,baseline_granger),baseline_granger));
% timebs_granger2(1,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, x2y2,baseline_granger),baseline_granger));
% timebs_granger2(2,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, y2x2,baseline_granger),baseline_granger));
% timebs_granger3(1,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, x2y3,baseline_granger),baseline_granger));
% timebs_granger3(2,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, y2x3,baseline_granger),baseline_granger));
% timebs_granger4(1,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, x2y4,baseline_granger),baseline_granger));
% timebs_granger4(2,:,:) = 100 * (bsxfun(@rdivide, bsxfun(@minus, y2x4,baseline_granger),baseline_granger));
% 
% 
% timebs_diff1 = squeeze(mean(timebs_granger1(1,:,:),2))-squeeze(mean(timebs_granger1(2,:,:),2));
% timebs_diff2 = squeeze(mean(timebs_granger2(1,:,:),2))-squeeze(mean(timebs_granger2(2,:,:),2));
% timebs_diff3 = squeeze(mean(timebs_granger3(1,:,:),2))-squeeze(mean(timebs_granger3(2,:,:),2));
% timebs_diff4 = squeeze(mean(timebs_granger4(1,:,:),2))-squeeze(mean(timebs_granger4(2,:,:),2));
% 
% timebs_Static= (timebs_diff1 + timebs_diff2)/2;
% timebs_Dynamic = (timebs_diff3 + timebs_diff4)/2;

save TGC_40+400+200  times2save bic1  bic2 bic3 bic4...
    data_station1 data_station2 data_station3 data_station4...
    x2y1 y2x1 x2y2 y2x2 x2y3 y2x3 x2y4 y2x4
%     timebs_granger1 timebs_granger2 timebs_granger3 timebs_granger4...
%     timebs_diff1 timebs_diff2 timebs_diff3 timebs_diff4...
%     timebs_Static timebs_Dynamic

%%
clear all
clc
%--------------------------------------------------Static
data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
%---------------------------------------innitialize parameter
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
% EEG = pop_resample( EEG, 100);
timewin = 400; % Granger prediction parametersin ms
timewin_points = round(timewin/(1000/EEG.srate));% convert parameters to indices

times2save = -400:200:2200; % temporal down-sample results (but not data!)in ms
times2saveidx = dsearchn(EEG.times' ,times2save');% convert requested times to indices

% initialize
bic_s = zeros(length(data1_name),length(times2save),40); % Bayes info criteria (hard-coded to order=40)
data_station_s = zeros (length(data1_name),length(times2save));

order = 40;% ms
order_points   = round(order/(1000/EEG.srate));
[x2y_s,y2x_s] = deal(zeros(length(data1_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
%     EEG = pop_resample( EEG, 100);
    data2go1 = zeros(2,EEG.pnts,EEG.trials);
    data2go1(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go1(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
        EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\');
    %     EEG = pop_resample( EEG, 250);
    data2go2 = zeros(2,EEG.pnts,EEG.trials);
    data2go2(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go2(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
    if size (data2go1,3) > size (data2go2,3)
       data2go1(:,:,size(data2go2,3)+1:size(data2go1,3)) = [];
    elseif size (data2go1,3) < size (data2go2,3)
       data2go2(:,:,size(data2go1,3)+1:size(data2go2,3)) = [];
    end
    
    data2go = (data2go1 + data2go2)./2;
    
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));% remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save);
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
        
        % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,size (data2go,3),timewin_points,[]);% channels * trials
        data_station_s (ii,timei) = sum (sum(unit_root));
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*size (data2go,3));
  % %----------------------------------------Get the order
        % test bic1 for optimal model order at each time point
        for bici=1:size(bic_s,3)
            % run model
            [~,E] = armorf(tempdata,size(data2go,3),timewin_points,bici);
            % compute Bayes Information Criteria
            bic_s(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
     % %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),size(data2go,3),timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),size(data2go,3),timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,size(data2go,3),timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x_s(ii,timei)=log(Ex/E(1,1));
        x2y_s(ii,timei)=log(Ey/E(2,2));    
    end
end
% ----------------------------------------------------Dynamic
data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
% initialize
bic_d = zeros(length(data1_name),length(times2save),40); % Bayes info criteria (hard-coded to order=40)
data_station_d = zeros (length(data1_name),length(times2save));
[x2y_d,y2x_d] = deal(zeros(length(data1_name),length(times2save))); % the function deal assigns inputs to all outputs

for ii = 1 :length(data1_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
    data2go1 = zeros(2,EEG.pnts,EEG.trials);
    data2go1(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go1(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
        EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
    %     EEG = pop_resample( EEG, 250);
    data2go2 = zeros(2,EEG.pnts,EEG.trials);
    data2go2(1,:,:) = squeeze (mean(EEG.CSD([27:30,64],:,:))); %occipital
    data2go2(2,:,:) = squeeze (mean(EEG.CSD([1,33,34],:,:))); %frontal
    if size (data2go1,3) > size (data2go2,3)
       data2go1(:,:,size(data2go2,3)+1:size(data2go1,3)) = [];
    elseif size (data2go1,3) < size (data2go2,3)
       data2go2(:,:,size(data2go1,3)+1:size(data2go2,3)) = [];
    end
    
    data2go = (data2go1 + data2go2)./2;
    
 % %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3)); % remove ERP from selected electrodes to improve stationarity
    for timei=1:length(times2save);
        % data from all trials in this time window
        tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
        
        % detrend and zscore all data
        for triali=1:size(tempdata)
            tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
            tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
        end
               % check covariance stationarity
        unit_root = cca_check_cov_stat_mtrial(tempdata,size (data2go,3),timewin_points,[]);% channels * trials
        data_station_d (ii,timei) = sum (sum(unit_root));
        
        % reshape tempdata for armorf
        tempdata = reshape(tempdata,2,timewin_points*size(data2go,3));
          % %----------------------------------------Get the order
        % test bic1 for optimal model order at each time point
        for bici=1:size(bic_d,3)
            % run model
            [~,E] = armorf(tempdata,size(data2go,3),timewin_points,bici);
            % compute Bayes Information Criteria
            bic_d(ii,timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
        end
  % %---------------------------------------main procedure-------------------------
        
        % fit AR models (model estimation from bsmart toolbox)
        [Ax,Ex] = armorf(tempdata(1,:),size(data2go,3),timewin_points,order_points);
        [Ay,Ey] = armorf(tempdata(2,:),size(data2go,3),timewin_points,order_points);
        [Axy,E] = armorf(tempdata     ,size(data2go,3),timewin_points,order_points);
        
        % time-domain causal  estimate
        y2x_d(ii,timei)=log(Ex/E(1,1));
        x2y_d(ii,timei)=log(Ey/E(2,2));
    end
end
save TGC_40+400+200_ds times2save bic_s  bic_d data_station_s data_station_d...
    x2y_s y2x_s x2y_d y2x_d
