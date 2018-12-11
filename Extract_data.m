%% Extract data 
% Created by M.-Y. Wang
% 25-10-2017
%% 
% clear all
% clc
% %% Extract ERP 
% cd ('F:\EEG\face-random\Preprocessing\Conditions')
% load ssVEP_ERP
% time2test = 0:100:2400;
% time2test_indx = dsearchn (EEG.times',time2test');
% 
% for  timewind = 1:length (time2test);
%     eval (['erp_all_' num2str(timewind) '= zeros (21,8);'])
% end
% 
% for timewind = 1:length(time2test);
%     eval (['erp_all_' num2str(timewind) '(:,1) = squeeze (mean(mean (Neutral_data([27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     eval (['erp_all_' num2str(timewind) '(:,2) = squeeze (mean(mean (Happy_data([27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     eval (['erp_all_' num2str(timewind) '(:,3) = squeeze (mean(mean (N2H_data([27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     eval (['erp_all_' num2str(timewind) '(:,4) = squeeze (mean(mean (H2N_data([27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     eval (['erp_all_' num2str(timewind) '(:,5) = squeeze (mean(mean (Neutral_data([1,33,34],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     eval (['erp_all_' num2str(timewind) '(:,6) = squeeze (mean(mean (Happy_data([1,33,34],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     eval (['erp_all_' num2str(timewind) '(:,7) = squeeze (mean(mean (N2H_data([1,33,34],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     eval (['erp_all_' num2str(timewind) '(:,8) = squeeze (mean(mean (H2N_data([1,33,34],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2)));'])
%     xlswrite ( ['erp_' num2str(timewind)] ,eval(['erp_all_' num2str(timewind)]) );
% end
% %% Extract ERP 
% clear all
% clc
% cd ('F:\EEG\face-random\Preprocessing\Conditions')
% load ssVEP_ERP
% % Emotion (1,[1:21, 64:21*4])= 1; Emotion (1,22:21*3)=2;
% % Dynamic (1,1:21*2)= 1; Dynamic (1,43:21*4) = 2;
% % Subs(1,1:1:21)=1:21; Subs(1,22:21*2)=1:21; Subs(1,43:21*3)=1:21; Subs(1,64:21*4)=1:21;
% time2test = 1100;
% time2test_indx = dsearchn (EEG.times',time2test);
%     erp_all = zeros (21,8);
%     erp_all (:,1) = squeeze (mean(mean (Neutral_data([27:30,64],time2test_indx+1:time2test_indx+1000,:),2)));
%     erp_all (:,2) = squeeze (mean(mean (Happy_data([27:30,64],time2test_indx+1:time2test_indx+1000,:),2)));
%     erp_all (:,3) = squeeze (mean(mean (N2H_data([27:30,64],time2test_indx+1:time2test_indx+1000,:),2)));
%     erp_all (:,4) = squeeze (mean(mean (H2N_data([27:30,64],time2test_indx+1:time2test_indx+1000,:),2)));
%     erp_all (:,5) = squeeze (mean(mean (Neutral_data([1 33 34],time2test_indx+1:time2test_indx+1000,:),2)));
%     erp_all (:,6) = squeeze (mean(mean (Happy_data([1 33 34],time2test_indx+1:time2test_indx+1000,:),2)));
%     erp_all (:,7) = squeeze (mean(mean (N2H_data([1 33 34],time2test_indx+1:time2test_indx+1000,:),2)));
%     erp_all (:,8) = squeeze (mean(mean (H2N_data([1 33 34],time2test_indx+1:time2test_indx+1000,:),2)));

%% Extract Z score of amplitude and power
clear all
clc

load fft_SNRZ

amp_all_10 = zeros (21,8);
amp_all_20 = zeros (21,8);
%--------------------------------10 Hz--------------
amp_all_10 (:,1) = (squeeze(mean (Neutral_Z([27:30,64],19,:))));
amp_all_10 (:,2) = (squeeze(mean (Happy_Z([27:30,64],19,:))));
amp_all_10 (:,3) = (squeeze(mean (N2H_Z([27:30,64],19,:))));
amp_all_10 (:,4) = (squeeze(mean (H2N_Z([27:30,64],19,:))));
amp_all_10 (:,5) = (squeeze(mean (Neutral_Z([1,33,34],19,:))));
amp_all_10 (:,6) = (squeeze(mean (Happy_Z([1,33,34],19,:))));
amp_all_10 (:,7) = (squeeze(mean (N2H_Z([1,33,34],19,:))));
amp_all_10 (:,8) = (squeeze(mean (H2N_Z([1,33,34],19,:))));
xlswrite ('amp_10',amp_all_10); 

%--------------------------------20 Hz--------------
amp_all_20 (:,1) = (squeeze(mean (Neutral_Z([27:30,64],37,:))));
amp_all_20 (:,2) = (squeeze(mean (Happy_Z([27:30,64],37,:))));
amp_all_20 (:,3) = (squeeze(mean (N2H_Z([27:30,64],37,:))));
amp_all_20 (:,4) = (squeeze(mean (H2N_Z([27:30,64],37,:))));
amp_all_20 (:,5) = (squeeze(mean (Neutral_Z([1,33,34],37,:))));
amp_all_20 (:,6) = (squeeze(mean (Happy_Z([1,33,34],37,:))));
amp_all_20 (:,7) = (squeeze(mean (N2H_Z([1,33,34],37,:))));
amp_all_20 (:,8) = (squeeze(mean (H2N_Z([1,33,34],37,:))));
xlswrite ('amp_20',amp_all_20); 
%% Extract time-frequency data 
% cd ('F:\EEG\face-random\Preprocessing\Conditions')
clear all
clc
load ssVEP_TF_10+20
time2test = 0:100:2100;
time2test_indx = dsearchn (time2save',time2test');
%----------------------------------------------
for  timewind = 1:length (time2test);
    eval (['tf_all' num2str(timewind) '= zeros (21,8);'])
end

for timewind = 1:length(time2test);
    eval (['tf_all' num2str(timewind) '(:,1) = squeeze (mean(mean (Neutral_dB(1,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    eval (['tf_all' num2str(timewind) '(:,2) = squeeze (mean(mean (Happy_dB(1,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    eval (['tf_all' num2str(timewind) '(:,3) = squeeze (mean(mean (N2H_dB(1,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    eval (['tf_all' num2str(timewind) '(:,4) = squeeze (mean(mean (H2N_dB(1,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    eval (['tf_all' num2str(timewind) '(:,5) = squeeze (mean(mean (Neutral_dB(2,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    eval (['tf_all' num2str(timewind) '(:,6) = squeeze (mean(mean (Happy_dB(2,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    eval (['tf_all' num2str(timewind) '(:,7) = squeeze (mean(mean (N2H_dB(2,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    eval (['tf_all' num2str(timewind) '(:,8) = squeeze (mean(mean (H2N_dB(2,[27:30,64],time2test_indx(timewind)+1:time2test_indx(timewind)+100,:),2),3));'])
    xlswrite ( ['tf_all_' num2str(timewind)] ,eval(['tf_all' num2str(timewind)]) );
end

%% Extract tf ssvep amplitude 
    clear all
    clc
    cd ('F:\EEG\face-random\Preprocessing\Conditions')
    load ssVEP_TF_10+20
    time2test = [300 1100 2100];
    time2test_indx(1) = dsearchn (time2save',time2test(1));
    time2test_indx(2) = dsearchn (time2save',time2test(2));
    time2test_indx(3) = dsearchn (time2save',time2test(3));

    tf_all = zeros (21,8);
    tf_all (:,1) = squeeze (mean(mean (Neutral_dB(1,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    tf_all (:,2) = squeeze (mean(mean (Happy_dB(1,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    tf_all (:,3) = squeeze (mean(mean (N2H_dB(1,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    tf_all (:,4) = squeeze (mean(mean (H2N_dB(1,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    tf_all (:,5) = squeeze (mean(mean (Neutral_dB(2,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    tf_all (:,6) = squeeze (mean(mean (Happy_dB(2,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    tf_all (:,7) = squeeze (mean(mean (N2H_dB(2,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    tf_all (:,8) = squeeze (mean(mean (H2N_dB(2,[27:30,64],time2test_indx(1):time2test_indx(3),:),2),3));
    xlswrite ( 'tf_10_20_300-2100',tf_all);
    
    tf_all1 = zeros (21,8);
    tf_all1 (:,1) = squeeze (mean(mean (Neutral_dB(1,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    tf_all1 (:,2) = squeeze (mean(mean (Happy_dB(1,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    tf_all1 (:,3) = squeeze (mean(mean (N2H_dB(1,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    tf_all1 (:,4) = squeeze (mean(mean (H2N_dB(1,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    tf_all1 (:,5) = squeeze (mean(mean (Neutral_dB(2,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    tf_all1 (:,6) = squeeze (mean(mean (Happy_dB(2,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    tf_all1 (:,7) = squeeze (mean(mean (N2H_dB(2,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    tf_all1 (:,8) = squeeze (mean(mean (H2N_dB(2,[27:30,64],time2test_indx(1):time2test_indx(2),:),2),3));
    xlswrite ( 'tf_10_20_300-1100',tf_all);
    
    
    
    
    tf_all2 = zeros (21,8);
    tf_all2 (:,1) = squeeze (mean(mean (Neutral_dB(1,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    tf_all2 (:,2) = squeeze (mean(mean (Happy_dB(1,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    tf_all2 (:,3) = squeeze (mean(mean (N2H_dB(1,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    tf_all2 (:,4) = squeeze (mean(mean (H2N_dB(1,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    tf_all2 (:,5) = squeeze (mean(mean (Neutral_dB(2,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    tf_all2 (:,6) = squeeze (mean(mean (Happy_dB(2,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    tf_all2 (:,7) = squeeze (mean(mean (N2H_dB(2,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    tf_all2 (:,8) = squeeze (mean(mean (H2N_dB(2,[27:30,64],time2test_indx(2):time2test_indx(3),:),2),3));
    xlswrite ( 'tf_10_20_1100-2100',tf_all);