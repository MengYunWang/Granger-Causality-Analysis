%% Transform the data into time-frequency domain
% define frequency parameter> define other wavelet parameters> initialize
% output TF data>  main function

% Created by M.-Y. Wang
% 12-10-2017

%%
clear all
clc
% -------------------------------------Initialize the parameter-------------

data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral');

% frequency parameters
frex =  [10,20]; % define the lowest freq

% other wavelet parameters
range_cycles = [4 6]; % define the cycles; can use the fixed number 3 or 6
s = range_cycles./(2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;
nWave = length(wavtime);

% --------------------------------------------------- Condition1-Neutral ------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral
% initialize output time-frequency data

Neutral_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name)); % freq * chan * time * subs
Neutral_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    ERP_data = squeeze(mean(EEG.CSD,3));
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv);
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
        Neutral_tfamp (fi,:,:,ii) = abs(as);
        Neutral_power (fi,:,:,ii) = abs(as).^2;
    end
end

% --------------------------------------------------- Condition2-Happy ------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy

Happy_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
Happy_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\');
    ERP_data = squeeze(mean(EEG.CSD,3));
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv);    
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
        Happy_tfamp (fi,:,:,ii) = abs(as);
        Happy_power (fi,:,:,ii) = abs(as).^2;
    end
end

% --------------------------------------------------- Condition3-N2H ------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H

N2H_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
N2H_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
    ERP_data = squeeze(mean(EEG.CSD,3)); 
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv); 
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
        N2H_tfamp (fi,:,:,ii) = abs(as);
        N2H_power (fi,:,:,ii) = abs(as).^2;
    end
end

% --------------------------------------------------- Condition4-H2N------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N

H2N_tfamp = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
H2N_power = zeros(length(frex),EEG.nbchan,EEG.pnts,length(data1_name));
for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
    ERP_data = squeeze(mean(EEG.CSD,3)); 
    % FFT parameters
    nConv = nWave + EEG.pnts - 1;
    nConv_pow2 = 2^nextpow2(nConv); 
    alldata = ERP_data;   
    dataX   = fft(alldata,nConv_pow2,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv_pow2);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv_pow2,2);
        as = as (:,1:nConv);
        as = as(:,half_wave+1:end-half_wave);
        H2N_tfamp (fi,:,:,ii) = abs(as);
        H2N_power (fi,:,:,ii) = abs(as).^2;
    end
end
time2save = EEG.times;

% ---------------------------------------- Compute the baseline corrected TF power
baseline_time = [-500, -300];
baseindex = dsearchn (time2save',baseline_time');

baseline_power =squeeze( median (cat(3,squeeze (median(Neutral_power(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(Happy_power(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(N2H_power(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(H2N_power(:,:,baseindex(1):baseindex(2),:),3))),3));

        Neutral_dB = 10.*log10((Neutral_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));

        Happy_dB = 10.*log10((Happy_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));
        
        N2H_dB = 10.*log10((N2H_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));
        
        H2N_dB = 10.*log10((H2N_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));

save  ssVEP_TF_10+20 frex time2save Neutral_dB Happy_dB N2H_dB  H2N_dB Neutral_tfamp  Happy_tfamp N2H_tfamp H2N_tfamp -v7.3

%% --------------------------Plot an averaged ssVEP amplitude of Occipital electrodes
% 27:30 O1 IZ OZ POZ 64 O2
figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_tfamp(1,[27:30,64],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_tfamp(1,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_tfamp(1,[27:30,64],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_tfamp(1,[27:30,64],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([1700,1700],[0,2],'color','k','linewidth',1,'linestyle','--');
%     set (gca,'xlim',[-50 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')

    figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_tfamp(2,[27:30,64],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_tfamp(2,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_tfamp(2,[27:30,64],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_tfamp(2,[27:30,64],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on

%     set (gca,'xlim',[-50 2500],'ylim',[0.15,1],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
%     legend Neutral  Happy N2H H2N

%% --------------------------Plot an averaged ssVEP (dB) of Occipital electrodes
% 27:30 O1 IZ OZ POZ 64 O2
figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_dB(1,[27:30,64],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_dB(1,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_dB(1,[27:30,64],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_dB(1,[27:30,64],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-10,20],'color','k','linewidth',1,'linestyle','--');
    
    set (gca,'xlim',[-200 2500],'ylim',[-4,10],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')

    figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_dB(2,[27:30,64],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_dB(2,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_dB(2,[27:30,64],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_dB(2,[27:30,64],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on
    line ([0,0],[-10,20],'color','k','linewidth',1,'linestyle','--');
    
    set (gca,'xlim',[-200 2500],'ylim',[-4,10],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
%     legend Neutral  Happy N2H H2N
%% --------------------------Plot an averaged ssVEP of Occipital electrodes across all conditions
% 27:30 O1 IZ OZ POZ 64 O2
figure, clf
 tf1 = zeros (4,4000); tf2 = zeros (4,4000);
 tf1(1,:) = squeeze (mean(mean(Neutral_tfamp(1,[27:30,64],:,:),4),2)); tf2(1,:) = squeeze (mean(mean(Neutral_tfamp(2,[27:30,64],:,:),4),2));
 tf1(2,:) = squeeze (mean(mean(Happy_tfamp(1,[27:30,64],:,:),4),2)); tf2(2,:) = squeeze (mean(mean(Happy_tfamp(2,[27:30,64],:,:),4),2));
 tf1(3,:)= squeeze (mean(mean(N2H_tfamp(1,[27:30,64],:,:),4),2));  tf2(3,:)= squeeze (mean(mean(N2H_tfamp(2,[27:30,64],:,:),4),2));
 tf1(4,:) = squeeze (mean(mean(H2N_tfamp(1,[27:30,64],:,:),4),2)); tf2(4,:) = squeeze (mean(mean(H2N_tfamp(2,[27:30,64],:,:),4),2));
set (gcf,'color','w')
    plot (time2save, mean (tf1),'-k','linewidth',3);
    hold on;
    plot (time2save, mean (tf2),':k','linewidth',3)

%     set (gca,'xlim',[-200 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
 legend 10Hz 20Hz
 %% --------------------------Plot an averaged ssVEP (dB) of Occipital electrodes across all conditions
% 27:30 O1 IZ OZ POZ 64 O2
    figure, clf
    tf1 = zeros (4,800); tf2 = zeros (4,800);
    tf1(1,:) = squeeze (mean(mean(Neutral_dB(1,[27:30,64],:,:),4),2)); tf2(1,:) = squeeze (mean(mean(Neutral_dB(2,[27:30,64],:,:),4),2));
    tf1(2,:) = squeeze (mean(mean(Happy_dB(1,[27:30,64],:,:),4),2)); tf2(2,:) = squeeze (mean(mean(Happy_dB(2,[27:30,64],:,:),4),2));
    tf1(3,:)= squeeze (mean(mean(N2H_dB(1,[27:30,64],:,:),4),2));  tf2(3,:)= squeeze (mean(mean(N2H_dB(2,[27:30,64],:,:),4),2));
    tf1(4,:) = squeeze (mean(mean(H2N_dB(1,[27:30,64],:,:),4),2)); tf2(4,:) = squeeze (mean(mean(H2N_dB(2,[27:30,64],:,:),4),2));
    set (gcf,'color','w')
    plot (time2save, mean (tf1),'-k','linewidth',3);
    hold on;
    plot (time2save, mean (tf2),':k','linewidth',3)

    set (gca,'xlim',[-700 2500],'ylim',[-7,10],'ytick',-5:5:10,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
%     xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
%     legend 10Hz 20Hz
    
%% --------------------------Plot an averaged ssVEP amplitude of Frontal electrodes
% 1 FP1 33 FPZ 34 FP2
figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_tfamp(1,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_tfamp(1,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_tfamp(1,[1,33,34],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_tfamp(1,[1,33,34],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
%     hold on
%     line ([1700,1700],[0,2],'color','k','linewidth',1,'linestyle','--');
%     set (gca,'xlim',[-50 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')

    figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_tfamp(2,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_tfamp(2,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_tfamp(2,[1,33,34],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_tfamp(2,[1,33,34],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on

%     set (gca,'xlim',[-50 2500],'ylim',[0.15,1],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
%     legend Neutral  Happy N2H H2N

%% --------------------------Plot an averaged ssVEP (dB) of Frontal electrodes
% 1 FP1 33 FPZ 34 FP2
figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_dB(1,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_dB(1,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_dB(1,[1,33,34],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_dB(1,[1,33,34],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
%     hold on
%     line ([1700,1700],[0,2],'color','k','linewidth',1,'linestyle','--');
%     set (gca,'xlim',[0 2500],'ylim',[-4,16],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')

    figure, clf
set (gcf,'color','w')
    plot (time2save, squeeze (mean(mean(Neutral_dB(2,[1,33,34],:,:),4),2)),'-k','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(Happy_dB(2,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(N2H_dB(2,[1,33,34],:,:),4),2)),'-b','linewidth',2.5);
    hold on
    plot (time2save, squeeze (mean(mean(H2N_dB(2,[1,33,34],:,:),4),2)),'-r','linewidth',2.5);
    hold on
    line ([-200,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
    hold on

%     set (gca,'xlim',[0 2500],'ylim',[-4,16],'ytick',-5:5:15,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
%     legend Neutral  Happy N2H H2N
%% --------------------------Plot an averaged ssVEP of Frontal electrodes across all conditions
% 1 FP1 33 FPZ 34 FP2
figure, clf
 tf1 = zeros (4,4000); tf2 = zeros (4,4000);
 tf1(1,:) = squeeze (mean(mean(Neutral_tfamp(1,[1,33,34],:,:),4),2)); tf2(1,:) = squeeze (mean(mean(Neutral_tfamp(2,[1,33,34],:,:),4),2));
 tf1(2,:) = squeeze (mean(mean(Happy_tfamp(1,[1,33,34],:,:),4),2)); tf2(2,:) = squeeze (mean(mean(Happy_tfamp(2,[1,33,34],:,:),4),2));
 tf1(3,:)= squeeze (mean(mean(N2H_tfamp(1,[1,33,34],:,:),4),2));  tf2(3,:)= squeeze (mean(mean(N2H_tfamp(2,[1,33,34],:,:),4),2));
 tf1(4,:) = squeeze (mean(mean(H2N_tfamp(1,[1,33,34],:,:),4),2)); tf2(4,:) = squeeze (mean(mean(H2N_tfamp(2,[1,33,34],:,:),4),2));
set (gcf,'color','w')
    plot (time2save, mean (tf1),'-k','linewidth',3);
    hold on;
    plot (time2save, mean (tf2),':k','linewidth',3)

%     set (gca,'xlim',[-200 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel (['\mu','V'],'FontSize',28,'fontweight','bold','fontname','arial black')
 legend 10Hz 20Hz
 %% --------------------------Plot an averaged ssVEP (dB) of Frontal electrodes across all conditions
% 1 FP1 33 FPZ 34 FP2
    figure, clf
    tf1 = zeros (4,4000); tf2 = zeros (4,4000);
    tf1(1,:) = squeeze (mean(mean(Neutral_dB(1,[1,33,34],:,:),4),2)); tf2(1,:) = squeeze (mean(mean(Neutral_dB(2,[1,33,34],:,:),4),2));
    tf1(2,:) = squeeze (mean(mean(Happy_dB(1,[1,33,34],:,:),4),2)); tf2(2,:) = squeeze (mean(mean(Happy_dB(2,[1,33,34],:,:),4),2));
    tf1(3,:)= squeeze (mean(mean(N2H_dB(1,[1,33,34],:,:),4),2));  tf2(3,:)= squeeze (mean(mean(N2H_dB(2,[1,33,34],:,:),4),2));
    tf1(4,:) = squeeze (mean(mean(H2N_dB(1,[1,33,34],:,:),4),2)); tf2(4,:) = squeeze (mean(mean(H2N_dB(2,[1,33,34],:,:),4),2));
    set (gcf,'color','w')
    plot (time2save, mean (tf1),'-k','linewidth',3);
    hold on;
    plot (time2save, mean (tf2),':k','linewidth',3)

%     set (gca,'xlim',[-200 2500],'ylim',[0.15,1.3],'ytick',0:0.2:1.4,'linewidth',3);
    set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
    xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
    ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
    legend 10Hz 20Hz    