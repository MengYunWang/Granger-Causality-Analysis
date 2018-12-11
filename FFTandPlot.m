%% Transform the time domain into frequency domain
%   and then plot the results

% Created by M.-Y. Wang
% 25-10-2017

%% FFT
clear all
clc
%--------------------------------------------------------------------Condition 1
data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral');
%define parameters
time2ft = dsearchn (EEG.times',[300;2100]);
data_n = size (EEG.CSD(:,time2ft(1):time2ft(2),:),2);
data_hz = linspace (0,EEG.srate/2,floor(data_n/2+1));
%initialization 
Neutral_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name)); %channels*data_points*sub_numbs
Neutral_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Happy_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
Happy_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
N2H_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
N2H_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));
H2N_amp = zeros (EEG.nbchan,length(data_hz),length(data1_name));
H2N_pow = zeros (EEG.nbchan,length(data_hz),length(data1_name));

for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
    data1_ft = fft (EEG.CSD(:,time2ft(1):time2ft(2),:),[],2);
    data1_apt1 = abs (data1_ft)./data_n;
    data1_apt2 = data1_apt1 (:,1:data_n/2+1,:);
    data1_apt2(:,2:end-1,:) = 2.*data1_apt2(:,2:end-1,:);
    data1_pow = data1_apt2.^2;
    Neutral_amp(:,:,ii) = squeeze (mean(data1_apt2,3));
    Neutral_pow(:,:,ii) = squeeze (mean(data1_pow,3));
end
%--------------------------------------------------------------------Condition 2
data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\');
    data2_ft = fft (EEG.CSD(:,time2ft(1):time2ft(2),:),[],2);
    data2_apt1 = abs (data2_ft)./data_n;
    data2_apt2 = data2_apt1 (:,1:data_n/2+1,:);
    data2_apt2(:,2:end-1,:) = 2.*data2_apt2(:,2:end-1,:);
    data2_pow = data2_apt2.^2;
    Happy_amp(:,:,ii) = squeeze (mean(data2_apt2,3));
    Happy_pow(:,:,ii) = squeeze (mean(data2_pow,3));
end
%--------------------------------------------------------------------Condition 3
data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
    data3_ft = fft (EEG.CSD(:,time2ft(1):time2ft(2),:),[],2);
    data3_apt1 = abs (data3_ft)./data_n;
    data3_apt2 = data3_apt1 (:,1:data_n/2+1,:);
    data3_apt2(:,2:end-1,:) = 2.*data3_apt2(:,2:end-1,:);
    data3_pow = data3_apt2.^2;
    N2H_amp(:,:,ii) = squeeze (mean(data3_apt2,3));
    N2H_pow(:,:,ii) = squeeze (mean(data3_pow,3));
end
%--------------------------------------------------------------------Condition 4
data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
    data4_ft = fft (EEG.CSD(:,time2ft(1):time2ft(2),:),[],2);
    data4_apt1 = abs (data4_ft)./data_n;
    data4_apt2 = data4_apt1 (:,1:data_n/2+1,:);
    data4_apt2(:,2:end-1,:) = 2.*data4_apt2(:,2:end-1,:);
    data4_pow = data4_apt2.^2;
    H2N_amp(:,:,ii) = squeeze (mean(data4_apt2,3));
    H2N_pow(:,:,ii) = squeeze (mean(data4_pow,3));
end
save fft_amp_pow Neutral_amp Neutral_pow Happy_amp  Happy_pow N2H_amp N2H_pow H2N_amp H2N_pow data_hz EEG

%% Compute the SNR and Z-score
    Neutral_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));  Neutral_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    Happy_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));    Happy_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    N2H_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));      N2H_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    H2N_SNR = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));      H2N_Z = zeros (EEG.nbchan,length(data_hz),size (Neutral_amp,3));
    for mm = 1:size (Neutral_amp,3);
        for ii = 11:length (data_hz)-11;
            neutral_temp = zeros (64,18);
            neutral_temp(:,1:9) = squeeze (Neutral_amp(:,ii-10:ii-2,mm)); neutral_temp (:,10:18) = squeeze (Neutral_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (neutral_temp,[],2); for jj = 1:64; neutral_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (neutral_temp,[],2); for jj = 1:64; neutral_temp (jj,indx2(jj)) = NaN; end
            Neutral_SNR (:,ii,mm) = squeeze(Neutral_amp(:,ii,mm)) ./ nanmean(neutral_temp,2);
            Neutral_Z (:,ii,mm)= (squeeze(Neutral_amp(:,ii,mm)) - nanmean(neutral_temp,2)) ./ nanstd(neutral_temp,[],2);
            
            happy_temp = zeros (64,18);
            happy_temp(:,1:9) = squeeze (Happy_amp(:,ii-10:ii-2,mm)); happy_temp (:,10:18) = squeeze (Happy_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (happy_temp,[],2); for jj = 1:64; happy_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (happy_temp,[],2); for jj = 1:64; happy_temp (jj,indx2(jj)) = NaN; end
            Happy_SNR (:,ii,mm) = squeeze (Happy_amp(:,ii,mm)) ./ nanmean(happy_temp,2);
            Happy_Z (:,ii,mm) = (squeeze (Happy_amp(:,ii,mm)) - nanmean(happy_temp,2)) ./ nanstd(happy_temp,[],2);
            
            N2H_temp = zeros (64,18);
            N2H_temp(:,1:9) = squeeze (N2H_amp(:,ii-10:ii-2,mm)); N2H_temp (:,10:18) = squeeze (N2H_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (N2H_temp,[],2);   for jj = 1:64; N2H_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (N2H_temp,[],2);   for jj = 1:64; N2H_temp (jj,indx2(jj)) = NaN; end;
            N2H_SNR (:,ii,mm) = squeeze (N2H_amp(:,ii,mm)) ./ nanmean(N2H_temp,2);
            N2H_Z (:,ii,mm) = (squeeze (N2H_amp(:,ii,mm)) - nanmean(N2H_temp,2)) ./ nanstd(N2H_temp,[],2);
            
            H2N_temp = zeros (64,18);
            H2N_temp(:,1:9) = squeeze (H2N_amp(:,ii-10:ii-2,mm)); H2N_temp (:,10:18) = squeeze (H2N_amp(:,ii+2:ii+10,mm));
            [~,indx1] = max (H2N_temp,[],2); for jj = 1:64; H2N_temp (jj,indx1(jj)) = NaN; end
            [~,indx2] = max (H2N_temp,[],2); for jj = 1:64; H2N_temp (jj,indx2(jj)) = NaN; end
            H2N_SNR (:,ii,mm) = squeeze (H2N_amp(:,ii,mm)) ./ nanmean(H2N_temp,2);
            H2N_Z (:,ii,mm) = (squeeze (H2N_amp(:,ii,mm)) - nanmean(H2N_temp,2)) ./ nanstd(H2N_temp,[],2);
        end
    end
save fft_SNRZ  Neutral_SNR  Neutral_Z  Happy_SNR Happy_Z  N2H_SNR N2H_Z H2N_SNR H2N_Z data_hz EEG

%% Plot amplitude
clear all
clc

load fft_SNRZ
figure, clf, 
set (gcf,'color','w')
for chani = 1:16;
    subplot (4, 4,chani)
    plot (data_hz, squeeze (mean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (mean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
        set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.5 3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure, clf
set (gcf,'color','w')
for chani = 17:32;
    subplot (4, 4,chani-16)
    plot (data_hz, squeeze (mean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (mean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.5 3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure, clf
set (gcf,'color','w')
for chani = 33:48;
    subplot (4, 4,chani-32)
    plot (data_hz, squeeze (mean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (mean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.5 3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N

figure, clf
set (gcf,'color','w')
for chani = 49:64;
    subplot (4, 4,chani-48)
    plot (data_hz, squeeze (mean(Neutral_SNR(chani,:,:),3)),'--k','linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(Happy_SNR(chani,:,:),3)),'-.','color',[.5,.5,.5],'linewidth',2.5);
    hold on
    plot (data_hz, squeeze (mean(N2H_SNR(chani,:,:),3)),'-b','linewidth',2.5');
    hold on
    plot (data_hz, squeeze (mean(H2N_SNR(chani,:,:),3)),'-r','linewidth',2.5);
    set (gca,'xlim',[5 35],'xtick',5:5:35,'ylim',[1.5 3])
    title ([num2str(chani), EEG.chanlocs(chani).labels ])
end
legend Neutral Happy N2H H2N
%% plot brian topology
temptdata = ( Neutral_SNR + Happy_SNR  + N2H_SNR + H2N_SNR )./4;
data2plot = squeeze (mean (temptdata,3));

temptdata_s = ( Neutral_SNR + Happy_SNR )./2;
data2plot_s = squeeze (mean (temptdata_s,3));
temptdata_d = ( N2H_SNR + H2N_SNR )./2;
data2plot_d = squeeze (mean (temptdata_d,3));

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
% topoplot (data2plot (:,19),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% topoplot (data2plot_d (:,19),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
topoplot (data2plot_s (:,19),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure, clf
set (gcf,'color','w')
% topoplot (data2plot (:,37),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
% topoplot (data2plot_d (:,37),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
topoplot (data2plot_s (:,37),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')
%% plot 10 and 20 brian topology under different conditions
%---------------------------------------10 hz under different conditions
figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (Neutral_SNR (:,19,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (Happy_SNR (:,19,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (N2H_SNR (:,19,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (H2N_SNR (:,19,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')
%---------------------------------------20 hz under different conditions
figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (Neutral_SNR (:,37,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (Happy_SNR (:,37,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (N2H_SNR (:,37,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

figure , clf
set (gcf,'color','w')
%clo = othercolor ('Paired7');
%topoplot (data2plot (:,21),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',200,'colormap',clo);
topoplot (squeeze (mean (H2N_SNR (:,37,:),3)),EEG.chanlocs,'style','map','electrodes','off','conv','on','gridscale',300);
set(gca,'clim',[0 2.5])
% colorbar([0.85 0.2 0.03 0.5],'ytick',0:0.5:2.5,'YTickLabel',{'0','','1','','2',''},'fontsize',18,'fontweight','bold','fontname','arial black')

%% Plot SNR of Occipital across all conditions and subjects
figure, clf
 SNR_all = zeros (4,226);
 SNR_all(1,:) = squeeze (mean(mean(Neutral_SNR([26:30,63,64],:,:),3),1));
 SNR_all(2,:) = squeeze (mean(mean(Happy_SNR([26:30,63,64],:,:),3),1));
 SNR_all(3,:)= squeeze (mean(mean(N2H_SNR([26:30,63,64],:,:),3),1));
 SNR_all(4,:) = squeeze (mean(mean(H2N_SNR([26:30,63,64],:,:),3),1));
 
 set (gcf,'color','w')
 plot (data_hz,mean(SNR_all),'-k','linewidth',3);
 set (gca,'xlim',[5 35],'xtick',0:10:40,'ylim',[0.75 2.5],'ytick',0:1:4,'linewidth',3)
 set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
 title ('Occipital','FontSize',28,'fontweight','bold','fontname','arial black')
%  xlabel ('frequency (Hz)','FontSize',28,'fontweight','bold','fontname','arial black')
 ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')

%% Plot SNR of Occipital under different conditions
figure, clf
set (gcf,'color','w')
plot (data_hz, squeeze (mean(mean(Neutral_SNR([26:30,63,64],:,:),3),1)),'-k','linewidth',2.5);
hold on
plot (data_hz, squeeze (mean(mean(Happy_SNR([26:30,63,64],:,:),3),1)),'-','color',[.6,.6,.6],'linewidth',2.5);
hold on
plot (data_hz, squeeze (mean(mean(N2H_SNR([26:30,63,64],:,:),3),1)),'-b','linewidth',2.5);
hold on
plot (data_hz, squeeze (mean(mean(H2N_SNR([26:30,63,64],:,:),3),1)),'-r','linewidth',2.5);
 set (gca,'xlim',[9.5 10.5],'xtick',0:10:40,'ylim',[1.6 2],'ytick',0:1:3,'linewidth',2.5,'box','off')
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Occipital','FontSize',16,'fontweight','bold','fontname','arial black')
% xlabel ('frequency (Hz)','FontSize',16,'fontweight','bold','fontname','arial black')
% ylabel (['Amplitude(','\mu','V)'],'FontSize',16,'fontweight','bold','fontname','arial black')
% legend Neutral Happy N2H H2N
 %% Plot SNR of Frontal across all conditions and subjects
figure, clf
 SNR_all = zeros (4,226);
 SNR_all(1,:) = squeeze (mean(mean(Neutral_SNR([1,33,34],:,:),3),1));
 SNR_all(2,:) = squeeze (mean(mean(Happy_SNR([1,33,34],:,:),3),1));
 SNR_all(3,:)= squeeze (mean(mean(N2H_SNR([1,33,34],:,:),3),1));
 SNR_all(4,:) = squeeze (mean(mean(H2N_SNR([1,33,34],:,:),3),1));
 
 set (gcf,'color','w')
 plot (data_hz,mean(SNR_all),'-k','linewidth',3);
 set (gca,'xlim',[5 35],'xtick',0:10:30,'ylim',[0.77 1.7],'ytick',0:0.5:1.5,'linewidth',3)
 set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
 title ('Frontal','FontSize',28,'fontweight','bold','fontname','arial black')
%  xlabel ('frequency (Hz)','FontSize',28,'fontweight','bold','fontname','arial black')
 ylabel ('SNR','FontSize',28,'fontweight','bold','fontname','arial black')
 %% Plot SNR of frontal under different conditions
figure, clf
set (gcf,'color','w')
plot (data_hz, squeeze (mean(mean(Neutral_SNR([1,33,34],:,:),3),1)),'-k','linewidth',2.5);
hold on
plot (data_hz, squeeze (mean(mean(Happy_SNR([1,33,34],:,:),3),1)),'-','color',[.6,.6,.6],'linewidth',2.5);
hold on
plot (data_hz, squeeze (mean(mean(N2H_SNR([1,33,34],:,:),3),1)),'-b','linewidth',2.5);
hold on
plot (data_hz, squeeze (mean(mean(H2N_SNR([1,33,34],:,:),3),1)),'-r','linewidth',2.5);
 set (gca,'xlim',[5 25],'xtick',4:2:26,'ylim',[0.75 1.5],'ytick',0:.25:1.5,'linewidth',2.5)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
% title ('Frontal','FontSize',16,'fontweight','bold','fontname','arial black')
% xlabel ('frequency (Hz)','FontSize',16,'fontweight','bold','fontname','arial black')
% ylabel (['Amplitude(','\mu','V)'],'FontSize',16,'fontweight','bold','fontname','arial black')
% legend Neutral Happy N2H H2N