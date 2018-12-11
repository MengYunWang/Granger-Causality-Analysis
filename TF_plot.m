%% plot the time frequency (power)

% Created by M.-Y. Wang
% 21-10-2017
clear all
clc
load ssVEP_TF

tf_dB (:,:,:,1) = squeeze (mean(tf1_dB(:,:,:,:),4));
tf_dB (:,:,:,2) = squeeze (mean(tf2_dB(:,:,:,:),4));
tf_dB (:,:,:,3) = squeeze (mean(tf3_dB(:,:,:,:),4));
tf_dB (:,:,:,4) = squeeze (mean(tf4_dB(:,:,:,:),4));
%% Image the time frequency of each channel, across all conditions and all subjects

figure, clf,  
set (gcf,'color','w')
for chani=1:16;
subplot (4,4,chani)
contourf(time2save,frex(),squeeze(mean (tf_dB(:,chani,:,:),4)),40,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-1000 3000])

title([num2str(chani),EEG.chanlocs(chani).labels])
end

figure, clf
set (gcf,'color','w')
for chani=17:32;
subplot (4,4,chani-16)
contourf(time2save,frex(),squeeze(mean (tf_dB(:,chani,:,:),4)),'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-1000 2000])

title([num2str(chani),EEG.chanlocs(chani).labels])
end

figure, clf
set (gcf,'color','w')
for chani=33:48;
subplot (4,4,chani-32)
contourf(time2save,frex(),squeeze(mean (tf_dB(:,chani,:,:),4)),40,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-1000 3000])

title([num2str(chani),EEG.chanlocs(chani).labels])
end

figure, clf
set (gcf,'color','w')
for chani=49:64;
subplot (4,4,chani-48)
contourf(time2save,frex(),squeeze(mean (tf_dB(:,chani,:,:),4)),40,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-1000 3000])

title([num2str(chani),EEG.chanlocs(chani).labels])
end

%% Image the time frequency of frontal and occipital channels, across all conditions and all subjects
close all
figure (1), clf,  
set (gcf,'color','w')
contourf(time2save,frex(),squeeze (mean(tf_dB(:,[1,33,34],:),2)),80,'linecolor','none')

set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600])
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Frontal','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')
colorbar ('Fontsize',16,'fontweight','bold','fontname','arial black');

figure (2), clf,  
set (gcf,'color','w')
contourf(time2save,frex(),squeeze(mean(tf_dB(:,[27:30,64],:),2)),80,'linecolor','none')

set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600])
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Occipital','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')
colorbar;

%% Image the time frequency of frontal and occipital channels, four conditions across all subjects
close all

figure (1), clf,  
set (gcf,'color','w')
subplot (2,2,1)
contourf(time2save,frex(),squeeze (mean(tf_dB(:,[1,33,34],:,1),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Neutral','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')

subplot (2,2,2)
contourf(time2save,frex(),squeeze (mean(tf_dB(:,[1,33,34],:,2),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Happy','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')

subplot (2,2,3)
contourf(time2save,frex(),squeeze (mean(tf_dB(:,[1,33,34],:,3),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('N2H','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')

subplot (2,2,4)
contourf(time2save,frex(),squeeze (mean(tf_dB(:,[1,33,34],:,4),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('H2N','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')


figure (2), clf,  
set (gcf,'color','w')
subplot (2,2,1)
contourf(time2save,frex(),squeeze(mean(tf_dB(:,[27:30,64],:,1),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Neutral','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')

subplot (2,2,2)
contourf(time2save,frex(),squeeze(mean(tf_dB(:,[27:30,64],:,2),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('Happy','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')

subplot (2,2,3)
contourf(time2save,frex(),squeeze(mean(tf_dB(:,[27:30,64],:,3),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('N2H','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')

subplot (2,2,4)
contourf(time2save,frex(),squeeze(mean(tf_dB(:,[27:30,64],:,4),2)),80,'linecolor','none')
set(gca,'clim',[-4 4],'ydir','normal','xlim',[-500 2600],'linewidth',3)
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title ('H2N','FontSize',16,'fontweight','bold','fontname','arial black')
xlabel ('Time(ms)','FontSize',16,'fontweight','bold','fontname','arial black')
ylabel ('Hz','FontSize',16,'fontweight','bold','fontname','arial black')
%% plot time course of ssVEP amplitude of each channels, four conditions at 10Hz OR 20Hz
% close all
% figure(11), clf
% set (gcf,'color','w')
% for chani=1:16;
% subplot (4,4,chani)
% plot (time2save,squeeze (mean(tf1_dB([6 ],chani,:,:),4)),'-k','linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf2_dB([6 ],chani,:,:),4)),'-','color',[.6,.6,.6],'linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf3_dB([6 ],chani,:,:),4)),'-b','linewidth',2.5')
% hold on
% plot (time2save,squeeze (mean(tf4_dB([6 ],chani,:,:),4)),'-r','linewidth',2.5)
% 
% set(gca,'ylim',[-2 3],'ydir','normal','xlim',[-200 2500])
% title([num2str(chani),EEG.chanlocs(chani).labels])
% end
% 
% figure(22), clf
% set (gcf,'color','w')
% for chani=17:32;
% subplot (4,4,chani-16)
% plot (time2save,squeeze (mean(tf1_dB([6 ],chani,:,:),4)),'-k','linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf2_dB([6 ],chani,:,:),4)),'-','color',[.6,.6,.6],'linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf3_dB([6 ],chani,:,:),4)),'-b','linewidth',2.5')
% hold on
% plot (time2save,squeeze (mean(tf4_dB([6 ],chani,:,:),4)),'-r','linewidth',2.5)
% 
% set(gca,'ylim',[-2 3],'ydir','normal','xlim',[-200 2500])
% title([num2str(chani),EEG.chanlocs(chani).labels])
% end
% 
% figure(33), clf
% set (gcf,'color','w')
% for chani=33:48;
% subplot (4,4,chani-32)
% plot (time2save,squeeze (mean(tf1_dB([6 ],chani,:,:),4)),'-k','linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf2_dB([6 ],chani,:,:),4)),'-','color',[.6,.6,.6],'linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf3_dB([6 ],chani,:,:),4)),'-b','linewidth',2.5')
% hold on
% plot (time2save,squeeze (mean(tf4_dB([6 ],chani,:,:),4)),'-r','linewidth',2.5)
% 
% set(gca,'ylim',[-2 3],'ydir','normal','xlim',[-200 2500])
% title([num2str(chani),EEG.chanlocs(chani).labels])
% end
% 
% figure(44), clf
% set (gcf,'color','w')
% for chani=49:64;
% subplot (4,4,chani-48)
% plot (time2save,squeeze (mean(tf1_dB([6 ],chani,:,:),4)),'-k','linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf2_dB([6 ],chani,:,:),4)),'-','color',[.6,.6,.6],'linewidth',2.5)
% hold on
% plot (time2save,squeeze (mean(tf3_dB([6 ],chani,:,:),4)),'-b','linewidth',2.5')
% hold on
% plot (time2save,squeeze (mean(tf4_dB([6 ],chani,:,:),4)),'-r','linewidth',2.5)
% 
% set(gca,'ylim',[-2 3],'ydir','normal','xlim',[-200 2500])
% title([num2str(chani),EEG.chanlocs(chani).labels])
% end
% %% plot time course of ssVEP dB of Frontal and Occipital sites, four conditions
% close all
% %------------------------------------------------Occipital 10Hz
% figure(1), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf_dB(6,[27:30,64],:,:),4),2)),'-k','linewidth',2.5)
% hold on
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% 
% set(gca,'ylim',[-3.5 0.5],'ydir','normal','xlim',[-200 2500],'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Occipital-10Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
% %----------------------------------------------Occipital 20Hz
% figure(2), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf_dB(16,[27:30,64],:,:),4),2)),'-k','linewidth',2.5)
% hold on
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% 
% set(gca,'ylim',[-1.5 2],'ydir','normal','xlim',[-200 2500], 'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Occipital-20Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
% %-----------------------------------------------Frontal 10Hz
% figure(3), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf_dB(6,[1,33,34],:,:),4),2)),'-k','linewidth',2.5)
% hold on
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% 
% set(gca,'ylim',[-2 1],'ydir','normal','xlim',[-200 2500],'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Frontal-10Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
% %----------------------------------------------Frontal 20Hz
% figure(4), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf_dB(16,[1,33,34],:,:),4),2)),'-k','linewidth',2.5)
% hold on
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% 
% set(gca,'ylim',[-1 0.8],'ydir','normal','xlim',[-200 2500], 'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Frontal-20Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
% %% plot time course of ssVEP amplitude of Frontal and Occipital sites, four conditions
% close all
% %------------------------------------------------Occipital 10Hz
% figure(1), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf1_dB(6,[27:30,64],:,:),4),2)),'-k','linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf2_dB(6,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf3_dB(6,[27:30,64],:,:),4),2)),'-b','linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf4_dB(6,[27:30,64],:,:),4),2)),'-r','linewidth',3)
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([-100,-100],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([200,200],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([1200,1200],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([1800,1800],[-5,12],'color','k','linewidth',1,'linestyle','--');
% set(gca,'ylim',[-3.5 0.5],'ydir','normal','xlim',[-200 2500],'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Occipital-10Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
% %----------------------------------------------Occipital 20Hz
% figure(2), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf1_dB(16,[27:30,64],:,:),4),2)),'-k','linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf2_dB(16,[27:30,64],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf3_dB(16,[27:30,64],:,:),4),2)),'-b','linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf4_dB(16,[27:30,64],:,:),4),2)),'-r','linewidth',3)
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([-100,-100],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([200,200],[-5,12],'color','k','linewidth',1,'linestyle','--');
% 
% set(gca,'ylim',[-1.5 2],'ydir','normal','xlim',[-200 2500], 'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Occipital-20Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
% %-----------------------------------------------Frontal 10Hz
% figure(3), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf1_dB(6,[1,33,34],:,:),4),2)),'-k','linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf2_dB(6,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf3_dB(6,[1,33,34],:,:),4),2)),'-b','linewidth',3')
% hold on
% plot (time2save,squeeze (mean(mean(tf4_dB(6,[1,33,34],:,:),4),2)),'-r','linewidth',3)
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([-100,-100],[-5,12],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([200,200],[-5,12],'color','k','linewidth',1,'linestyle','--');
% 
% set(gca,'ylim',[-2 1],'ydir','normal','xlim',[-200 2500],'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Frontal-10Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
% %----------------------------------------------Frontal 20Hz
% figure(4), clf
% set (gcf,'color','w')
% 
% plot (time2save,squeeze (mean(mean(tf1_dB(16,[1,33,34],:,:),4),2)),'-k','linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf2_dB(16,[1,33,34],:,:),4),2)),'-','color',[.6,.6,.6],'linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf3_dB(16,[1,33,34],:,:),4),2)),'-b','linewidth',3)
% hold on
% plot (time2save,squeeze (mean(mean(tf4_dB(16,[1,33,34],:,:),4),2)),'-r','linewidth',3)
% line ([-500,2500],[0,0],'color','k','linewidth',1,'linestyle','--');
% hold on
% line ([0,0],[-5,12],'color','k','linewidth',1,'linestyle','--');
% 
% set(gca,'ylim',[-1 0.8],'ydir','normal','xlim',[-200 2500], 'linewidth',3)
% set (gca,'FontSize',28,'fontweight','bold','fontname','arial black')
% title ('Frontal-20Hz','FontSize',28,'fontweight','bold','fontname','arial black')
% xlabel ('Time(ms)','FontSize',28,'fontweight','bold','fontname','arial black')
% ylabel ('dB','FontSize',28,'fontweight','bold','fontname','arial black')
