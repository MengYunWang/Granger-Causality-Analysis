%% Plot GC-time domain

figure,clf
set (gcf,'color','w')

value = spcrv([times2save;squeeze(mean(x2y1))]); % smooth the line----function 'spcrv'
plot(value(1,:),value(2,:),'k','linewidth',3)
hold on
 value = spcrv([times2save;squeeze(mean(x2y2))]);
plot(value(1,:),value(2,:),'color',[0.6 0.6 0.6],'linewidth',3)
hold on;
value = spcrv([times2save;squeeze(mean(x2y3))]);
plot(value(1,:),value(2,:),'r','linewidth',3)
hold on;
value = spcrv([times2save;squeeze(mean(x2y4))]);
plot(value(1,:),value(2,:),'b','linewidth',3)
% plot (times2save,squeeze(mean(x2y1)))
% plot (times2save,squeeze(mean(y2x1)))
set (gca,'xlim',[300 2000],'xtick',0:500:2000,'ylim',[0.075 0.115],'ytick',0.06:0.01:0.2,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
% title('Neutral','FontSize',18,'fontweight','bold','fontname','arial black') 

figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(mean(y2x1))]);
plot(value(1,:),value(2,:),'k','linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(y2x2))]);
plot(value(1,:),value(2,:),'color',[0.6 0.6 0.6],'linewidth',3)
hold on
 value = spcrv([times2save;squeeze(mean(y2x3))]);
plot(value(1,:),value(2,:),'r','linewidth',3)
hold on
 value = spcrv([times2save;squeeze(mean(y2x4))]);
plot(value(1,:),value(2,:),'b','linewidth',3)

set (gca,'xlim',[300 2000],'xtick',0:500:2000,'ylim',[0.035 0.065],'ytick',0.04:0.01:0.06,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
% title('Happy','FontSize',18,'fontweight','bold','fontname','arial black')

% title([ 'Window length: ' num2str(timewin) ' ms, order: ' num2str(order) ' ms' ],'FontSize',20,'fontweight','bold','fontname','arial black')
% xlabel('Time (ms)','FontSize',20,'fontweight','bold','fontname','arial black')
% ylabel('Granger Causality estimate','FontSize',20,'fontweight','bold','fontname','arial black')
% legend({ 'GC: Occipital -> Frontal'; 'GC:Frontal -> Occipital'},'FontSize',16,'fontweight','bold','fontname','arial black')
%% Plot the difference between O2F and F2O
figure,clf
set (gcf,'color','w')

value = spcrv([times2save;squeeze(mean(x2y_d))-squeeze(mean(y2x_d))]); % smooth the line----function 'spcrv'
plot(value(1,:),value(2,:),'k','linewidth',3)
hold on;
value = spcrv([times2save;squeeze(mean(x2y_s))-squeeze(mean(y2x_s))]);
plot(value(1,:),value(2,:),':k','linewidth',3)
set (gca,'xlim',[-100 2100],'xtick',0:500:2000,'ylim',[0.015 0.065],'ytick',0.01:0.01:0.08,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
legend( 'Dynamic','Static','FontSize',16,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(mean(x2y2))-squeeze(mean(y2x2))]);
plot(value(1,:),value(2,:),'color',[0.6 0.6 0.6],'linewidth',3)
hold on;
value = spcrv([times2save;squeeze(mean(x2y3))-squeeze(mean(y2x3))]);
plot(value(1,:),value(2,:),'r','linewidth',3)
hold on;
value = spcrv([times2save;squeeze(mean(x2y4))-squeeze(mean(y2x4))]);
plot(value(1,:),value(2,:),'b','linewidth',3)

set (gca,'xlim',[-100 2100],'xtick',0:500:2000,'ylim',[0.015 0.075],'ytick',0.01:0.01:0.08,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')

% title('O->F - F->O','FontSize',18,'fontweight','bold','fontname','arial black')
% title([ 'Window length: ' num2str(timewin) ' ms, order: ' num2str(order) ' ms' ],'FontSize',20,'fontweight','bold','fontname','arial black')
% xlabel('Time (ms)','FontSize',20,'fontweight','bold','fontname','arial black')
ylabel('Granger Causality','FontSize',20,'fontweight','bold','fontname','arial black')
% legend('Happy', 'N2H','H2N','FontSize',16,'fontweight','bold','fontname','arial black')

%% Plot GC-time domain    Baseline Correction

figure,clf
set (gcf,'color','w')

subplot(2,2,1)
value = spcrv([times2save;squeeze(mean(timebs_granger1(1,:,:),2))']); % smooth the line----function 'spcrv'
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(timebs_granger1(2,:,:),2))']);
plot(value(1,:),value(2,:),'r','linewidth',3)

set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-50 -5],'ytick',-50:10:-10,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Neutral','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,2)
 value = spcrv([times2save;squeeze(mean(timebs_granger2(1,:,:),2))']);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
 value = spcrv([times2save;squeeze(mean(timebs_granger2(2,:,:),2))']);
plot(value(1,:),value(2,:),'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-15 50],'ytick',-20:10:50,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Happy','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,3)
 value = spcrv([times2save;squeeze(mean(timebs_granger3(1,:,:),2))']);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
 value = spcrv([times2save;squeeze(mean(timebs_granger3(2,:,:),2))']);
plot(value(1,:),value(2,:),'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-15 35],'ytick',-20:10:50,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('N2H','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,4)
 value = spcrv([times2save;squeeze(mean(timebs_granger4(1,:,:),2))']);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
 value = spcrv([times2save;squeeze(mean(timebs_granger4(2,:,:),2))']);
plot(value(1,:),value(2,:),'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-5 35],'ytick',-20:10:50,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('H2N','FontSize',18,'fontweight','bold','fontname','arial black')

% title([ 'Window length: ' num2str(timewin) ' ms, order: ' num2str(order) ' ms' ],'FontSize',20,'fontweight','bold','fontname','arial black')
% xlabel('Time (ms)','FontSize',20,'fontweight','bold','fontname','arial black')
% ylabel('Granger Causality estimate','FontSize',20,'fontweight','bold','fontname','arial black')

legend({ 'GC: Occipital -> Frontal'; 'GC:Frontal -> Occipital'},'FontSize',16,'fontweight','bold','fontname','arial black')

%% Plot GC-time domain    Baseline Correction

figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(mean(timebs_granger3(1,:,:),2))']);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(timebs_granger4(1,:,:),2))']);
plot(value(1,:),value(2,:),'r','linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(timebs_granger2(1,:,:),2))']); % smooth the line----function 'spcrv'
plot(value(1,:),value(2,:),'color',[0.5 0.5 0.5],'linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(timebs_granger1(1,:,:),2))']);
plot(value(1,:),value(2,:),':','color',[0.5 0.5 0.5],'linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-45 45],'ytick',-50:10:40,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Occipital -> Frontal','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'N2H'; 'H2N';'Happy';'Neutral'},'FontSize',16,'fontweight','bold','fontname','arial black')


figure,clf
set (gcf,'color','w')
value = spcrv([times2save;squeeze(mean(timebs_granger3(2,:,:),2))']);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(timebs_granger4(2,:,:),2))']);
plot(value(1,:),value(2,:),'r','linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(timebs_granger2(2,:,:),2))']); % smooth the line----function 'spcrv'
plot(value(1,:),value(2,:),'color',[0.5 0.5 0.5],'linewidth',3)
hold on
value = spcrv([times2save;squeeze(mean(timebs_granger1(2,:,:),2))']);
plot(value(1,:),value(2,:),':','color',[0.5 0.5 0.5],'linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-45 50],'ytick',-50:10:50,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Frontal -> Occipital','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'N2H'; 'H2N';'Happy';'Neutral'},'FontSize',16,'fontweight','bold','fontname','arial black')
%% Plot GC-time domain  bs and (O2F SUB F2O)

figure,clf
set (gcf,'color','w')
plot(times2save,timebs_Dynamic,'k','linewidth',3); hold on
plot(times2save,timebs_Static,':k','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-15 22],'ytick',-20:10:20,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Dynamic vs. Static','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'Dynamic'; 'Static'},'FontSize',16,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
plot(times2save,timebs_diff3,'b','linewidth',3);hold on
plot(times2save,timebs_Static,':k','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-15 22],'ytick',-20:10:20,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('N2H vs. Static','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'N2H'; 'Static'},'FontSize',16,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')

plot(times2save,timebs_diff4,'r','linewidth',3); hold on
plot(times2save,timebs_Static,':k','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-15 22],'ytick',-20:10:20,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('H2N vs. Static','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'H2N'; 'Static'},'FontSize',16,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
plot(times2save,timebs_diff3,'b','linewidth',3)
hold on
plot(times2save,timebs_diff4,'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-15 35],'ytick',-20:10:30,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('N2H vs. H2N','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'N2H'; 'H2N'},'FontSize',16,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
plot(times2save,timebs_diff2,'color',[0.5 0.5 0.5],'linewidth',3); hold on
plot(times2save,timebs_diff1,':','color',[0.5 0.5 0.5], 'linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',0:500:2000,'ylim',[-35 35],'ytick',-30:10:30,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Happy vs. Neutral','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'Happy'; 'Neutral'},'FontSize',16,'fontweight','bold','fontname','arial black')

%% Plot GC- frequency

figure,clf
set (gcf,'color','w')

subplot(2,2,1)
 value = spcrv([frequencies;squeeze(mean(freq_granger1(1,:,:),3))],3);
plot(value(1,:),value(2,:),'b','linewidth',3)
hold on
 value = spcrv([frequencies;squeeze(mean(freq_granger1(2,:,:),3))],3);
plot(value(1,:),value(2,:),'r','linewidth',3)
% plot(frequencies,squeeze(mean(freq_granger1(2,:,:),3)),'r','linewidth',3)
set (gca,'xlim',[5 40],'xtick',0:5:40,'ylim',[0.002 0.045],'ytick',0:0.01:0.05);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
title('Neutral','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,2)
plot(frequencies,squeeze(mean(freq_granger2(1,:,:),3)),'b','linewidth',3)
hold on
plot(frequencies,squeeze(mean(freq_granger2(2,:,:),3)),'r','linewidth',3)
set (gca,'xlim',[5 40],'xtick',0:5:40,'ylim',[0.002 0.045],'ytick',0:0.01:0.05);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
title('Happy','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,3)
plot(frequencies,squeeze(mean(freq_granger3(1,:,:),3)),'b','linewidth',3)
hold on
plot(frequencies,squeeze(mean(freq_granger3(2,:,:),3)),'r','linewidth',3)
set (gca,'xlim',[5 40],'xtick',0:5:40,'ylim',[0.002 0.045],'ytick',0:0.01:0.05);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
title('N2H','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,4)
plot(frequencies,squeeze(mean(freq_granger4(1,:,:),3)),'b','linewidth',3)
hold on
plot(frequencies,squeeze(mean(freq_granger4(2,:,:),3)),'r','linewidth',3)
set (gca,'xlim',[5 40],'xtick',0:5:40,'ylim',[0.002 0.045],'ytick',0:0.01:0.05);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
title('H2N','FontSize',18,'fontweight','bold','fontname','arial black')

% title([ 'Window length: ' num2str(timewin) ' ms, order: ' num2str(order) ' ms' ],'FontSize',20,'fontweight','bold','fontname','arial black')
% xlabel('Time (ms)','FontSize',20,'fontweight','bold','fontname','arial black')
% ylabel('Granger Causality estimate','FontSize',20,'fontweight','bold','fontname','arial black')

legend({ 'GC: Occipital -> Frontal'; 'GC: Frontal -> Occipital'},'FontSize',16,'fontweight','bold','fontname','arial black')

%% Plot GC- frequency difference (O2F SUB F2O)


% figure,clf; % Spline the line
% set (gcf,'color','w')
% value1 = spcrv([frequencies;diff1],3);
% value2 = spcrv([frequencies;diff2],3);
% value3 = spcrv([frequencies;diff3],3);
% value4 = spcrv([frequencies;diff4],3);
% plot(value1(1,:),value1(2,:),'color',[0.5,0.5,0.5],'linewidth',3)
% hold on
% plot(value2(1,:),value2(2,:),'k','linewidth',3);
% hold on
% plot(value3(1,:),value3(2,:),'b','linewidth',3);
% hold on
% plot(value4(1,:),value4(2,:),'r','linewidth',3)
% set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)

% figure,clf; % The original value
% set (gcf,'color','w')
% plot(frequencies,diff1,'k','linewidth',3);
% hold on
% plot(frequencies,diff2,'color',[0.5 0.5 0.5],'linewidth',3);
% hold on
% plot(frequencies,diff3,'b','linewidth',3);
% hold on
% plot(frequencies,diff4,'r','linewidth',3);
% set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
%--------------------------------------------dynamic vs static
figure,clf; % The original value
set (gcf,'color','w')
plot(frequencies,Dynamic,'-k','linewidth',3); 
hold on 
plot(frequencies,Static,':k','linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','ylim',[-0.005 0.005],'ytick',0:0.001:0.006,'linewidth',3)
title('Dynamic vs. Static','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'Dynamic'; 'Static'},'FontSize',16,'fontweight','bold','fontname','arial black')
%--------------------------------------------N2H vs static
figure,clf; % The original value
set (gcf,'color','w')
plot(frequencies,diff3,'-b','linewidth',3);
hold on
plot(frequencies,Static,':k','linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','ylim',[0 0.05],'ytick',0:0.1:0.3,'linewidth',3)
title('N2H vs. Static','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'N2H'; 'Static'},'FontSize',16,'fontweight','bold','fontname','arial black')
%--------------------------------------------H2N vs static
figure,clf; % The original value
set (gcf,'color','w')
plot(frequencies,diff4,'-r','linewidth',3);
hold on
plot(frequencies,Static,':k','linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','ylim',[0 0.05],'ytick',0:0.1:0.3,'linewidth',3)
title('H2N vs. Static','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'H2N'; 'Static'},'FontSize',16,'fontweight','bold','fontname','arial black')
%--------------------------------------------N2H vs H2N
figure,clf; % The original value
set (gcf,'color','w')
plot(frequencies,diff3,'-b','linewidth',3);
hold on
plot(frequencies,diff4,'-r','linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','ylim',[0 0.05],'ytick',0:0.1:0.25,'linewidth',3)
title('N2H vs. H2N','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'N2H'; 'H2N'},'FontSize',16,'fontweight','bold','fontname','arial black')
%--------------------------------------------N vs H
figure,clf; % The original value
set (gcf,'color','w')
plot(frequencies,diff2,'-','color',[0.5 0.5 0.5],'linewidth',3);
hold on
plot(frequencies,diff1,':','color',[0.5 0.5 0.5],'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','ylim',[0 0.05],'ytick',0:0.1:0.25,'linewidth',3)
title('Happy vs. Neutral','FontSize',18,'fontweight','bold','fontname','arial black')
legend({ 'Happy'; 'Neutral'},'FontSize',16,'fontweight','bold','fontname','arial black')

%% Plot GC- time-frequency
% clear all
% clc
% load TF_granger_100-400
figure,clf
set (gcf,'color','w')

subplot(2,4,1)
contourf (time2save,frequencies,squeeze(mean(tf_granger1(1,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])
title('Neutral','FontSize',18,'fontweight','bold','fontname','arial black')
ylabel('Occipital -> Frontal','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,4,5)
contourf (time2save,frequencies,squeeze(mean(tf_granger1(2,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])
ylabel('Frontal -> Occipital','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,4,2)
contourf (time2save,frequencies,squeeze(mean(tf_granger2(1,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])
title('Happy','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,4,6)
contourf (time2save,frequencies,squeeze(mean(tf_granger2(2,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])


subplot(2,4,3)
contourf (time2save,frequencies,squeeze(mean(tf_granger3(1,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])
title('N2H','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,4,7)
contourf (time2save,frequencies,squeeze(mean(tf_granger3(2,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])

subplot(2,4,4)
contourf (time2save,frequencies,squeeze(mean(tf_granger4(1,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])
title('H2N','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,4,8)
contourf (time2save,frequencies,squeeze(mean(tf_granger4(2,:,:,:),4)),50,'linecolor','none')
set (gca,'clim',[0 0.05]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3,'ylim',[5 40])
colorbar ('Fontsize',16,'fontweight','bold','fontname','arial black','linewidth',2)
%% Plot GC- time-frequency difference between occipital to frontal and frontal to occipital

figure,clf
set (gcf,'color','w')
contourf (time2save,frequencies,(tfbs_Dynamic+tfbs_Static)/2,50,'linecolor','none')
set (gca,'clim',[-100 100]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
set (gca,'xlim',[-200 2300],'ylim',[3 35],'ytick',0:5:35)
title('Overall','FontSize',18,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
contourf (time2save,frequencies,tfbs_Dynamic,50,'linecolor','none')
set (gca,'clim',[-100 100]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
set (gca,'xlim',[-200 2300],'ylim',[3 35],'ytick',0:5:35)
title('Dynamic','FontSize',18,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
contourf (time2save,frequencies,tfbs_Static,50,'linecolor','none')
set (gca,'clim',[-100 100]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
set (gca,'xlim',[-200 2300],'ylim',[3 35],'ytick',0:5:35)
title('Static','FontSize',18,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
contourf (time2save,frequencies,tfbs_diff1,50,'linecolor','none')
set (gca,'clim',[-100 100]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
set (gca,'xlim',[-200 2300],'ylim',[3 35],'ytick',0:5:35)
title('Neutral','FontSize',18,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
contourf (time2save,frequencies,tfbs_diff2,50,'linecolor','none')
set (gca,'clim',[-100 100]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
set (gca,'xlim',[-200 2300],'ylim',[3 35],'ytick',0:5:35)
title('Happy','FontSize',18,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
contourf (time2save,frequencies,tfbs_diff3,50,'linecolor','none')
set (gca,'clim',[-100 100]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
set (gca,'xlim',[-200 2300],'ylim',[3 35],'ytick',0:5:35)
title('N2H','FontSize',18,'fontweight','bold','fontname','arial black')

figure,clf
set (gcf,'color','w')
contourf (time2save,frequencies,tfbs_diff4,50,'linecolor','none')
set (gca,'clim',[-100 100]);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black','linewidth',3)
set (gca,'xlim',[-200 2300],'ylim',[3 35],'ytick',0:5:35)
title('H2N','FontSize',18,'fontweight','bold','fontname','arial black')


%% GC- time-frequency at 10 or 20 or 30 HZ
figure (7),clf
set (gcf,'color','w')
frex2plot = 30;
frexindx = dsearchn (frequencies',frex2plot);

subplot(2,2,1)
plot(time2save,squeeze(mean(tf_granger1(1,frexindx,:,:),4)),'b','linewidth',3)
hold on
plot(time2save,squeeze(mean(tf_granger1(2,frexindx,:,:),4)),'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',[-200,0:500:2000,2300],'ylim',[0.2 0.8],'ytick',0:0.1:0.8,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Neutral','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,2)
plot(time2save,squeeze(mean(tf_granger2(1,frexindx,:,:),4)),'b','linewidth',3)
hold on
plot(time2save,squeeze(mean(tf_granger2(2,frexindx,:,:),4)),'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',[-200,0:500:2000,2300],'ylim',[0.2 0.8],'ytick',0:0.1:0.8,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('Happy','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,3)
plot(time2save,squeeze(mean(tf_granger3(1,frexindx,:,:),4)),'b','linewidth',3)
hold on
plot(time2save,squeeze(mean(tf_granger3(2,frexindx,:,:),4)),'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',[-200,0:500:2000,2300],'ylim',[0.2 0.8],'ytick',0:0.1:0.8,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('N2H','FontSize',18,'fontweight','bold','fontname','arial black')

subplot(2,2,4)
plot(time2save,squeeze(mean(tf_granger4(1,frexindx,:,:),4)),'b','linewidth',3)
hold on
plot(time2save,squeeze(mean(tf_granger4(2,frexindx,:,:),4)),'r','linewidth',3)
set (gca,'xlim',[-200 2300],'xtick',[-200,0:500:2000,2300],'ylim',[0.2 0.8],'ytick',0:0.1:0.8,'linewidth',3);
set (gca,'FontSize',16,'fontweight','bold','fontname','arial black')
title('H2N','FontSize',18,'fontweight','bold','fontname','arial black')

% title([ 'Window length: ' num2str(timewin) ' ms, order: ' num2str(order) ' ms' ],'FontSize',20,'fontweight','bold','fontname','arial black')
% xlabel('Time (ms)','FontSize',20,'fontweight','bold','fontname','arial black')
% ylabel('Granger Causality estimate','FontSize',20,'fontweight','bold','fontname','arial black')

legend({ 'GC: Occipital -> Frontal'; 'GC:Frontal -> Occipital'},'FontSize',16,'fontweight','bold','fontname','arial black')

