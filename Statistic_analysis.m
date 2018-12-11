%% ------------------------------------------------------Statistic Analysis
clear all
load('TGC_40+400+100.mat')
load('TGC_40+400+100_ds.mat')

% load('TGC_40+500+100.mat')
% load('TGC_40+500+100_ds.mat')
%% ---------------------------------------------

diff_d = x2y_d-y2x_d;
diff_s = x2y_s-y2x_s;
diff_N = x2y1-y2x1;
diff_H = x2y2-y2x2;
diff_N2H = x2y3-y2x3;
diff_H2N = x2y4-y2x4;
x1 = mean (diff_d(:,8:12),2);
x2 = mean (diff_s(:,8:12),2);

[H,P,CI,STATS] = ttest(x1,x2,0.05);

H1 = mean (diff_H (:,6:9),2);
H2 = mean (diff_H (:,10:13),2); 
N2H1 = mean (diff_N2H (:,6:9),2);
N2H2 = mean ( diff_N2H (:,10:13),2);
H2N1 = mean (diff_H2N (:,6:9),2);
H2N2 = mean (diff_H2N (:,10:13),2);
