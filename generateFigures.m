
function generateFigures();

X = SPexam1; %simulate example 1 of Stokes and Purdon's PNAS
%--- generate Figure 1+ ------------
%fs = 120; [Mn,Sn] = compute_allnpCGCvar3(X,fs); %NonParametric CGC Calculation
fs = 120; [Mn]=computeCGCvar3(X,fs); %NonParametric CGC Calculation with a unified program
%Parametric CGC Calculation with the unified program, model order, p = 3
%Please note that the parametric calculations are better with finer
%frequency resolution, e.g. para.freq = 0:0.001:fs/2; or, resolution from
%the nonparametric method for the given length if the length is sufficiently long
%p = 3;freq = Mn.freq;[Mp,Sp] = compute_allpCGCvar3(X,fs,freq,p); %matrix - partitioning + spectral factorization
para.p = 3; para.freq = Mn.freq; [Mp] = computeCGCvar3(X,fs,para); %spectral factorization-based
displayGCcohPow(Mp,Mn); %plotting Granger causality, coherence and power from both methods
%figure(1); print -depsc -r500 nfig1.eps; print -dtiff -r500 nfig1.tiff;

%------- generate Figure with VAR GC results for p = 6 and p = 20 
Mp3 = Mp; % with p = 3
para.p = 6; para.freq = Mn.freq; [Mp6] = computeCGCvar3(X,fs,para); 
para.p = 20; para.freq = Mn.freq; [Mp20] = computeCGCvar3(X,fs,para); 
freq = para.freq; 
figure; 
plot(freq,Mp3.gc(:,1,2),'b',freq,Mp6.gc(:,1,2),'g--', freq, Mp20.gc(:,1,2),'r--','linewidth',2); 
set(gca,'linewidth',1,'fontsize',18); grid on;
title('VAR-based Estimates of GC: 1\rightarrow{2}|3','fontsize',18,'interpreter','tex');
xlabel('frequency','fontsize',18); ylabel('Granger causality','fontsize',18);
h = legend('p=3','p = 6', 'p = 20');
%figure(3); print -depsc -r500 supfig.eps; print -dtiff -r500 supfig.tiff;

end

%-----------------------------------%-------------------------------------
function displayGCcohPow(Mp,Mn);
%Usage: displayGCcohPow(Mp,Mn);
% M.freq = frequencies, M.gc = Granger causality
% M.coh = coherence, M.pow = power spectral density
% if you don't have M, do the following: 
%         M.freq = f; M.pow = p; M.coh = c; M.gc=gc; M.fs = fs;
% Here, power in db = 10*log10(p*fs/2);
%Written by M. Dhamala

%---------for Mp -------------------------
M = Mp;
f = M.freq; p = M.pow; c = M.coh; g = M.gc; 
fs = M.fs;
p = 10*log10(p*fs/2); % power in Decibel

M = Mn;
f1 = M.freq; p1 = M.pow; c1 = M.coh; g1 = M.gc; 
p1 = 10*log10(p1*fs/2); % power in Decibel

ylimsp = [min(p(:))*2.1 1.1*max(p(:))]; %for power plot
ylimsc = [min(c(:))*2.1 1.1*max(c(:))]; %for coherence plot
%ylimsg = [min(g(:))*2.1 1.1*max(g(:))]; %for Granger causality plot
ylimsg = [-0.5 0.5]; %for this particular plot
n = size(g,2);
cm1 = 'b';
cm2 = 'g';

figure(1); k = 0; 
for i = 1:n
    for j = 1:n
        k = k+1;
        if i ~= j
            subplot(n,n,k);
            plot(f,g(:,i,j),cm1,f1,g1(:,i,j),[cm2 '--'],'LineWidth',2); grid on; 
            axis('square'); 
            xlim([f(1) f(end)]); ylim(ylimsg);
            if k==2, ylim([-0.5 1.1*max(g(:))]); end;
            xlabel('frequency');
            rr = setdiff([1:n],[i j]);
            ylabel(['GC: ' num2str(i) '\rightarrow' num2str(j) '|' num2str(rr)],'interpreter','tex'); 
            %ylabel(sprintf('GC: %s \rightarrow %s',num2str(i),num2str(j)));
            set(gca,'fontsize',15,'linewidth',1);
            
        end
        if  i==j
            subplot(n,n,k); 
            plot(f,p(:,i),cm1,f1,p1(:,i),[cm2 '--'],'LineWidth',2); grid on; 
            xlim([f(1) f(end)]); ylim(ylimsp);
            xlabel('frequency');ylabel(['Power (' num2str(i) ')']);
            %ylabel(['power-' num2str(i) ' (dB)']);
            axis('square'); 
            set(gca,'fontsize',15,'linewidth',1);
       end
    end
end
ax = suptitle('Power and Granger Causality (GC) Spectra (blue: parametric with p = 3, green: nonparametric)'); 
set(ax,'fontsize',18); 
set(gcf,'color','white'); hold on;

figure(2); k = 0; 
for i = 1:n
    for j = 1:n
        k = k+1;
        if i ~= j
            subplot(n,n,k);
            plot(f,c(:,i,j),cm1,f1,c1(:,i,j),[cm2 '--'],'LineWidth',2); grid on; hold on;
            axis('square');
            xlim([f(1) f(end)]);ylim(ylimsc);
            xlabel('frequency');
            %rr = setdiff([1:n],[i j]);
            ylabel(sprintf('coh (%s,%s)',num2str(i),num2str(j)));
            set(gca,'fontsize',15,'linewidth',1);
            
        end
        if  i==j
            sb = subplot(n,n,k); 
            plot(f,p(:,i),cm1,f1,p1(:,i),[cm2 '--'],'LineWidth',2); grid on; hold on;
            xlim([f(1) f(end)]); ylim(ylimsp);
            xlabel('frequency'); ylabel(['Power (' num2str(i) ')']);
            %ylabel(['power-' num2str(i) ' (dB)']);
            axis('square'); 
            set(gca,'fontsize',15,'linewidth',1);
       end
    end
end
ax = suptitle('Power and Coherence Spectra (blue: Parametric)'); set(ax,'fontsize',18);  
set(gcf,'color','white'); hold on;

end
%--------------------------------------------------------------------------
function [X,x] = SPexam1(Nt,Ntr);
% This function implements a 3-node interacting autoregressive dynamical system
%Driving: 1--> 2 --> 3
%There is no instantaneous causality because of completely independent noise rd's;
if nargin<1,
   Nt = 500;
   Ntr = 1000;
end
% system parameters -------------
r1 = 0.9;
r2 = 0.7;
r3 = 0.8;

f1 = 40; f2 = 10; f3 = 50;

dt = 1/120;
th1 = f1*dt*2*pi;
th2 = f2*dt*2*pi;
th3 = f3*dt*2*pi;

A1 =[2*r1*cos(th1) 0 0; -0.356 2*r2*cos(th2) 0; 0 -0.3098 2*r3*cos(th3)];
A2 = [-r1^2 0 0; 0.7136 -r2^2 0; 0 0.5 -r3^2];
A3 = [0 0 0; -0.356 0 0; 0 -0.3098 0];



ntra=3000; % number of transient points
% noise terms with zero mean and unit variances (sigmas)------

sigma1=1;
sigma2=1;
sigma3=1;

%initialize ----


x=[];

for l=1:Ntr,
    rd1=sigma1*randn(Nt+ntra,1);
    rd2=sigma2*randn(Nt+ntra,1);
    rd3=sigma3*randn(Nt+ntra,1);
    wt = [rd1'; rd2'; rd3'];

    x1(1:4) = rand(4,1);
    x2(1:4) = rand(4,1);
    x3(1:4) = rand(4,1);
    
    for k=4:Nt+ntra,

        %Xk = [x1(k); x2(k); x3(k)];

        Xkm1 =[x1(k-1); x2(k-1); x3(k-1)];
        Xkm2 =[x1(k-2); x2(k-2); x3(k-2)];
        Xkm3 =[x1(k-3); x2(k-3); x3(k-3)];

        Xk = A1*Xkm1 + A2*Xkm2 + A3*Xkm3+ wt(:,k);
        x1(k) = Xk(1); x2(k) = Xk(2); x3(k) = Xk(3);
    end
    
    xvar=[x1(ntra+1:end);x2(ntra+1:end);x3(ntra+1:end)];%
    clear x1 x2 x3;
    x=[x,xvar];

end

X = reshape(x',Nt,Ntr,size(x,1));

end

