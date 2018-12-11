%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)
clear all
clc
%%
regmode   = [];
morder    = 25;
acmaxlags = []; 

f1 = cell(1,21);
data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
for ii = 1:length(data1_name)
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
%     EEG = pop_resample( EEG, 250);
    data2go = zeros(2,EEG.pnts,EEG.trials);
    data2go(1,:,:) = squeeze (mean(EEG.data([27:30,64],:,:))); %occipital
    data2go(2,:,:) = squeeze (mean(EEG.data([1,33,34],:,:))); %frontal
    %---------------------------------------Preprocessing
    eegdata = bsxfun(@minus,data2go(:,:,:),mean(data2go(:,:,:),3));
    tempdata = eegdata(:,1001:3101,:);
      % detrend and zscore all data
    
    for triali=1:size(tempdata)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
    end  
    [A,SIG] = tsdata_to_var(tempdata,morder,regmode);
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    % Calculate spectral pairwise-conditional causalities at given frequency
    % resolution - again, this only requires the autocovariance sequence.
    
    f = autocov_to_spwcgc(G,[]);
    f1{1,ii} = f;
    % Check for failed spectral GC calculation
    assert(~isbad(f,false),'spectral GC calculation failed');

end
%% Plot spectral causal graph.
fs  = EEG.srate;
xlims = [1 50];
for ii= 1:21
    P = f1{1,ii};
    ylims = [min(P(:)) 1.1*max(P(:))];
    h = size(P,3);
    fres = h-1;
    lam = sfreqs(fres,fs)';
    [~,b] = min (abs(xlims(1)-lam));
    [~,d] = min (abs(xlims(2)-lam));
    
    
    figure (ii),clf;
    set (gcf,'color','w');
    subplot (1,2,1)
    plot(lam,squeeze(P(2,1,:)));
    xlim([lam(b) lam(d)]); ylim(ylims); ylabel('Occipital -> Frontal');    
        
    subplot (1,2,2)
    plot(lam,squeeze(P(1,2,:)));
    xlim([lam(b) lam(d)]); ylim(ylims); ylabel('Frontal -> Occipital');
end






