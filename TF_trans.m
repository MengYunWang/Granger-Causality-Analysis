%% Transform the data into time-frequency domain
% define frequency parameter> define other wavelet parameters> initialize
% output TF data>  main function

% Created by M.-Y. Wang
% 12-10-2017
clear all
clc
%% -------------------------------------Initialize the parameter-------------

data1_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\*.set');
data2_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\*.set');
data3_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\*.set');
data4_name = dir ('F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\*.set');
EEG = pop_loadset('filename',data1_name(1).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral');

% frequency parameters
min_freq =  5; % define the lowest freq
max_freq = 30; % define the highest freq
num_frex = 26; % define the interval between two frex
frex = linspace(min_freq,max_freq,num_frex);

time2save = -1000:20:3000; % down sample the results
time2saveindx = dsearchn (EEG.times',time2save'); % search the index

% other wavelet parameters
range_cycles = [4 10]; % define the cycles; can use the fixed number 3 or 6
s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex)  ./ (2*pi*frex); 
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;
nWave = length(wavtime);

%% --------------------------------------------------- Condition1-Neutral ------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral
% initialize output time-frequency data

tf1_amp = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name)); % freq * chan * time * subs
tf1_power = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name));
ITPC1 = zeros (length(frex),EEG.nbchan,length(time2save),length(data1_name));

for ii = 1:length(data1_name);
    EEG = pop_loadset('filename',data1_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition1_Neutral\');
 
    % FFT parameters
    nData = EEG.pnts*EEG.trials;
    nConv = nWave + nData - 1;
    alldata = reshape (EEG.data,EEG.nbchan,[]);   
    dataX   = fft(alldata,nConv,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv,2);
        as = as(:,half_wave+1:end-half_wave);
        as = reshape (as,EEG.nbchan,EEG.pnts,EEG.trials);
        tf1_amp (fi,:,:,ii) = mean (abs(as(:,time2saveindx,:)),3);
        tf1_power (fi,:,:,ii) = mean (abs(as(:,time2saveindx,:)).^2,3);
        ITPC1 (fi,:,:,ii) = abs(mean(exp(1i*angle(as(:,time2saveindx,:))),3));
    end
end

%% --------------------------------------------------- Condition2-Happy ------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy

tf2_amp = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name));
tf2_power = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name));
ITPC2 = zeros (length(frex),EEG.nbchan,length(time2save),length(data1_name));

for ii = 1:length(data2_name);
    EEG = pop_loadset('filename',data2_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition2_Happy\');
 
    % FFT parameters
    nData = EEG.pnts*EEG.trials;
    nConv = nWave + nData - 1;
    alldata = reshape (EEG.data,EEG.nbchan,[]);   
    dataX   = fft(alldata,nConv,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv,2);
        as = as(:,half_wave+1:end-half_wave);
        as = reshape (as,EEG.nbchan,EEG.pnts,EEG.trials);
        tf2_amp (fi,:,:,ii) = mean (abs(as(:,time2saveindx,:)),3);
        tf2_power (fi,:,:,ii) = mean (abs(as(:,time2saveindx,:)).^2,3);
        ITPC2 (fi,:,:,ii) = abs(mean(exp(1i*angle(as(:,time2saveindx,:))),3));
    end
end

%% --------------------------------------------------- Condition3-N2H ------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H

tf3_amp = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name));
tf3_power = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name));
ITPC3 = zeros (length(frex),EEG.nbchan,length(time2save),length(data1_name));

for ii = 1:length(data3_name);
    EEG = pop_loadset('filename',data3_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition3_N2H\');
 
    % FFT parameters
    nData = EEG.pnts*EEG.trials;
    nConv = nWave + nData - 1;
    alldata3 = reshape (EEG.data,EEG.nbchan,[]);   
    data3X   = fft(alldata3,nConv,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX./max(waveletX);
        as3 = ifft(bsxfun(@times, data3X, waveletX),nConv,2);
        as3 = as3(:,half_wave+1:end-half_wave);
        as3 = reshape (as3,EEG.nbchan,EEG.pnts,EEG.trials);
        tf3_amp (fi,:,:,ii) = mean (abs(as3(:,time2saveindx,:)),3);
        tf3_power (fi,:,:,ii) = mean (abs(as3(:,time2saveindx,:)).^2,3);
        ITPC3 (fi,:,:,ii) = abs(mean(exp(1i*angle(as3(:,time2saveindx,:))),3));
    end
end

%% --------------------------------------------------- Condition4-Neutral ------------------
cd F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N

tf4_amp = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name));
tf4_power = zeros(length(frex),EEG.nbchan,length(time2save),length(data1_name));
ITPC4 = zeros (length(frex),EEG.nbchan,length(time2save),length(data1_name));

for ii = 1:length(data4_name);
    EEG = pop_loadset('filename',data4_name(ii).name,'filepath','F:\EEG\face-random\Preprocessing\Conditions\Condition4_H2N\');
 
    % FFT parameters
    nData = EEG.pnts*EEG.trials;
    nConv = nWave + nData - 1;
    alldata = reshape (EEG.data,EEG.nbchan,[]);   
    dataX   = fft(alldata,nConv,2);
    
    for fi=1:length(frex);
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX./max(waveletX);
        as = ifft(bsxfun(@times, dataX, waveletX),nConv,2);
        as = as(:,half_wave+1:end-half_wave);
        as = reshape (as,EEG.nbchan,EEG.pnts,EEG.trials);
        tf4_amp (fi,:,:,ii) = mean (abs(as(:,time2saveindx,:)),3);
        tf4_power (fi,:,:,ii) = mean (abs(as(:,time2saveindx,:)).^2,3);
        ITPC4 (fi,:,:,ii) = abs(mean(exp(1i*angle(as(:,time2saveindx,:))),3));
    end
end

%% ---------------------------------------- Compute the baseline corrected TF power and bc ITPC
baseline_time = [-700, -600];
baseindex = dsearchn (time2save',baseline_time');

baseline_power =squeeze( median (cat(3,squeeze (median(tf1_power(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(tf2_power(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(tf3_power(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(tf4_power(:,:,baseindex(1):baseindex(2),:),3))),3));

baseline_itpc =  squeeze( median (cat(3,squeeze (median(ITPC1(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(ITPC2(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(ITPC3(:,:,baseindex(1):baseindex(2),:),3)),...
    squeeze (median(ITPC4(:,:,baseindex(1):baseindex(2),:),3))),3));

% basline_itpc =  (squeeze(mean(mean(ITPC1(:,:,baseindex(1):baseindex(2),:),3),4))+...
%     squeeze(mean(mean(ITPC2(:,:,baseindex(1):baseindex(2),:),3),4))+...
%     squeeze(mean(mean(ITPC3(:,:,baseindex(1):baseindex(2),:),3),4))+...
%     squeeze(mean(mean(ITPC4(:,:,baseindex(1):baseindex(2),:),3),4)))./4;


        tf1_dB = 10.*log10((tf1_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));
        ITPC1_bc = ITPC1 - repmat(baseline_itpc,1,1,length(time2save),length(data1_name));

        tf2_dB = 10.*log10((tf2_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));
        ITPC2_bc = ITPC2 - repmat(baseline_itpc,1,1,length(time2save),length(data1_name));
        
        tf3_dB = 10.*log10((tf3_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));
        ITPC3_bc = ITPC3 - repmat(baseline_itpc,1,1,length(time2save),length(data1_name));
        
        tf4_dB = 10.*log10((tf4_power)./repmat(baseline_power,1,1,length(time2save),length(data1_name)));
        ITPC4_bc = ITPC4 - repmat(baseline_itpc,1,1,length(time2save),length(data1_name));


save  ssVEP_TF frex time2save tf1_amp tf1_power ITPC1 tf1_dB ITPC1_bc tf2_amp tf2_power ITPC2 tf2_dB ITPC2_bc...
    tf3_amp tf3_power ITPC3 tf3_dB ITPC3_bc tf4_amp tf4_power ITPC4  tf4_dB ITPC4_bc -v7.3
 
