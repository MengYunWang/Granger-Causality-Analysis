%% -----------------------sig2mTspect_nv.m------------------------
function [S,f]= sig2mTspect_nv(X,fs,fRes);
%Usage: [S, f] = sig2mTspect_nv(X,fs,fRes);
%This function computes auto- & cross- spectra by using multitapers
%Inputs: X is multichannel data (a 3D-matrix in the form of time x trial x channel)
%               fs = sampling rate in Hz
%               nv stands for 'not vectorized program'
%               fRes = desired (lower) frequency resolution (e.g. 1), achieved with zero-padding 
%              default frequency resolution is fs/datalength
%Outputs: S = 3D matrix: m by m spectral matrix at each frequency point of f
%Note:     One can change nw (half the number of tapers) below and see the effect
%Written by M. Dhamala, UF, August 2006.
%Revised by M. Dhamala, GSU, June, 2017

[N,Ntr,m] = size(X); % N = timepoints, Ntr = trials, m = channels 

if nargin<3|fRes>fs/N,
       npad = 0;  fRes = fs/N;  
end 
fRes0 = fs/N;
if (nargin==3) & (fRes<=fs/N), 
    npad = round((fs/fRes-N)/2);  %These many zeros will be padded on each side of the data
end
f = fs*(0:fix((N+2*npad)/2))/(N+2*npad);% upto Nyquist-f

nw = 2; % number of tapers = 2*nw .......good nw are 1.5, 2, 3, 4, 5, 6, or 7...
[tapers,v] = dpss(N+2*npad, nw);

S = zeros(m,m,N+2*npad);

for itrial = 1: Ntr,
     for ii = 1: m, Xft(:,:,ii) = mtfft(squeeze(X(:,itrial,ii)),tapers,fs,npad); end
     for ii = 1:m,
           for jj = 1:m,
                s(ii,jj,:) = squeeze(mean(Xft(:,:,ii).*conj(Xft(:,:,jj)),2));
                %averaging over tapers 
          end
     end
     S = S + s; 
end
S = S/Ntr; %averaging over trials
S = S(:,:,1:fix(end/2)+1)/fs;%half part of two-sided spectra
S = S*fRes0/fRes;%multiplying by a factor to adjust distribution of power from zero-padding
end
%-------------------------
function xf  = mtfft(data,tapers,fs,npad);
%Usage: xf = mtfft(data,tapers,fs,npad);
%Written by M. Dhamala (August 2006)
x0 = zeros(npad,size(data,2));
data = cat(1,x0,data); data= cat(1,data,x0);
data = data(:,ones(1,size(tapers,2)));
data = data.*tapers;
xf = fft(data,[],1);
end

