function [M] = computeCGCvar3(X,fs,para)
%Usage: [M] = computevarCGC(X,fs,para); 
%This routine computes conditional Granger causality (CGC), 
% coherence (coherence) and power (power)
% given the trivariate data (var3) X in the form of time x trial x channel.
% If para.p and para.freq (model order (p) and frequencies (freq)) are supplied, then it
% uses the VAR parametric method. If para is empty or nargin<3, it uses the
% nonparametric approach. 
%-------------------------------------------------------------------------------------------------
%Inputs: X = trivariate data in the form of 3D matrix  for time. trial. channel.
%               fs = data sampling rate in Hz
%               para.freq = frequencies at which GC is computed, eg. para.freq = 0:0.1:fs/2;
%               para.p = model order (e.g. 3), or
%Outputs: M.freq = freq, M.gc = conditional Granger causality (i to j conditional on k), 
%                                                       e.g., 1 to 2 conditional on 3 
%M.coh = coherence, M.pow = power spectral density, M.freq = f, M.fs = fs;             
%-------------------------------------------------------------------------------------------------
%Ref : M. Dhamala, et al. NeuroImage and PRL (2008). Written by M. Dhamala
%Revised by M. Dhamala on March, 2018 
%-------------------------------------------------------------------------------------------------
[Nt, Ntr,Nc] = size(X); %Nt = number of timepoints, Ntr = trials, Nc = channels 

if nargin<3
    para = []; 
end

if nargin<3||isempty(para)==1 %nonparametric approach
    fRes = fs/Nt; 
    [S,freq]= sig2mTspect_nv(X,fs,fRes); 
    spectra = permute(S,[3 1 2]);coherence = S2coh(spectra);
    for ichan = 1: Nc
        power(:,ichan) = 2*spectra(:,ichan,ichan);%one-sided power
    end
    cgc = getCGC(S,fs,freq);
else                 % parametric approach 
   x = reshape(X,Nt*Ntr,Nc); x = x'; % x in the form of channel. (time x trials)
   p = para.p; freq = para.freq;
   [A, Z]=armorf(x,Ntr,Nt,p); %parameters by autoregressive fitting
   [S,H] = AZ2spectra(A,Z,p,freq,fs);
   spectra = permute(S,[3 1 2]);coherence = S2coh(spectra);
   for ichan = 1: Nc
        power(:,ichan) = 2*spectra(:,ichan,ichan);%one-sided power
   end
   cgc = getCGC(S,fs,freq);
end
M.freq = freq; M.pow = power; 
M.coh = coherence; M.gc=cgc; 
M.fs = fs;

end
%----------------------------------sig2mTspect_nv.m------------------------------------------

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
%Written and revised by M. Dhamala


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
%-------------------------------------------------------------------------------------
function xf  = mtfft(data,tapers,fs,npad);
%Usage: xf = mtfft(data,tapers,fs,npad);
%Written by M. Dhamala 
x0 = zeros(npad,size(data,2));
data = cat(1,x0,data); data= cat(1,data,x0);
data = data(:,ones(1,size(tapers,2)));
data = data.*tapers;xf = fft(data,[],1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------
% Functions below : armorf.m, AZ2spectra.m, spectrum.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------
%---------------------------------------------------------------------
function varargout = armorf(x,ntrls,npts,p)
% Script performs AR parameter estimation via LWR method by Morf modified.

%   X is a matrix whose every row is one variable's time series
%   ntrls is the number of realizations, npts is the length of every realization
%   If the time series are stationary long, just let ntrls=1, npts=length(x)
%
%   A = ARMORF(X,NR,NL,ORDER) returns the polynomial coefficients A corresponding to
%   the AR model estimate of matrix X using Morf's method.
%   ORDER is the order of the AR model.
%
%   [A,E] = ARMORF(...) returns the final prediction error E (the variance
%   estimate of the white noise input to the AR model).
%
%   [A,E,K] = ARMORF(...) returns the vector K of reflection coefficients (parcor coefficients).
%
%   Ref: M. Morf, etal, Recursive Multichannel Maximum Entropy Spectral Estimation,
%              IEEE trans. GeoSci. Elec., 1978, Vol.GE-16, No.2, pp85-94.
%        S. Haykin, Nonlinear Methods of Spectral Analysis, 2nd Ed.
%              Springer-Verlag, 1983, Chapter 2
%
%   finished on Aug.9, 2002 by Yonghong Chen
%   bug fixes and revised on April, 2005 by Rajasimhan

% Initialization
[L,N]=size(x);
R0=zeros(L,L);
R0f=R0;
R0b=R0;
pf=R0;
pb=R0;
pfb=R0;
ap(:,:,1)=R0;
bp(:,:,1)=R0;
En=R0;

for i=1:ntrls
    En=En+x(:,(i-1)*npts+1:i*npts)*x(:,(i-1)*npts+1:i*npts)';
    ap(:,:,1)=ap(:,:,1)+x(:,(i-1)*npts+2:i*npts)*x(:,(i-1)*npts+2:i*npts)';
    bp(:,:,1)=bp(:,:,1)+x(:,(i-1)*npts+1:i*npts-1)*x(:,(i-1)*npts+1:i*npts-1)';
end

ap(:,:,1) = inv((chol(ap(:,:,1)/ntrls*(npts-1)))');
bp(:,:,1) = inv((chol(bp(:,:,1)/ntrls*(npts-1)))');

for i=1:ntrls
    efp = ap(:,:,1)*x(:,(i-1)*npts+2:i*npts);
    ebp = bp(:,:,1)*x(:,(i-1)*npts+1:i*npts-1);
    pf = pf + efp*efp';
    pb = pb + ebp*ebp';
    pfb = pfb + efp*ebp';
end
En = chol(En/N)'; % Covariance of the noise

% Initial output variables
coeff = [];%  Coefficient matrices of the AR model
kr=[];  % reflection coefficients

for m=1:p
    % Calculate the next order reflection (parcor) coefficient
    ck = inv((chol(pf))')*pfb*inv(chol(pb));
    kr=[kr,ck];
    % Update the forward and backward prediction errors
    ef = eye(L)- ck*ck';
    eb = eye(L)- ck'*ck;
    
    % Update the prediction error
    En = En*chol(ef)';
    E = (ef+eb)./2;
    
    % Update the coefficients of the forward and backward prediction errors
    ap(:,:,m+1) = zeros(L);
    bp(:,:,m+1) = zeros(L);
    pf = zeros(L);
    pb = zeros(L);
    pfb = zeros(L);
    
    for i=1:m+1
        a(:,:,i) = inv((chol(ef))')*(ap(:,:,i)-ck*bp(:,:,m+2-i));
        b(:,:,i) = inv((chol(eb))')*(bp(:,:,i)-ck'*ap(:,:,m+2-i));
    end
    for k=1:ntrls
        efp = zeros(L,npts-m-1);
        ebp = zeros(L,npts-m-1);
        for i=1:m+1
            k1=m+2-i+(k-1)*npts+1;
            k2=npts-i+1+(k-1)*npts;
            efp = efp+a(:,:,i)*x(:,k1:k2);
            ebp = ebp+b(:,:,m+2-i)*x(:,k1-1:k2-1);
        end
        pf = pf + efp*efp';
        pb = pb + ebp*ebp';
        pfb = pfb + efp*ebp';
    end
    ap = a;
    bp = b;
end
for j=1:p
    coeff = [coeff,inv(a(:,:,1))*a(:,:,j+1)];
end

varargout{1} = coeff;
if nargout >= 2
    varargout{2} = En*En';
end
if nargout >= 3
    varargout{3} = kr;
end
end
%-------------------------------armorf.m-------------------------------------------------------

function varargout = armorf0(x,Nr,Nl,p);
%ARMORF   AR parameter estimation via LWR method by Morf modified.
%   x is a matrix whose every row is one variable's time series
%   Nr is the number of realizations, Nl is the length of every realization
%   If the time series are stationary long, just let Nr=1, Nl=length(x)
%   p is the order of AR model
%
%   A = ARMORF(X,NR,NL,P) returns the polynomial coefficients A corresponding to 
%     the AR model estimate of matrix X using Morf's method.
%
%   [A,E] = ARMORF(...) returns the final prediction error E (the
%   covariance matrix of the white noise of the AR model).
%
%   [A,E,K] = ARMORF(...) returns the vector K of reflection 
%     coefficients (parcor coefficients).
%
%   Ref: M. Morf, etal, Recursive Multichannel Maximum Entropy Spectral Estimation,
%              IEEE trans. GeoSci. Elec., 1978, Vol.GE-16, No.2, pp85-94.
%        S. Haykin, Nonlinear Methods of Spectral Analysis, 2nd Ed.
%              Springer-Verlag, 1983, Chapter 2

% Initialization
[L,N]=size(x);
R0=zeros(L,L);
R0f=R0;
R0b=R0;
pf=R0;
pb=R0;
pfb=R0;
ap(:,:,1)=R0;
bp(:,:,1)=R0;
En=R0;
for i=1:Nr
    En=En+x(:,(i-1)*Nl+1:i*Nl)*x(:,(i-1)*Nl+1:i*Nl)';
    ap(:,:,1)=ap(:,:,1)+x(:,(i-1)*Nl+2:i*Nl)*x(:,(i-1)*Nl+2:i*Nl)';        
    bp(:,:,1)=bp(:,:,1)+x(:,(i-1)*Nl+1:i*Nl-1)*x(:,(i-1)*Nl+1:i*Nl-1)';
end
ap(:,:,1) = inv((chol(ap(:,:,1)/Nr*(Nl-1)))');
bp(:,:,1) = inv((chol(bp(:,:,1)/Nr*(Nl-1)))');
for i=1:Nr
    efp = ap(:,:,1)*x(:,(i-1)*Nl+2:i*Nl);
    ebp = bp(:,:,1)*x(:,(i-1)*Nl+1:i*Nl-1);
    pf = pf + efp*efp';
    pb = pb + ebp*ebp';
    pfb = pfb + efp*ebp';
end
En = chol(En/N)'; % Covariance of the noise

% Initial output variables
coeff = [];%  Coefficient matrices of the AR model
kr=[];  % reflection coefficients

for m=1:p
   % Calculate the next order reflection (parcor) coefficient
   ck = inv((chol(pf))')*pfb*inv(chol(pb));
   kr=[kr,ck];
   % Update the forward and backward prediction errors
   ef = eye(L)- ck*ck';
   eb = eye(L)- ck'*ck;
     
   % Update the prediction error
   En = En*chol(ef)';
   E = (ef+eb)./2;   
   
   % Update the coefficients of the forward and backward prediction errors
   ap(:,:,m+1) = zeros(L);
   bp(:,:,m+1) = zeros(L);
   pf = zeros(L);
   pb = zeros(L);
   pfb = zeros(L);

   for i=1:m+1       
       a(:,:,i) = inv((chol(ef))')*(ap(:,:,i)-ck*bp(:,:,m+2-i));
       b(:,:,i) = inv((chol(eb))')*(bp(:,:,i)-ck'*ap(:,:,m+2-i));
   end
   for k=1:Nr
       efp = zeros(L,Nl-m-1);
       ebp = zeros(L,Nl-m-1);
       for i=1:m+1
           k1=m+2-i+(k-1)*Nl+1;
           k2=Nl-i+1+(k-1)*Nl;
           efp = efp+a(:,:,i)*x(:,k1:k2);
           ebp = ebp+b(:,:,m+2-i)*x(:,k1-1:k2-1);
       end
       pf = pf + efp*efp';
       pb = pb + ebp*ebp';
       pfb = pfb + efp*ebp';
   end
   ap = a;
   bp = b;
end
for j=1:p
    coeff = [coeff,inv(a(:,:,1))*a(:,:,j+1)];
end

varargout{1} = coeff;
if nargout >= 2
    varargout{2} = En*En';
end
if nargout >= 3
    varargout{3} = kr;
end

end
%----------------------------- AZ2spectra.m --------------------------------    

function [S,H] = AZ2spectra(A,Z,p,freq,fs);
f_ind = 0;
for f = freq,
      f_ind = f_ind+1;
      [Stmp,Htmp] = spectrum(A,Z,p,f,fs);       
      H(:,:,f_ind) = Htmp;
      S(:,:,f_ind) = Stmp; %auto-& cross-spectra 
end
for ind = 1:size(Z,1),
      S(ind,ind,:) = real(S(ind,ind,:)); %avoiding numerical errors
end
end
%------------------------------spectrum.m-----------------------------------
function [S,H] = spectrum(A,Z,M,f,fs);
% Get the coherence spectrum
N = size(Z,1);
H = eye(N,N); % identity matrix
for m = 1 : M
    H = H + A(:,(m-1)*N+1:m*N)*exp(-i*m*2*pi*f/fs);
    % Multiply f in the exponent by sampling interval (=1/fs). See Shiavi
end
H = inv(H);
S = H*Z*H'/fs;
% One has to multiply HZH' by sampling interval (=1/fs)
%to get the properly normalized spectral density. See Shiavi.
%To get 1-sided power spectrum, multiply S by 2.
end
%%%%%%%%%-------------------------S2coh.m ------------------------
function coh = S2coh(S); 
%Input: S auto-& cross pectra in the form: frequency. channel. channel 
%Output: coh (Coherence) in the form: frequency. channel. channel 
%Written by M. Dhamala

Nc = size(S,2);
for ii = 1: Nc,
   for jj = 1: Nc,
       coh(:,ii,jj) = real(abs(S(:,ii,jj)).^2./(S(:,ii,ii).*S(:,jj,jj)));
   end
end
end
%----------------------- compute CGC given the spectra----------------------
function causality = getCGC(S,fs,freq);
%Input: S auto-& cross pectra in the form: frequency. channel. channel 
%Output: causality in the form: frequency. channel. channel
%Written and revised by M. Dhamala

[H0, SIGMA0] = wilson_sf(S, fs); Nc = size(S,1);
        for y = 1:Nc
            y_omit = [1:y-1 y+1:Nc];	
            [tmpHred0, tmpSIGMAred0] = wilson_sf(S(y_omit,y_omit,:), fs);
            Hred0     = zeros(Nc,Nc,length(freq)); 
            SIGMAred0 = zeros(Nc,Nc);
            Hred0(y_omit,y_omit,:)   = tmpHred0;
            SIGMAred0(y_omit,y_omit) = tmpSIGMAred0;
            for x = 1:Nc
                if x~=y
                    z = 1:Nc;
                    z([x y]) = []; 
                    %-- FULL model ----------------
                    xyz = [x y z];
                    SIGMA = SIGMA0(xyz,xyz);
                    H     = H0(xyz,xyz,:);
                    b1  = 1;
                    b2  = b1+1;
                    b3  = (b2+1):length(xyz);
                    b   = [numel(b1) numel(b2) numel(b3)];
                    tmp1 = -SIGMA(b2,b1)/SIGMA(b1,b1);
                    tmp2 = -SIGMA(b3,b1)/SIGMA(b1,b1);
                    tmp3 = -(SIGMA(b3,b2)+tmp2*SIGMA(b1,b2))/(SIGMA(b2,b2)+tmp1*SIGMA(b1,b2));
                    p1_1 = [eye(b(1))    zeros(b(1),b(2))    zeros(b(1),b(3));
                            tmp1         eye(b(2))           zeros(b(2),b(3));
                            tmp2         zeros(b(3),b(2))    eye(b(3))];
                    p1_2 = [eye(b(1))            zeros(b(1),b(2))    zeros(b(1),b(3));
                            zeros(b(2),b(1))     eye(b(2))           zeros(b(2),b(3));
                            zeros(b(3),b(1))     tmp3                eye(b(3))];
                    P1   = p1_2*p1_1;
                    %-- REDUCED model --------------
                    xz   = [x z];
                    SIGMAred = SIGMAred0(xz,xz);
                    Hred     = Hred0(xz,xz,:);
                    bx1  = 1;
                    bx2  = (bx1+1):length(xz);
                    bx   = [numel(bx1) numel(bx2)];
                    P2   = [eye(bx(1))                              zeros(bx(1),bx(2));
                            -SIGMAred(bx2,bx1)/SIGMAred(bx1,bx1)    eye(bx(2))];
                    %-- conditional GC -------------
                    for kk = 1:size(S,3)
                        HH = H(:,:,kk)/P1;
                        B  = P2/Hred(:,:,kk);
                        BB = [B(bx1,bx1)        zeros(b(1),b(2))	B(bx1,bx2);
                              zeros(b(2),b(1))  eye(b(2))           zeros(b(2),b(3));
                              B(bx2,bx1)        zeros(b(3),b(2))    B(bx2,bx2)];
                        FF = BB*HH;
                        numer = abs(det(SIGMAred(bx1,bx1)));
                        denom = abs(det(FF(b1,b1)*SIGMA(b1,b1)*conj(FF(b1,b1))));
                        causality(kk,y,x) = log(numer./denom);
                    end 
                elseif x==y
                    causality(:,y,x) = 0;   
                end
            end
        end
end
%--------------------------------------------------------------------------------      
%% wilson's method of spectral factorization
%-------------------------------------------------------------------------------
% [H Z, ps, ps0, converged] = wilson_sf(S, fs, tol)
% Performs a numerical inner-outer factorization of a spectral matrix, using
% Wilsons method. This implementation here is a slight modification of the
% original implemention by M. Dhamala (mdhamala@bme.ufl.edu) & G. Rangarajan
% (rangaraj@math.iisc.ernet.in), UF, Aug 3-4, 2006.
%
% modified by S K Mody (modysk@gmail.com), 22.Sept.2016
%revised by M. Dhamala, Oct, 2016
% REF:
% The Factorization of Matricial Spectral Densities, SIAM J. Appl. Math,
% Vol. 23, No. 4, pgs 420-426 December 1972 by G T Wilson).
%
% ARGS:
% S:
%	Spectral matrix function. This should be specified as a (k x k x m)
%	array for the frequencies in the closed range [0, 0.5], divided
%	into equal intervals. The spectrum for negative frequencies is assumed
%	to be symmetric, ie:-
%		S(-f) = transpose(S(f))
%
% fs:
%	Sampling rate. This is required for normalization.
%	***IMPORTANT: Please ensure that the spectral matrix input, S, has
%	been normalized by the same value of fs, otherwise the the output
%	Z will be off by a factor.
%
% tol [default: 1e-9]:
%	The tolerance with which to check for convergence. Iterations stop
%	either when the number of iterations reaches a prespecified maximum
%	or when all of the following conditions are satisfied:-
%		|(ps0 - ps0_prev)./ps0_prev| < tol
%		|(ps - ps_prev)./ps_prev| < tol
%		|(S - ps*ps')./S| < tol
%	where |.| is the max norm.
%
% OUTPUT:
% H, Z:
%	H is complex array of the same size as S, and Z is real symmetric
%	positive definite matrix such that for each i:
%		S(:,:,i) = H(:,:,i)*Z*H(:,:,i)'
%
% ps, ps0:
%	(k x k x m) complex array. Theoretically Ps is a function defined on
%	the on the boundary of the unit circle in the complex plane, such that:
%		S(:,:,i) = ps(:,:,i)*ps(:,:,i)'
%	Theoretically, Ps has a holomorphic extension in the complex plane to
%	all |z| < 1.). ps0 is the upper triangular matrix that is the value of
%	ps at the origin. Z is related to ps0 by:
%		Z = ps0*ps0'
%	
% converged:
%	Boolean value indicating whether the iteration converged to within the
%	specified tolerance.
%
% relerr:
%	The relative Cauchy error of the convergence of the spectrum or Ps.
%
function [H, Z, ps, ps0, converged, relerr] = wilson_sf(S, fs, tol)
	if (nargin < 3) || isempty(tol), tol = 1e-9; end
	assert(isscalar(fs) && (fs > 0), ...
		'fs must be a positive scalar value representing the sampling rate. ');
	
	[k, ~, N] = size(S);

	Sarr = cat(3, S, conj(S(:, :, N-1:-1:2)));
	ps0 = ps0_initial__(Sarr);
	ps = repmat(ps0, [1,1,N]);
	ps = cat(3, ps, conj(ps(:,:,N-1:-1:2)));
	M = size(Sarr, 3);

	I = eye(k);
	maxiter = min( 500, floor(sqrt(10/tol)) );

	U = zeros(size(Sarr));
	for j = 1 : M
		U(:,:,j) = chol(Sarr(:,:,j));
	end

	niter = 0;
	converged = false;
	g = zeros(k,k,M);
	while ( (niter < maxiter) && ~converged )
		for i = 1 : M
			% Equivalent to:
			% g(:,:,i) = ps(:,:,i)\Sarr(:,:,i)/ps(:,:,i)' + I;
			V = ps(:,:,i)\U(:,:,i)';
			g(:,:,i) = V*V' + I;
		end

		[gp, gp0] = PlusOperator(g);
		T = -tril(gp0, -1);
		T = T - T';

		ps_prev = ps;
		for i = 1 : M,
			ps(:,:,i) = ps(:,:,i)*(gp(:,:,i) + T);
		end

		ps0_prev = ps0;
		ps0 = ps0*(gp0 + T);

		% Relative cauchy error. Check on S is expensive, so check Ps0 first, then Ps and only then S.
		[converged relerr] = check_converged_ps__(ps, ps_prev, ps0, ps0_prev, tol);
		if converged
			% Uncomment this next line to check for relative cauchy error in spectrum.
			%[converged relerr] = check_converged_S__(Sarr, ps, tol);
		end

		niter = niter + 1;
	end

	H = zeros(k,k,N);
	for i = 1 : N
		H(:,:,i) = ps(:,:,i)/ps0;
	end
	
	ps = sqrt(fs)*ps(:,:,1:N);
	ps0 = sqrt(fs)*ps0;
	Z = ps0*ps0';
end

function ps0 = ps0_initial__(Sarr)
	[k, ~, M] = size(Sarr);
	
	% perform ifft to obtain gammas.
	Sarr = reshape(Sarr, [k*k, M]);
	gamma = ifft(transpose(Sarr));
	gamma0 = gamma(1,:);
	gamma0 = reshape(gamma0, [k k]);
	
	% Remove any assymetry due to rounding error.
	% This also will zero out any imaginary values
	% on the diagonal - real diagonals are required for cholesky.
	gamma0 = real((gamma0 + gamma0')/2);

	ps0 = chol(gamma0);
	
end

%% This function is for [ ]+operation
function [gp, gp0] = PlusOperator(g)

	[k, ~, M] = size(g);
	N = ceil( (M+1)/2 );
	
	g = reshape(g, [k*k, M]);
	gammma = real(ifft(transpose(g)));
	gammma = reshape(transpose(gammma), [k,k,M]);
	
	% Take half of the zero lag
	gammma(:,:,1) = 0.5*gammma(:,:,1);
	gp0 = gammma(:,:,1);
	
	% Zero out negative powers.
	gammma(:, :, N+1:end) = 0;

	% Reconstitute
	gammma = reshape(gammma, [k*k, M]);
	gp = fft(transpose(gammma));
	gp = reshape(transpose(gp), [k,k,M]);
	
end
%%

function [converged_ps relerr] = check_converged_ps__(ps, ps_prev, ps0, ps0_prev, tol)

	[converged_ps relerr] = CheckRelErr__(ps0, ps0_prev, tol);
	if converged_ps
		[converged_ps RelErr2] = CheckRelErr__(ps, ps_prev, tol);
		relerr = max(relerr, RelErr2);
	end
	
end
%%

function [converged_S relerr] = check_converged_S__(S, ps, tol)

	FX = zeros(size(ps));
	parfor j = 1 : size(ps,3)
		FX(:,:,j) = ps(:,:,j)*ps(:,:,j)';
	end
	[converged_S relerr] = CheckRelErr__(FX, S, tol);
	
end
%%

function [ok, relerr] = CheckRelErr__(A,B,reltol)
	D = abs(B - A);
	
	A = abs(A);
	A(A <= 2*eps) = 1; % Minimum detectable difference between
							 % x and a value close to x is O(x)*eps.
	E = D./A;
	relerr = max(E(:));
	
	ok = (relerr <= reltol);
end
%%

