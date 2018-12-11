function [M,S] = compute_allpCGCvar3(X,fs,f,p);
%Usage: [M] = compute_pCGCvar3(X,fs,freq,p); This is a
%parametric method. This routine computes conditional 
%Granger causality (gc), coherence (coherence) and power (power)
% given the trivariate data (var3) X in the form of time x trial x channel
%-------------------------------------------------------------------------------------------------
%Inputs: X = trivariate data in the form of 3D matrix  for time. trial. channel.
%               fs = data sampling rate in Hz
%               f = frequencies at which GC is computed, eg. f = 0:1:fs/2;
%               p = model order (e.g. 3)             
%Outputs: M.freq = freq, M.gc = conditional Granger causality (i to j on k)
% M.coh = coherence, M.pow = power spectral density, M.freq = f, M.fs = fs;
              
%-------------------------------------------------------------------------------------------------
%Ref : M. Dhamala, et al. NeuroImage (2008). Written by M. Dhamala
%Revised by M. Dhamala on Aug, 2017 
%-------------------------------------------------------------------------------------------------
[Nt, Ntr,Nc] = size(X); %Nt = number of timepoints, Ntr = trials, Nc = channels 

x = reshape(X,Nt*Ntr,Nc); x = x'; % x in the form of channel. (time x trials)

[A, Z]=armorf(x,Ntr,Nt,p); %parameters by autoregressive fitting

[S,H] = AZ2spectra(A,Z,p,f,fs);

spectra = permute(S,[3 1 2]);coherence = S2coh(spectra);
for ichan = 1: Nc,
    power(:,ichan) = 2*spectra(:,ichan,ichan);%one-sided power
end

[H, Z] = wilson_sf(S,fs);

cgc = hz2cgcAll(S,H,Z,f,fs);

M.freq = f; M.pow = power; 
M.coh = coherence; M.gc=cgc; 
M.fs = fs;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------
% Functions below : armorf.m, AZ2spectra.m, spectrum.m, S2coh, hz2cdcausality.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------

%% ----------------- armorf.m ------------------------------------
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
%% ----------------------------- AZ2spectra.m --------------------

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
%------------------------------spectrum.m----------------------------------
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
%% ----------------------S2coh.m -----------
function coh = S2coh(S); 
%Input: S auto-& cross pectra in the form: frequency. channel. channel 
%Output: coh (Coherence) in the form: frequency. channel. channel 
%M. Dhamala, UF, August 2006.

Nc = size(S,2);
for ii = 1: Nc,
   for jj = 1: Nc,
       coh(:,ii,jj) = real(abs(S(:,ii,jj)).^2./(S(:,ii,ii).*S(:,jj,jj)));
   end
end
end
%% ------------------------ hz2cgcALL.m ------------------
function cgc = hz2cgcAll(S,H,Z,freq,fs);
%Revised by M. Dhamala, Aug 2017   
% Get Granger causality for all pairs
index = 0; nc = size(H,1);
cgc=zeros(nc,nc); 
 for i = 1:nc,
     for j = 1:nc,
         if i~=j
         [Zij,Zijk,Zkij,Zkk] = pttmatrx(Z,i,j);
         Z3=[Zij,Zijk;Zkij,Zkk];
         cc=Z3(1,1)*Z3(2,2)-Z3(1,2)^2;
         P=[1,0,zeros(1,nc-2);-Z3(1,2)/Z3(1,1),1,zeros(1,nc-2);(Z3(1,2)*Z3(2,3:end)-Z3(2,2)*Z3(1,3:end))'./cc,(Z3(1,2)*Z3(1,3:end)-Z3(1,1)*Z3(2,3:end))'./cc,eye(nc-2)];
         ZZ=P*Z3*P';
         f_ind = 0;
         for f = freq
             f_ind = f_ind+1;
             [Sij,Sijk,Skij,Skk]=pttmatrx(squeeze(S(:,:,f_ind)),i,j);
             [Hij,Hijk,Hkij,Hkk]=pttmatrx(squeeze(H(:,:,f_ind)),i,j);
             S3(:,:,f_ind)=[Sij,Sijk;Skij,Skk];
             S2(:,:,f_ind)=[Sij(1,1),Sijk(1,:);Skij(:,1),Skk];
             H3(:,:,f_ind)=[Hij,Hijk;Hkij,Hkk];
         end
         [H2,Z2] = wilson_sf(S2,fs);
         f_ind = 0;
         for f = freq
             f_ind = f_ind+1;
             %HH=squeeze(H3(:,:,f_ind))*inv(P);
             HH=squeeze(H3(:,:,f_ind))/P;%right slash instead of inv
             Q=[1,zeros(1,nc-2);-Z2(1,2:end)'./Z2(1,1),eye(nc-2)];
             %B=Q*inv(squeeze(H2(:,:,f_ind)));
             B=Q/squeeze(H2(:,:,f_ind));% right slash instead of inv
             BB=[B(1,1),0,B(1,2:end);0,1,zeros(1,nc-2);B(2:end,1),zeros(nc-2,1),B(2:end,2:end)];
             FF=BB*HH;
             cgc(j,i,f_ind)=log(abs(Z2(1,1))/abs(FF(1,1)*Z3(1,1)*conj(FF(1,1))));
         end
         end
     end
 end
 cgc = permute(cgc,[3 1 2]);
end

%---------------------------- hz2cdcausality.m ----------------------------
function Cj2i  = hz2cdcausality(H,Z,freq);
%Usage: Cj2i = hz2cdcausality(H,Z,freq);
%Inputs: H (transfer function of channel 1 to 3);
%      : Z (noise covariance of channel 1 to 3);
%Output: Cj2i, causality value of channel 2 driving 1 conditional on 3
%For example, if i = 1, j = 2, k = 3, then Cj2i is causality from 2 to 1 conditional on 3.
%Reference: Y.Chen, S.L. Bressler and M. Ding, J. Neurosc. Methods 150, 228 (2006).
% Written by M. Dhamala, UF, August 2006.
%--------------------------------------------------------------------------
f_ind  = 0; Z3 = Z; 
for f  = freq,
        f_ind = f_ind + 1;
        H3 = squeeze(H(:,:,f_ind));
        cc=Z3(1,1)*Z3(2,2)-Z3(1,2)^2;
        P=[1,0,0;-Z3(1,2)/Z3(1,1),1,0;(Z3(1,2)*Z3(2,3)-Z3(2,2)*Z3(1,3))/cc,(Z3(1,2)*Z3(1,3)-Z3(1,1)*Z3(2,3))/cc,1];
        ZZ=P*Z3*P';
        HH=H3*inv(P);%inv(P*inv(H3));
        [g1,g2,g3,g4]=pttmatrx(H3,1,3); %using partition matrix to get fitted model for two time series
        [e1,e2,e3,e4]=pttmatrx(Z3,1,3);
        gg=inv(g1)*g2;
        Z2=e1+2*real(gg)*e3+gg*e4*gg';
        H2=g1;
        S2=H2*Z2*H2';
        Q=[1,0;-Z2(1,2)/Z2(1,1),1];
        B=Q*inv(H2);
        BB=[B(1,1),0,B(1,2);0,1,0;B(2,1),0,B(2,2)];
        FF=BB*HH;
        Cj2i(f_ind)=log(abs(Z2(1,1))/abs(FF(1,1)*Z3(1,1)*conj(FF(1,1)))); %Geweke's original measure
        Ij2i(f_ind)=abs((Z2(1,1)-FF(1,1)*Z3(1,1)*conj(FF(1,1))))/abs(Z2(1,1)); %measure within [0,1]  
        %if causality values are small, Cj2i is almost equal to Ij2i.
end
end
%----------------------------- pttmatrx.m ---------------------------------
function [B11,B12,B21,B22]=pttmatrx(H,i,j)
% partition matrix into 2*2 blocks
L=length(H);
B1=H;
B2=H;
B11=reshape([H(i,i),H(j,i),H(i,j),H(j,j)],2,2);
B1(i,:)=[];
if j<i
    B1(j,:)=[];
else
    B1(j-1,:)=[];
end
B21=[B1(:,i),B1(:,j)];
B1(:,i)=[];
if j<i
    B1(:,j)=[];
else
    B1(:,j-1)=[];
end
B2(:,i)=[];
if j<i
    B2(:,j)=[];
else
    B2(:,j-1)=[];
end
B12=[B2(i,:);B2(j,:)];
B22=B1;
end
%% --------------------- wilson_sf.m -----------------------------
%-------- updated code for Wilson's method of spectral factorization ------
% [H Z, ps, ps0, converged] = wilson_sf(S, fs, tol)
% Performs a numerical inner-outer factorization of a spectral matrix, using
% Wilsons method. This implementation here is a slight modification of the
% original implemention by M. Dhamala (mdhamala@bme.ufl.edu) & G. Rangarajan
% (rangaraj@math.iisc.ernet.in), UF, Aug 3-4, 2006.
%
% modified by S K Mody (modysk@gmail.com), 22.Sept.2016
% revised by M. Dhamala, June, 2017
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

%----------------------------------- This function is for [ ]+operation
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

%--------------------------------
function [converged_ps relerr] = check_converged_ps__(ps, ps_prev, ps0, ps0_prev, tol)

	[converged_ps relerr] = CheckRelErr__(ps0, ps0_prev, tol);
	if converged_ps
		[converged_ps RelErr2] = CheckRelErr__(ps, ps_prev, tol);
		relerr = max(relerr, RelErr2);
	end
	
end

%-------------------------------
function [converged_S relerr] = check_converged_S__(S, ps, tol)

	FX = zeros(size(ps));
	parfor j = 1 : size(ps,3)
		FX(:,:,j) = ps(:,:,j)*ps(:,:,j)';
	end
	[converged_S relerr] = CheckRelErr__(FX, S, tol);
	
end

%-----------------------------
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
