function [S,H] = spectrum(A,Z,M,f,fs);% Get the coherence spectrum
N = size(Z,1);
H = eye(N,N); % identity matrix
for m = 1 : M
    H = H + A(:,(m-1)*N+1:m*N)*exp(-i*m*2*pi*f/fs);
    % Multiply f in the exponent by sampling interval (=1/fs). See Shiaviend
H = inv(H);S = H*Z*H'/fs; 
% One has to multiply HZH' by sampling interval (=1/fs)
%to get the properly normalized spectral density. See Shiavi.
%To get 1-sided power spectrum, multiply S by 2.
return
