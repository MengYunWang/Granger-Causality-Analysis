%% ------------------------S2coh.m -----------------------

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
