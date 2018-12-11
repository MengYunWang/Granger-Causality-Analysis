%% ------------------------ hz2cgcALL.m --------------------------
% the main code for computing the gc values.

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
             %HH=squeeze(H3(:,:,f_ind))*inv(P);%inv(P*inv(H3));
             HH=squeeze(H3(:,:,f_ind))/P;% right slash instead of inv
             Q=[1,zeros(1,nc-2);-Z2(1,2:end)'./Z2(1,1),eye(nc-2)];
             %B=Q*inv(squeeze(H2(:,:,f_ind)));
             B=Q/squeeze(H2(:,:,f_ind));%right slash instead of inv
             BB=[B(1,1),0,B(1,2:end);0,1,zeros(1,nc-2);B(2:end,1),zeros(nc-2,1),B(2:end,2:end)];
             FF=BB*HH;
             cgc(j,i,f_ind)=log(abs(Z2(1,1))/abs(FF(1,1)*Z3(1,1)*conj(FF(1,1))));
         end
         end
     end
 end
 cgc = permute(cgc,[3 1 2]);
end

%--------------------------------------------------------------------------
function Cj2i  = hz2cdcausality(H,Z,freq);
%Usage: Cj2i = hz2cdcausality(H,Z,freq);
%Inputs: H (transfer function of channel 1 to 3);
%      : Z (noise covariance of channel 1 to 3);
%Output: Cj2i, causality value of channel 2 driving 1 conditional on 3
%For example, if i = 1, j = 2, k = 3, then Cj2i is causality from 2 to 1 conditional on 3.
%Reference: Y.Chen, S.L. Bressler and M. Ding, J. Neurosc. Methods 150, 228 (2006).
% Written by M. Dhamala, UF, August 2006.
%
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

%--------------------------------------------------------------------------
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