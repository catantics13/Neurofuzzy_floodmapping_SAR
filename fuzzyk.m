%% Read files

val=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulm3.asc');
sar=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmsar3.asc');
ctex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmctex3.asc');
optex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmoptex3.asc');

sar=sar(:,3:2278);
ctex=ctex(:,3:2278);
optex=optex(:,3:2278);

%% create fuzzy vectors
[m,n]=size(val);
val_kf=zeros(m,n,2);
val_kf(:,:,1)=val;
val_kf(:,:,2)=1-val;

sar_kf=zeros(m,n,2);
sar_kf(:,:,1)=sar;
sar_kf(:,:,2)=1-sar;

ctex_kf=zeros(m,n,2);
ctex_kf(:,:,1)=ctex;
ctex_kf(:,:,2)=1-ctex;

optex_kf=zeros(m,n,2);
optex_kf(:,:,1)=optex;
optex_kf(:,:,2)=1-optex;


%% preallocate vectors
val_crisp=zeros(m,n,2);
sar_crisp=zeros(m,n,2);
ctex_crisp=zeros(m,n,2);
optex_crisp=zeros(m,n,2);

val_nbh=zeros(m,n,2);
sar_nbh=zeros(m,n,2);
ctex_nbh=zeros(m,n,2);
optex_nbh=zeros(m,n,2);

p_val_nbh=zeros(m,n,2);
p_sar_nbh=zeros(m,n,2);
p_ctex_nbh=zeros(m,n,2);
p_optex_nbh=zeros(m,n,2);

sar_sa=zeros(m,n);
sar_sb=zeros(m,n);
sar_ss=zeros(m,n);

ctex_sa=zeros(m,n);
ctex_sb=zeros(m,n);
ctex_ss=zeros(m,n);

optex_sa=zeros(m,n);
optex_sb=zeros(m,n);
optex_ss=zeros(m,n);

c_delta_val=zeros(m,n,2);
c_delta_sar=zeros(m,n,2);
c_delta_ctex=zeros(m,n,2);
c_delta_optex=zeros(m,n,2);

sar_fcr=zeros(m,n);
ctex_fcr=zeros(m,n);
optex_fcr=zeros(m,n);

p_val=zeros(2,1);
p_sar=zeros(2,1);
p_ctex=zeros(2,1);
p_optex=zeros(2,1);

p_flood_val=zeros(2,1);
p_flood_sar=zeros(2,1);
p_flood_ctex=zeros(2,1);
p_flood_optex=zeros(2,1);

sar_e=zeros(m,n);
ctex_e=zeros(m,n);
optex_e=zeros(m,n);


%% fuzzy map comparison


%nbh=2^(-d/2)---> exp. decay

%%Hagen 2003

for j=1:n
    for i=1:m
        % calculation of cell-by-cell similarity
        sar_fcr(i,j)=max(min(val_kf(i,j,1),sar_kf(i,j,1)),min(val_kf(i,j,2),sar_kf(i,j,2)));
        ctex_fcr(i,j)=max(min(val_kf(i,j,1),ctex_kf(i,j,1)),min(val_kf(i,j,2),ctex_kf(i,j,2)));
        optex_fcr(i,j)=max(min(val_kf(i,j,1),optex_kf(i,j,1)),min(val_kf(i,j,2),optex_kf(i,j,2)));
        % calculation of neighbourhood vectors
         
        % calculation of crisp sets
        if val_kf(i,j,1)>val_kf(i,j,2)
            val_crisp(i,j,1)=1;
            val_crisp(i,j,2)=0;
        else
            val_crisp(i,j,1)=0;
            val_crisp(i,j,2)=1;
        end
        
        if sar_kf(i,j,1)>sar_kf(i,j,2)
            sar_crisp(i,j,1)=1;
            sar_crisp(i,j,2)=0;
        else
            sar_crisp(i,j,1)=0;
            sar_crisp(i,j,2)=1;
        end
        
        if ctex_kf(i,j,1)>ctex_kf(i,j,2)
            ctex_crisp(i,j,1)=1;
            ctex_crisp(i,j,2)=0;
        else
            ctex_crisp(i,j,1)=0;
            ctex_crisp(i,j,2)=1;
        end
        
        if optex_kf(i,j,1)>optex_kf(i,j,2)
            optex_crisp(i,j,1)=1;
            optex_crisp(i,j,2)=0;
        else
            optex_crisp(i,j,1)=0;
            optex_crisp(i,j,2)=1;
        end

    end
end
        
for j=1:n-2
    for i=1:m-2
        A=[val_kf(i+1,j+2,1)*0.5,val_kf(i,j+2,1)*0.25,val_kf(i,j+1,1)*0.5,...
            val_kf(i,j,1)*0.25,val_kf(i+1,j,1)*0.5,val_kf(i+2,j,1)*0.25,val_kf(i+2,j+1,1)*0.5,...
            val_kf(i+2,j+2,1)*0.25,val_kf(i+1,j+1,1)*1];
        val_nbh(i+1,j+1,1)=max(A);
        B=[val_kf(i+1,j+2,2)*0.5,val_kf(i,j+2,2)*0.25,val_kf(i,j+1,2)*0.5,...
            val_kf(i,j,2)*0.25,val_kf(i+1,j,2)*0.5,val_kf(i+2,j,2)*0.25,val_kf(i+2,j+1,2)*0.5,...
            val_kf(i+2,j+2,2)*0.25,val_kf(i+1,j+1,2)*1];
        val_nbh(i+1,j+1,2)=max(B);
        
        C=[sar_kf(i+1,j+2,1)*0.5,sar_kf(i,j+2,1)*0.25,sar_kf(i,j+1,1)*0.5,...
            sar_kf(i,j,1)*0.25,sar_kf(i+1,j,1)*0.5,sar_kf(i+2,j,1)*0.25,sar_kf(i+2,j+1,1)*0.5,...
            sar_kf(i+2,j+2,1)*0.25,sar_kf(i+1,j+1,1)*1];
        sar_nbh(i+1,j+1,1)=max(C);
        
        D=[sar_kf(i+1,j+2,2)*0.5,sar_kf(i,j+2,2)*0.25,sar_kf(i,j+1,2)*0.5,...
            sar_kf(i,j,2)*0.25,sar_kf(i+1,j,2)*0.5,sar_kf(i+2,j,2)*0.25,sar_kf(i+2,j+1,2)*0.5,...
            sar_kf(i+2,j+2,2)*0.25,sar_kf(i+1,j+1,2)*1];
        sar_nbh(i+1,j+1,2)=max(D);
        
        
        E=[ctex_kf(i+1,j+2,1)*0.5,ctex_kf(i,j+2,1)*0.25,ctex_kf(i,j+1,1)*0.5,...
            ctex_kf(i,j,1)*0.25,ctex_kf(i+1,j,1)*0.5,ctex_kf(i+2,j,1)*0.25,ctex_kf(i+2,j+1,1)*0.5,...
            ctex_kf(i+2,j+2,1)*0.25,ctex_kf(i+1,j+1,1)*1];
        ctex_nbh(i+1,j+1,1)=max(E);
        F=[ctex_kf(i+1,j+2,2)*0.5,ctex_kf(i,j+2,2)*0.25,ctex_kf(i,j+1,2)*0.5,...
            ctex_kf(i,j,2)*0.25,ctex_kf(i+1,j,2)*0.5,ctex_kf(i+2,j,2)*0.25,ctex_kf(i+2,j+1,2)*0.5,...
            ctex_kf(i+2,j+2,2)*0.25,ctex_kf(i+1,j+1,2)*1];
        ctex_nbh(i+1,j+1,2)=max(F);
        
        G=[optex_kf(i+1,j+2,1)*0.5,optex_kf(i,j+2,1)*0.25,optex_kf(i,j+1,1)*0.5,...
            optex_kf(i,j,1)*0.25,optex_kf(i+1,j,1)*0.5,optex_kf(i+2,j,1)*0.25,optex_kf(i+2,j+1,1)*0.5,...
            optex_kf(i+2,j+2,1)*0.25,optex_kf(i+1,j+1,1)*1];
        optex_nbh(i+1,j+1,1)=max(G);
        
        H=[optex_kf(i+1,j+2,2)*0.5,optex_kf(i,j+2,2)*0.25,optex_kf(i,j+1,2)*0.5,...
            optex_kf(i,j,2)*0.25,optex_kf(i+1,j,2)*0.5,optex_kf(i+2,j,2)*0.25,optex_kf(i+2,j+1,2)*0.5,...
            optex_kf(i+2,j+2,2)*0.25,optex_kf(i+1,j+1,2)*1];
        optex_nbh(i+1,j+1,2)=max(H);
    end
end

for j=1:n
    for i=1:m
        % Calculation of spatial similarity
        sar_sa(i,j)=max(min(val_nbh(i,j,1),sar_kf(i,j,1)),min(val_nbh(i,j,2),sar_kf(i,j,2)));
        sar_sb(i,j)=max(min(val_kf(i,j,1),sar_nbh(i,j,1)),min(val_kf(i,j,2),sar_nbh(i,j,2)));
        sar_ss(i,j)=min(sar_sa(i,j),sar_sb(i,j));
        
        ctex_sa(i,j)=max(min(val_nbh(i,j,1),ctex_kf(i,j,1)),min(val_nbh(i,j,2),ctex_kf(i,j,2)));
        ctex_sb(i,j)=max(min(val_kf(i,j,1),ctex_nbh(i,j,1)),min(val_kf(i,j,2),ctex_nbh(i,j,2)));
        ctex_ss(i,j)=min(ctex_sa(i,j),ctex_sb(i,j));
        
        optex_sa(i,j)=max(min(val_nbh(i,j,1),optex_kf(i,j,1)),min(val_nbh(i,j,2),optex_kf(i,j,2)));
        optex_sb(i,j)=max(min(val_kf(i,j,1),optex_nbh(i,j,1)),min(val_kf(i,j,2),optex_nbh(i,j,2)));
        optex_ss(i,j)=min(optex_sa(i,j),optex_sb(i,j));
    end
end
   
 
%% fuzzy kappa calculation


for k=1:2
    p_val(k)=(sum(reshape(val_crisp(:,:,k),m*n,1)))/(m*n);
    p_sar(k)=(sum(reshape(sar_crisp(:,:,k),m*n,1)))/(m*n);
    p_ctex(k)=(sum(reshape(ctex_crisp(:,:,k),m*n,1)))/(m*n);
    p_optex(k)=(sum(reshape(optex_crisp(:,:,k),m*n,1)))/(m*n);
end


for j=1:n
    for i=1:m
        for k=1:2
            if abs(val_kf(i,j,k)-p_val(k))<0.001
                c_delta_val(i,j,k)=1;
            else
                c_delta_val(i,j,k)=0;
            end
            if abs(sar_kf(i,j,k)-p_sar(k))<0.001
                c_delta_sar(i,j,k)=1;
            else
                c_delta_sar(i,j,k)=0;
            end
            if abs(ctex_kf(i,j)-p_ctex(k))<0.001
                c_delta_ctex(i,j,k)=1;
            else
                c_delta_ctex(i,j,k)=0;
            end
            if abs(optex_kf(i,j,k)-p_optex(k))<0.001
                c_delta_optex(i,j,k)=1;
            else
                c_delta_optex(i,j,k)=0;
            end
        end
    end
end


for k=1:2
    p_flood_val(k)=sum(reshape(c_delta_val(:,:,k),m*n,1))/(m*n);
    p_flood_sar(k)=sum(reshape(c_delta_sar(:,:,k),m*n,1))/(m*n);
    p_flood_ctex(k)=sum(reshape(c_delta_ctex(:,:,k),m*n,1))/(m*n);
    p_flood_optex(k)=sum(reshape(c_delta_optex(:,:,k),m*n,1))/(m*n);
end

for j=1:n-2
    for i=1:m-2
        for k=1:2
            Q=[c_delta_val(i+1,j+2,k),c_delta_val(i,j+2,k),c_delta_val(i,j+1,k),...
            c_delta_val(i,j,k),c_delta_val(i+1,j,k),c_delta_val(i+2,j,k),c_delta_val(i+2,j+1,k),...
            c_delta_val(i+2,j+2,k),c_delta_val(i+1,j+1,k)];
            p_val_nbh(i,j,k)=sum(Q)/size(Q,1);
            A=[c_delta_sar(i+1,j+2,k),c_delta_sar(i,j+2,k),c_delta_sar(i,j+1,k),...
            c_delta_sar(i,j,k),c_delta_sar(i+1,j,k),c_delta_sar(i+2,j,k),c_delta_sar(i+2,j+1,k),...
            c_delta_sar(i+2,j+2,k),c_delta_sar(i+1,j+1,k)];
            p_sar_nbh(i,j,k)=(sum(A)/size(A,1))*p_val_nbh(i,j,k);
            B=[c_delta_ctex(i+1,j+2,k),c_delta_ctex(i,j+2,k),c_delta_ctex(i,j+1,k),...
            c_delta_ctex(i,j,k),c_delta_ctex(i+1,j,k),c_delta_ctex(i+2,j,k),c_delta_ctex(i+2,j+1,k),...
            c_delta_ctex(i+2,j+2,k),c_delta_ctex(i+1,j+1,k)];
            p_ctex_nbh(i,j,k)=(sum(B)/size(B,1))*p_val_nbh(i,j,k);
            C=[c_delta_optex(i+1,j+2,k),c_delta_optex(i,j+2,k),c_delta_optex(i,j+1,k),...
            c_delta_optex(i,j,k),c_delta_optex(i+1,j,k),c_delta_optex(i+2,j,k),c_delta_optex(i+2,j+1,k),...
            c_delta_optex(i+2,j+2,k),c_delta_optex(i+1,j+1,k)];
            p_optex_nbh(i,j,k)=(sum(C)/size(C,1))*p_val_nbh(i,j,k);
        end
    end
end
y=((n+1)/2);
x=((m+1)/2);
d=zeros(x,1);
% for jj=1:x
%     for ii=1:x
%         for k=1:2
%         Q=[c_delta_val(x,y+jj,k),c_delta_val(x-ii,y+jj,k),c_delta_val(x-ii,y,k),...
%             c_delta_val(x-ii,y-jj,k),c_delta_val(x,y-jj,k),c_delta_val(x+ii,y-jj,k),...
%             c_delta_val(x+ii,y,k), c_delta_val(x+ii,y+jj,k),c_delta_val(x,y,k)];
%         d(ii)=(-2^(ii/2));
%         p_sar(ii)=sum
%         
%         end
%     end
% end
val_v=reshape(val_kf,m*n,2);
sar_v=reshape(sar_kf,m*n,2);
ctex_v=reshape(ctex_kf,m*n,2);
optex_v=reshape(optex_kf,m*n,2);

[idx,D]= rangesearch(val_v,sar_v(x*y,:),20,'Distance','euclidean');
indexes=[idx{1,1}]';
distances=[D{1,1}]';
c_delta_val_v=reshape(c_delta_val,m*n,2);
c_delta_sar_v=reshape(c_delta_sar,m*n,2);
c_delta_ctex_v=reshape(c_delta_ctex,m*n,2);
c_delta_optex_v=reshape(c_delta_optex,m*n,2);

f=inline(vectorize('exp(log(1/2)*d/2)'),'d');

d_weights=f(distances);

count=1;

p=zeros(length(distances),1);
num_cells=zeros(length(distances),1);
cum_num_cells=zeros(length(distances),1);

nbh_sim_val=zeros(length(distances),2);
nbh_sim_sar=zeros(length(distances),2);
nbh_sim_ctex=zeros(length(distances),2);
nbh_sim_optex=zeros(length(distances),2);

loc_sim_sar=zeros(length(distances),1);
loc_sim_ctex=zeros(length(distances),1);
loc_sim_optex=zeros(length(distances),1);

prob_sar=zeros(length(distances),1);
prob_ctex=zeros(length(distances),1);
prob_optex=zeros(length(distances),1);


syms d

for ii=1:length(distances)-1
    if distances(ii)~=distances(ii+1)
       index=rangesearch(val_v,sar_v(x*y,:),distances(ii),'Distance',...
           'euclidean');
       index=[index{1,1}]';
       
       num_cells(ii)=length(index);
       if distances(ii)==0
           cum_num_cells(ii)=0;
       else
           cum_num_cells(ii)=sum(num_cells);
       end
           for k=1:2
               nbh_sim_val(ii,k)=max(val_v(index,k).*d_weights(ii));
               nbh_sim_sar(ii,k)=max(sar_v(index,k).*d_weights(ii));
               nbh_sim_ctex(ii,k)=max(ctex_v(index,k).*d_weights(ii));
               nbh_sim_optex(ii,k)=max(optex_v(index,k).*d_weights(ii));
               loc_sim_sar(ii,k)=min(nbh_sim_val(ii,k),nbh_sim_sar(ii,k));
               loc_sim_ctex(ii,k)=min(nbh_sim_val(ii,k),nbh_sim_sar(ii,k)); 
               loc_sim_optex(ii,k)=min(nbh_sim_val(ii,k),nbh_sim_sar(ii,k));
           end
        
       prob_sar(ii)=((sum(c_delta_sar_v(index,1))/length(index))*...
           (sum(c_delta_sar_v(index,2))/length(index)))...
           .*((sum(c_delta_val_v(index,1))/length(index))*...
           (sum(c_delta_val_v(index,2))/length(index)));
       prob_ctex(ii)=((sum(c_delta_ctex_v(index,1))/length(index))*...
           (sum(c_delta_ctex_v(index,2))/length(index)).*...
           (sum(c_delta_val_v(index,1))/length(index))*...
           (sum(c_delta_val_v(index,2))/length(index)));
       prob_optex(ii)=((sum(c_delta_optex_v(index,1))/length(index))*...
           (sum(c_delta_optex_v(index,2))/length(index)).*...
           (sum(c_delta_val_v(index,1))/length(index))*...
           (sum(c_delta_val_v(index,2))/length(index)));
    end
end
                
               
       

exp_sim_sar=sum(loc_sim_sar.*prob_sar);
exp_sim_ctex=sum(loc_sim_ctex.*prob_ctex);
exp_sim_optex=sum(loc_sim_optex.*prob_optex);

fuzzy_kappa_sar=(mean(nanmean(sar_ss))-exp_sim_sar)/(1-exp_sim_sar);
fuzzy_kappa_ctex=(mean(nanmean(ctex_ss))-exp_sim_ctex)/(1-exp_sim_ctex);
fuzzy_kappa_optex=(mean(nanmean(optex_ss))-exp_sim_optex)/(1-exp_sim_optex);

fuzzy_kappa_sar1=(0.879-exp_sim_sar)/(1-exp_sim_sar);
fuzzy_kappa_ctex1=(0.89-exp_sim_ctex)/(1-exp_sim_ctex);
fuzzy_kappa_optex1=(0.903-exp_sim_optex)/(1-exp_sim_optex);
