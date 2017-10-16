%% Read files

val=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhil3.asc');
sar=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilsar3.asc');
ctex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex3.asc');
optex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhiloptex3.asc');
% 
sar=sar(:,3:2278);
ctex=ctex(:,3:2278);
optex=optex(:,3:2278);

%% create fuzzy vectors
[m,n]=size(val);
val_kf=zeros(m,n,2);
val_kf(:,:,1)=val;
val_kf(:,:,2)=1-val;
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulm3_kf.asc',val_kf);
sar_kf=zeros(m,n,2);
sar_kf(:,:,1)=sar;
sar_kf(:,:,2)=1-sar;
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmsar3_kf.asc',sar_kf);
ctex_kf=zeros(m,n,2);
ctex_kf(:,:,1)=ctex;
ctex_kf(:,:,2)=1-ctex;
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmctex3_kf.asc',ctex_kf);
optex_kf=zeros(m,n,2);
optex_kf(:,:,1)=optex;
optex_kf(:,:,2)=1-optex;
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhiloptex3_kf.asc',optex_kf);
% imshow(val_kf(:,:,2));
% info=geotiffinfo('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex3_ss.tif');
% [~,R]=geotiffread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex3_ss.tif');
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilval_kf.tif', ...
%      val_kf, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilsar_kf.tif', ...
%      sar_kf, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex_kf.tif', ...
%      ctex_kf, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhiloptex_kf.tif', ...
%      optex_kf, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

%% preallocate vectors
val_crisp=zeros(m,n,2);
sar_crisp=zeros(m,n,2);
ctex_crisp=zeros(m,n,2);
optex_crisp=zeros(m,n,2);
val_nbh=zeros(m,n,2);
sar_nbh=zeros(m,n,2);
ctex_nbh=zeros(m,n,2);
optex_nbh=zeros(m,n,2);

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
%% fuzzy map comparison
sar_fcr=zeros(m,n);
ctex_fcr=zeros(m,n);
optex_fcr=zeros(m,n);

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
        sar_sa(i,j)=max(min(val_nbh(i,j,1),sar_crisp(i,j,1)),min(val_nbh(i,j,2),sar_crisp(i,j,2)));
        sar_sb(i,j)=max(min(val_crisp(i,j,1),sar_nbh(i,j,1)),min(val_crisp(i,j,2),sar_nbh(i,j,2)));
        sar_ss(i,j)=min(sar_sa(i,j),sar_sb(i,j));
        
        ctex_sa(i,j)=max(min(val_nbh(i,j,1),ctex_crisp(i,j,1)),min(val_nbh(i,j,2),ctex_crisp(i,j,2)));
        ctex_sb(i,j)=max(min(val_crisp(i,j,1),ctex_nbh(i,j,1)),min(val_crisp(i,j,2),ctex_nbh(i,j,2)));
        ctex_ss(i,j)=min(ctex_sa(i,j),ctex_sb(i,j));
        
        optex_sa(i,j)=max(min(val_nbh(i,j,1),optex_crisp(i,j,1)),min(val_nbh(i,j,2),optex_crisp(i,j,2)));
        optex_sb(i,j)=max(min(val_crisp(i,j,1),optex_nbh(i,j,1)),min(val_crisp(i,j,2),optex_nbh(i,j,2)));
        optex_ss(i,j)=min(optex_sa(i,j),optex_sb(i,j));
    end
end
   
        
% figure; imshow(sar_fcr);
% figure; imshow(ctex_fcr);
% figure; imshow(optex_fcr);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmsar3_fcr.asc',sar_fcr);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmctex3_fcr.asc',ctex_fcr);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmoptex3_fcr.asc',optex_fcr);

% figure; imshow(sar_ss); colorbar;
% figure; imshow(ctex_ss); colorbar;
% figure; imshow(optex_ss); colorbar;
% 
% info=geotiffinfo('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmval2.tif');
% [~,R]=geotiffread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmval2.tif');
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmsar3_ss.tif', ...
%      sar_ss, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmctex3_ss.tif', ...
%      ctex_ss, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmoptex3_ss.tif', ...
%      optex_ss, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmsar3_fcr.tif', ...
%      sar_fcr, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmctex3_fcr.tif', ...
%      ctex_fcr, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmoptex3_fcr.tif', ...
%      optex_fcr, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
 
%  
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilsar3_ss.asc',sar_ss);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex3_ss.asc',ctex_ss);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhiloptex3_ss.asc',optex_ss);



%% fuzzy kappa calculation
p_val=zeros(2,1);
p_sar=zeros(2,1);
p_ctex=zeros(2,1);
p_optex=zeros(2,1);

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

p_flood_val=zeros(2,1);
p_flood_sar=zeros(2,1);
p_flood_ctex=zeros(2,1);
p_flood_optex=zeros(2,1);
for k=1:2
    p_flood_val(k)=sum(reshape(c_delta_val(:,:,k),m*n,1))/(m*n);
    p_flood_sar(k)=sum(reshape(c_delta_sar(:,:,k),m*n,1))/(m*n);
    p_flood_ctex(k)=sum(reshape(c_delta_ctex(:,:,k),m*n,1))/(m*n);
    p_flood_optex(k)=sum(reshape(c_delta_optex(:,:,k),m*n,1))/(m*n);
end
% %% pixel wise error stats
% 
% x=val-sar;
% dimx=size(x);
% y=val-ctex;
% dimy=size(y);
% z=val-optex;
% dimz=size(z);
% 
% mae_sar=mae(x);
% mae_ctex=mae(y);
% mae_optex=mae(z);
% rmse_sar=sqrt(sumsqr(x)/(dimx(1)*dimx(2)));
% rmse_ctex=sqrt(sumsqr(y)/(dimy(1)*dimy(2)));
% rmse_optex=sqrt(sumsqr(z)/(dimz(1)*dimz(2)));
% 
% display(rmse_sar);
% display(rmse_ctex);
% display(rmse_optex);
% 
% figure; imshow(abs(x)); colormap; colorbar;
% figure; imshow(abs(y)); colormap; colorbar;
% figure; imshow(abs(z)); colormap; colorbar;
% 
% %% correlation
% 
% val=GRID_CODE;
% x_corr=corr(val,sar);
% y_corr=corr(val,stdtex);
% z_corr=corr(val,optex);
% valsar=zeros(size(val));
% valstdtex=zeros(size(val));
% valoptex=zeros(size(val));
% for i=1:size(sar,1)
%     if val(i)==1
%             if val(i)==sar(i)
%                 valsar(i)=[];
%             else
%                 valsar(i)=sar(i);
%             end
%             if val(i)==stdtex(i)
%                 valstdtex(i)=[];
%             else
%                 valstdtex(i)=stdtex(i);
%             end
%             if val(i)==optex(i)
%                 valoptex(i)=[];
%             else
%                 valoptex(i)=optex(i);
%             end
%     end
% end
% for i=1:size(valsar,1)
%     if val(i)==0
%             if val(i)==sar(i)
%                 valsar(i)=[];
%             else
%                 valsar(i)=sar(i);
%             end
%     end
% end
% for i=1:size(valstdtex,1)
%     if val(i)==0
%             if val(i)==stdtex(i)
%                 valstdtex(i)=[];
%             else
%                 valstdtex(i)=stdtex(i);
%             end
%     end
% end
% for i=1:size(valoptex,1)
%     if val(i)==0
%             if val(i)==optex(i)
%                 valoptex(i)=[];
%             else
%                 valoptex(i)=optex(i);
%             end
%     end
% end
%     
% %% ecdf plots
% figure; cdfplot(valsar); hold on; cdfplot(sar); hold on; cdfplot(stdtex);
% hold on; cdfplot(optex);
% 
% figure; cdfplot(val); hold on; cdfplot(sar); hold on; cdfplot(stdtex);
% hold on; cdfplot(optex);