% Read files

val=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhil3.asc');
sar=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilsar3.asc');
ctex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex3.asc');
optex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhiloptex3.asc');

% sar=sar(:,3:2278);
% ctex=ctex(:,3:2278);
% optex=optex(:,3:2278);

% create fuzzy vectors
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
%  preallocate vectors
val_crisp=zeros(m,n,2);
sar_crisp=zeros(m,n,2);
ctex_crisp=zeros(m,n,2);
optex_crisp=zeros(m,n,2);
val_nbh=zeros(m,n,2);
sar_nbh=zeros(m,n,2);
ctex_nbh=zeros(m,n,2);
optex_nbh=zeros(m,n,2);

sar_sa=zeros(m,n,1);
sar_sb=zeros(m,n,1);
sar_ss=zeros(m,n,1);

ctex_sa=zeros(m,n,1);
ctex_sb=zeros(m,n,1);
ctex_ss=zeros(m,n,1);

optex_sa=zeros(m,n,1);
optex_sb=zeros(m,n,1);
optex_ss=zeros(m,n,1);

% fuzzy map comparison
sar_fcr=zeros(m,n);
ctex_fcr=zeros(m,n);
optex_fcr=zeros(m,n);

%nbh=2^(-d/2)---> exp. decay

%%Hagen 2003

for j=1:n
    for i=1:m
        % calculation of cell-by-cell similarity
        sar_fcr(i,j)=fuzzy_centralcell(val_kf,sar_kf,i,j);
        ctex_fcr(i,j)=fuzzy_centralcell(val_kf,ctex_kf,i,j);
        optex_fcr(i,j)=fuzzy_centralcell(val_kf,ctex_kf,i,j);
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
% calculation of neighbourhood influences      
for j=1:n-2
    for i=1:m-2
        for k=1:2
        A=nbh(val_kf,i,j,k);
        val_nbh(i+1,j+1,k)=max(A);
        B=nbh(sar_kf,i,j,k);    
        sar_nbh(i+1,j+1,k)=max(B);
        C=nbh(ctex_kf,i,j,k);
        ctex_nbh(i+1,j+1,k)=max(C);
        D=nbh(optex_kf,i,j,k);
        optex_nbh(i+1,j+1,k)=max(D);
        end
    end
end

for j=1:n
    for i=1:m
        % Calculation of spatial similarity
        sar_ss(i,j)=nbh_cc(val_nbh,sar_crisp,val_crisp,sar_nbh,i,j);
        ctex_ss(i,j)=nbh_cc(val_nbh,ctex_crisp,val_crisp,ctex_nbh,i,j);
        optex_ss(i,j)=nbh_cc(val_nbh,optex_crisp,val_crisp,optex_nbh,i,j);
    end
end
    
%         
% figure; imshow(sar_fcr);
% figure; imshow(ctex_fcr);
% figure; imshow(optex_fcr);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmsar3_fcr.asc',sar_fcr);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmctex3_fcr.asc',ctex_fcr);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\ulmoptex3_fcr.asc',optex_fcr);
% 
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
%  
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilsar3_ss.asc',sar_ss);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex3_ss.asc',ctex_ss);
% dlmwrite('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhiloptex3_ss.asc',optex_ss);



% fuzzy kappa calculation

ring_num=[0,1,2,3,4,5,6,7,8,9];
num_cells=[1,4,4,4,8,4,4,8,8,4];
cum_num_cells=[0,4,8,12,20,24,28,36,44,48];
dm=[1,0.71,0.61,0.5,0.46,0.38,0.35,0.33,0.30,0.25];
p=zeros(10);
C=zeros(10);
Ec=zeros(10);
for j=10:n-10
    for i=10:m-10
        if val(i,j)==sar(i,j)
            delta=1;
        else
            delta=0;
        end
        for ii=1:size(ring_num)-1
            p(cum_num_cells(ii))=((1-(1-val(i,j)))^cum_num_cells(ii))*((1-(1-sar(i,j)))^cum_num_cells(i,j));
            C(ii)=((1-delta)*val(i,j)*sar(i,j)*(p(ii+1)-p(ii)));
            Ec(ii)=val(i,j)*sar(i,j);
            
        end
    
    end
    
end


% pixel wise error stats

x=val-sar;
dimx=size(x);
y=val-ctex;
dimy=size(y);
z=val-optex;
dimz=size(z);

mae_sar=mae(x);
mae_ctex=mae(y);
mae_optex=mae(z);
rmse_sar=sqrt(sumsqr(x)/(dimx(1)*dimx(2)));
rmse_ctex=sqrt(sumsqr(y)/(dimy(1)*dimy(2)));
rmse_optex=sqrt(sumsqr(z)/(dimz(1)*dimz(2)));

display(rmse_sar);
display(rmse_ctex);
display(rmse_optex);

figure; imshow(abs(x)); colormap; colorbar;
figure; imshow(abs(y)); colormap; colorbar;
figure; imshow(abs(z)); colormap; colorbar;

% correlation

val=GRID_CODE;
x_corr=corr(val,sar);
y_corr=corr(val,stdtex);
z_corr=corr(val,optex);
valsar=zeros(size(val));
valstdtex=zeros(size(val));
valoptex=zeros(size(val));
for i=1:size(sar,1)
    if val(i)==1
            if val(i)==sar(i)
                valsar(i)=[];
            else
                valsar(i)=sar(i);
            end
            if val(i)==stdtex(i)
                valstdtex(i)=[];
            else
                valstdtex(i)=stdtex(i);
            end
            if val(i)==optex(i)
                valoptex(i)=[];
            else
                valoptex(i)=optex(i);
            end
    end
end
for i=1:size(valsar,1)
    if val(i)==0
            if val(i)==sar(i)
                valsar(i)=[];
            else
                valsar(i)=sar(i);
            end
    end
end
for i=1:size(valstdtex,1)
    if val(i)==0
            if val(i)==stdtex(i)
                valstdtex(i)=[];
            else
                valstdtex(i)=stdtex(i);
            end
    end
end
for i=1:size(valoptex,1)
    if val(i)==0
            if val(i)==optex(i)
                valoptex(i)=[];
            else
                valoptex(i)=optex(i);
            end
    end
end
    
% ecdf plots
figure; cdfplot(valsar); hold on; cdfplot(sar); hold on; cdfplot(stdtex);
hold on; cdfplot(optex);

figure; cdfplot(val); hold on; cdfplot(sar); hold on; cdfplot(stdtex);
hold on; cdfplot(optex);

