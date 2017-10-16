%% read image

% image_info=read_envihdr(...
% 'E:\Clarence\Floodmap_journal\Resampled\3m-tex_topo\StdSAR_hand.hdr');
% image=multibandread(...
% 'E:\Clarence\Floodmap_journal\Resampled\3m-tex_topo\StdSAR_hand',...
% image_info.size,image_info.format,image_info.header_offset,....
% image_info.interleave,image_info.machine);
image=imread('E:\Clarence\Floodmap_journal\Resampled\3m-texture\ctex3.tif');
[r,c,b]=size(image);
inputdata=double(reshape(image,r*c,b));

%% SAR data training %%
info_flood=read_envihdr(...
'E:\Clarence\Floodmap_journal\Resampled\3m-texture\ctexfloodtrain3.hdr');
flood_training=multibandread(...
'E:\Clarence\Floodmap_journal\Resampled\3m-texture\ctexfloodtrain3',...
info_flood.size,info_flood.format,info_flood.header_offset,....
info_flood.interleave,info_flood.machine);
info_nonflood=read_envihdr(...
'E:\Clarence\Floodmap_journal\Resampled\3m-texture\ctexnonfloodtrain3.hdr');
nonflood_training=multibandread(...
'E:\Clarence\Floodmap_journal\Resampled\3m-texture\ctexnonfloodtrain3',...
info_nonflood.size,info_nonflood.format,info_nonflood.header_offset,....
info_nonflood.interleave,info_nonflood.machine);
[imax,jmax,k]=size(flood_training);
[mmax,nmax,p]=size(nonflood_training);

%% reshape data as compatible training inputs 
cnt = 1;
for i = 1:imax
    disp(i)
for j = 1:jmax
    data(1,:)=flood_training(i,j,:);
if sum(data)>0
data11(cnt,:) = flood_training(i,j,:);  cnt=cnt+1;
end
end
end
cnt = 1;
for i = 1:mmax
disp(i)
for j = i:nmax
    data(1,:)=nonflood_training(i,j,:);
if sum(data)>0
data12(cnt,:) = nonflood_training(i,j,:);  cnt=cnt+1;
end
end
end
dims=size(data11);
dim=size(data12);
flood(1:dims(1),1)=1;
nonflood(1:dim(1),1)=0;
response=vertcat(flood,nonflood);
trainingdata=horzcat(vertcat(data11,data12),response);
[rows,cols]=size(response);

%% divide into training, validation and testing

[trainsar_data,valsar_data,testsar_data] = dividerand(rows,0.7,0.15,0.15);
train=trainingdata(trainsar_data,:);
val=trainingdata(valsar_data,:);
test=trainingdata(testsar_data,:);
%% assess input data
results=evalfis(inputdata,fuzzydemctex);
out=reshape(results,r,c);
figure; imshow(out); colormap; colorbar;

%% save outputs

info=geotiffinfo('E:\Clarence\Floodmap_journal\Resampled\3m-texture\ctex3.tif');
[~,R]=geotiffread('E:\Clarence\Floodmap_journal\Resampled\3m-texture\ctex3.tif');
geotiffwrite('E:\Clarence\Floodmap_journal\Resampled\3m-texture\fuzzystdem2.tif', ...
    out, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

