val=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilrel3.asc');
sar=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilsar3.asc');
ctex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhilctex3.asc');
optex=arcgridread('F:\Clarence\Floodmap_journal\Resampled\3m-texture\jhiloptex3.asc');

x=val-sar(:,1:2276);
dimx=size(x);
y=val-ctex(:,1:2276);
dimy=size(y);
z=val-optex(:,1:2276);
dimz=size(z);

rmse_sar=sqrt(sumsqr(x)/(dimx(1)*dimx(2)));
rmse_ctex=sqrt(sumsqr(y)/(dimy(1)*dimy(2)));
rmse_optex=sqrt(sumsqr(z)/(dimz(1)*dimz(2)));

display(rmse_sar);
display(rmse_ctex);
display(rmse_optex);

