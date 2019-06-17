close all; clearvars; clc;

%% INPUTS

% NOTE: both images in the same projection and in geotiff format

% folder
fldr='DATA/BAR/20170808';

% band to process {'red','green','blue','nir','thermal'}; use exact word
band='red';

% for processing 
hr_band_order=1; %aggieair tpically is rgbnir, so red is 1, green 2...
lr_band_order=4; %red in L8 is band 4,green 3, blue 2, nir, 5

%% PROCESSING
 mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);

% subfolders:
fldr_hr=dir([newdir,'\',fldr,'\Aggie*']); %high resolution e.g AggieAir
fldr_lr=dir([newdir,'\',fldr,'\L*']);      %low resolution e.g Landsat

% low resolution image
image_lr=dir([fldr_lr.folder,'\',fldr_lr.name,'\*B',num2str(lr_band_order),'*.tif']);
% high resolution image
image_hr=dir([fldr_hr.folder,'\',fldr_hr.name,'\*NIR*.tif']);

% spatial info from both 
info_lr=geotiffinfo([image_lr.folder,'\',image_lr.name]);
info_hr=geotiffinfo([image_hr.folder,'\',image_hr.name]);

% pixel centers for low resolution image
[x_lr,y_lr] = pixcenters(info_lr.RefMatrix, [info_lr.SpatialRef.RasterSize],'makegrid');
x_lr=single(x_lr);y_lr=single(y_lr);

% pixel centers for high resolution image
[x_hr,y_hr] = pixcenters(info_hr.RefMatrix, [info_hr.SpatialRef.RasterSize],'makegrid');
x_hr=single(x_hr); y_hr=single(y_hr);

% finding high resolution image pixels within coarse resolution image
[rows,cols] = map2pix(info_lr.SpatialRef,x_hr,y_hr);
rows=round(rows); cols=round(cols);rows=unique(rows);cols=unique(cols);
rows=rows(rows>0); cols=cols(cols>0);

% low resolution pixels coordinates that overlap with high resolution image
x_lr=x_lr(rows,cols);y_lr=y_lr(rows,cols);

% reading high resolution image band
hires=imread([image_hr.folder,'\',image_hr.name]);
hires=hires(:,:,hr_band_order);
    Xmax=size(hires,1); Ymax=size(hires,2);

% load psf information for Landsat 8 
values=xlsread('Landsat_LCDM_OLI_MSF.xlsx');
psf_range_lr=size(values,1); %in meters
pixelsize_hr=info_hr.PixelScale(1);
psf_for_hr=fix(psf_range_lr/pixelsize_hr);

psf_AT = values(:,lr_band_order+1);
psf_XT = values(:,lr_band_order+1+7);

PSF=psf_AT*psf_XT';

PSF = imresize(PSF,(psf_for_hr+1)/psf_range_lr);
PSF = PSF/sum(PSF(:));

im_lr=zeros(size(x_lr),'single');
var_lr=im_lr;

for i=1:numel(x_lr)
    [row_i,col_i] = map2pix(info_hr.RefMatrix,x_lr(i),y_lr(i));
    row_i=round(row_i); col_i=round(col_i);
    
    xg= row_i-psf_for_hr/2:row_i+psf_for_hr/2;
    yg= col_i-psf_for_hr/2:col_i+psf_for_hr/2;
    
    [X,Y] = meshgrid(xg,yg);
    
    dd=X>0&Y>0 &X<Xmax+1 &Y<Ymax+1;
    mrow=X(dd); mcol=Y(dd);
    linearInd = sub2ind(size(hires), mrow, mcol);
    linearInd=sort(linearInd);
    
    I=hires(linearInd);
    I(I<=0)=NaN;
    psf_i=PSF(dd);
    
    im_lr(i)=nansum(I.*psf_i);
    var_lr(i)=sum(~isnan(I))/numel(PSF);
end
im_lr(im_lr<=0)=NaN;
% if size(ff,2)==1,ff=reshape(ff,size(x));end
figure, 
subplot(1,2,1),imagesc(im_lr); axis image off; colorbar;
subplot(1,2,2),imagesc(var_lr); axis image off; colorbar;
