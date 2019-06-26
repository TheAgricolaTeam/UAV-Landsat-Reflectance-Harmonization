% CODE TO HARMONIZE UAV AND LANDSAT 8 BANDS INFORMATION
% USING LANDSAT 8 POINT SPREAD FUNCTION INFORMATION
% By Alfonso Torres, 2016-2019 alfonso.torres@usu.edu

% This code produces a Landsat derived image from UAV and then produces a
% scatter comparison from Landsat derived and true Landsat.
% The code also provides the parameters to "correct" UAV band towards
% Landsat radiometric response. The result values are a slope and bias that
% can be used in ArcGIS Raster calculator or python code

% Asumptions:
% UAV and Landsat information were collected the same date/time (within
% 20min window)
% Both, UAV and Landsat are using the same projection system (WGS64 UTM) so
% this means the images are geotiff
% It is assumed that user is using Landsat Surface Reflectance product,
% also know as "Landsat Collection 1 Level 2"
% UAV/Landsat data may be scaled (by 10,000). Make sure this value is
% known.

% Code take care of pixel locations in UAV and Landsat bands.
% Code does not know which bands are being processed. Know beforehand the
% order or band names in both UAV and Landsat
% This code does not process THERMAL. Thermal requires a different approach.





close all; clearvars; clc;

%% INPUTS

% NOTE: both images in the same projection and in geotiff format

% folder where both UAV and Landsat are located
fldr='DATA/BAR/20170808';

% band to process {'red','green','blue','nir'}; use exact word
% band='red';

% Landsat scaling value (from Landsat metadata/handbook)
landsat_scale=10000;
% HIgh resolution scale value
uav_scale=1; %one means no scale (so data goes from zero to one)

% for processing the band name
uav_band_order=4; %UAV multiband tpically is rgbnir, so red is 1, green 2, blue is 3, nir iis 4
landsat_band_order=5; %red in L8 is band 4,green 3, blue 2, nir, 5

%% PROCESSING
% locating images in directory
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end-1)-1);

% subfolders:
fldr_uav=dir([newdir,'\',fldr,'\Aggie*']); %high resolution file name e.g AggieAir
fldr_landsat=dir([newdir,'\',fldr,'\L*SR']);   %low resolution file name e.g Landsat

% locating low resolution image
image_landsat=dir([fldr_landsat.folder,'\',fldr_landsat.name,'\*B*',num2str(landsat_band_order),'*.tif']);
% locating high resolution image
image_uav=dir([fldr_uav.folder,'\',fldr_uav.name,'\*NIR*.tif']);

% spatial info from both images
info_landsat=geotiffinfo([image_landsat.folder,'\',image_landsat.name]);
info_uav=geotiffinfo([image_uav.folder,'\',image_uav.name]);
% reading landsat image
img_landsat=imread([image_landsat.folder,'\',image_landsat.name]);
% Landsat SR has negative values, thus eliminating them
img_landsat(img_landsat<=0)=NaN;

% pixel centers for low resolution image
[x_landsat,y_landsat] = pixcenters(info_landsat.RefMatrix, [info_landsat.SpatialRef.RasterSize],'makegrid');
x_landsat=single(x_landsat);
y_landsat=single(y_landsat);

% pixel centers for high resolution image
[x_uav,y_uav] = pixcenters(info_uav.RefMatrix, [info_uav.SpatialRef.RasterSize],'makegrid');
x_uav=single(x_uav); 
y_uav=single(y_uav);

% finding high resolution image pixels within coarse resolution image
[rows,cols] = map2pix(info_landsat.SpatialRef,x_uav,y_uav);
rows=round(rows); 
cols=round(cols);
rows=unique(rows);
cols=unique(cols);
% selecting landsat pixels that ovelap UAV area
img_landsat=img_landsat(rows,cols);

% low resolution pixels coordinates that overlap with high resolution image
x_landsat=x_landsat(rows,cols);
y_landsat=y_landsat(rows,cols);

% reading high resolution image band
uav=imread([image_uav.folder,'\',image_uav.name]);
uav=uav(:,:,uav_band_order)/uav_scale;
Xmax=size(uav,1); Ymax=size(uav,2);

% load psf information for Landsat 8
values=xlsread('Landsat_LCDM_OLI_MSF.xlsx');
psf_range_landsat=size(values,1); %in meters
pixelsize_uav=info_uav.PixelScale(1); % pixel size of UAV band
psf_for_uav=fix(psf_range_landsat/pixelsize_uav); %number of UAV pixels within PSF 

psf_AT = values(:,landsat_band_order+1); %PSF along Landsat track direction 
psf_XT = values(:,landsat_band_order+1+7); % PSF across Landsat track direction

PSF=psf_AT*psf_XT';
PSF(PSF<0)=0;

PSF = imresize(PSF,(psf_for_uav+1)/psf_range_landsat);
PSF = PSF/sum(PSF(:)); % ensuring 2D PSF sums one

uav_landsat=zeros(size(x_landsat),'single'); %empty matrix for derived landsat from uav
var_landsat=uav_landsat; % matrix uav pixel count ratio within a PSF (quality measurement)

parfor i=1:numel(x_landsat)
%     finding the location in UAV image that correspond to landsat pixel
%     center
    [row_i,col_i] = map2pix(info_uav.RefMatrix,x_landsat(i),y_landsat(i));
    row_i=round(row_i); col_i=round(col_i);
%     finding UAV pixel locations for corresponding PSF
    xg= row_i-psf_for_uav/2:row_i+psf_for_uav/2;
    yg= col_i-psf_for_uav/2:col_i+psf_for_uav/2;
    
    [X,Y] = meshgrid(xg,yg); % preparing UAV pixel locations for matrix operation
    
    X=reshape(X,[],1); %matrix to colunm
    Y=reshape(Y,[],1);%matrix to colunm
    psf_i=reshape(PSF,[],1);%matrix to colunm
    
    dd=X>0&Y>0 &X<Xmax+1 &Y<Ymax+1; %logical filter for image edges
    mrow=X(dd); 
    mcol=Y(dd);
    psfa=psf_i(dd);
    
    linearInd = sub2ind(size(uav), mrow, mcol);
    
    I=uav(linearInd);
    I(I<=0)=NaN; %empty UAV pixels forced to NaN
    
    uav_landsat(i)=nansum(I.*psfa); %uav psf product per pixel
    var_landsat(i)=sum(~isnan(I))/numel(PSF); 
end
uav_landsat(uav_landsat<=0)=NaN;

img_landsat=single(img_landsat);
img_landsat(isnan(uav_landsat))=NaN;
img_landsat=img_landsat/landsat_scale;

% plot
uav_landsat(var_landsat<1)=NaN;
img_landsat(var_landsat<1)=NaN;

% % if size(ff,2)==1,ff=reshape(ff,size(x));end
% figure,
% subplot(1,3,1),imagesc(im_landsat); axis image off; colorbar;
% subplot(1,3,2),imagesc(var_landsat); axis image off; colorbar;
% subplot(1,3,3),imagesc(landsat_im); axis image off; colorbar;
%% regression
mdl = fitlm(uav_landsat(:),img_landsat(:),'RobustOpts','on');
disp(mdl);
%% results
result{1,1}='Landsat file name';
result{1,2}= 'band number';
result{1,3}='UAV file name';
result{1,4}='band number';
result{1,5}='correction slope';
result{1,6}='correction bias';
result{1,7}='approach';
result{2,1}=image_landsat.name;
result{2,2}=landsat_band_order;
result{2,3}=image_uav.name;
result{2,4}=uav_band_order;
result{2,7}='Point Spread Function';

figure
subplot(1,2,1),plot(uav_landsat(:),img_landsat(:),'.'); grid on;
xlim([0 0.5]);ylim([0 0.5]);
refline(1,0);
axis square;
title('UAV - Landsat before comparison');
ylabel('Landsat reflectance');
xlabel('UAV reflectance');
coeff=mdl.Coefficients.Estimate;
subplot(1,2,2),plot(uav_landsat(:)*coeff(2)+coeff(1),img_landsat(:),'.'); grid on;
xlim([0 0.5]);ylim([0 0.5]);
refline(1,0)
axis square
title('UAV - Landsat after comparison');
ylabel('Landsat reflectance');
xlabel('UAV reflectance');
sgtitle('PSF results')

result{2,5}=coeff(2);
result{2,6}=coeff(1);


disp('double click on result variable')
