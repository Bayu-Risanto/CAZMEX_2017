clear all
close all
clc

% This script is to the RMSE between GPS gauges and satellites GPS
%  station pixels. Using values and normalized values 
% OUTPUT IS SCATTER PLOT OF STATION ON A MAP !!!

%% load GPS and satellites
load '/ha1/bayu/data/GPS_2017/LatLon_V1_2017.mat'
load '/ha1/bayu/data/SATELLITE_VALIDATION/GPS/extract/sta_hourly_GPS_V1.mat';
load '/ha1/bayu/data/SATELLITE_VALIDATION/GPME/extract/sta_hourly_PrGPME_V1.dat';
load '/ha1/bayu/data/SATELLITE_VALIDATION/GPMF/extract/sta_hourly_PrGPMF_V1.dat';
load '/ha1/bayu/data/SATELLITE_VALIDATION/PERSIANN/extract/sta_hourly_PrPERS_V1.dat';
load '/ha1/bayu/data/SATELLITE_VALIDATION/CMORPH/extract/sta_hourly_PrCMOR_V1.dat';

load 'NAMcoast.mat'
HGT = double(ncread('/ha1/jmoker/bayu/HGT_d01.nc','HGT'));
Xlo = double(ncread('/ha1/jmoker/bayu/HGT_d01.nc','XLONG'));
Xla = double(ncread('/ha1/jmoker/bayu/HGT_d01.nc','XLAT'));
lat = double(ncread('/ha1/bayu/WRF/NAM_runs/Jul27_2017/determ_no_RAP/2D_variables_d03.nc','XLAT'));
lon = double(ncread('/ha1/bayu/WRF/NAM_runs/Jul27_2017/determ_no_RAP/2D_variables_d03.nc','XLONG'));

list = {'GPM-E','GPM-F','CMORPH','PERSIANN'};

layer = Jun2Sep(:,:,1);

%% process
precip = Jun2Sep(:,6,1:25);
sta_hourly_GPS_V1 = permute(precip,[3 1 2]);

%% calculate RMSE. Use normalized values
mxE = max(max(vertcat(sta_hourly_GPS_V1,sta_hourly_PrGPME_V1(:,2:end))));
mxF = max(max(vertcat(sta_hourly_GPS_V1,sta_hourly_PrGPMF_V1(:,2:end))));
mxC = max(max(vertcat(sta_hourly_GPS_V1,sta_hourly_PrCMOR_V1)));

% set length without NaN for each station
for i = 1:length(Jun2Sep(1,1,:))
    Length(i) = length(sta_hourly_GPS_V1(i,:)) - sum(isnan(sta_hourly_GPS_V1(i,:)));
end 

% calculate RMSE
for i = 1:length(sta_hourly_PrGPME_V1(:,1))
    
    resGPME(i,:) = nansum( ( sta_hourly_PrGPME_V1(i,2:end) - sta_hourly_GPS_V1(i,:) ).^ 2 );
    RMSE(1,i) = sqrt( resGPME(i,:) / ( Length(i) - 1 ) );
    
    resGPMF(i,:) = nansum( ( sta_hourly_PrGPMF_V1(i,2:end) - sta_hourly_GPS_V1(i,:) ).^ 2 );
    RMSE(2,i) = sqrt( resGPMF(i,:) / ( Length(i) - 1 ) );
    
    resCMOR(i,:) = nansum( ( sta_hourly_PrCMOR_V1(i,:) - sta_hourly_GPS_V1(i,:) ).^ 2 );
    RMSE(3,i) = sqrt( resCMOR(i,:) / ( Length(i) - 1 ) );
    
end 
%% calculate for PERSIANN only
% we need to cut 02 August because PERSIANN misses this date.
% since the date is always in the same index, we need only one sample
for i = 1
idx = find(Jun2Sep(:,2,i) < 214); 
jdx = find(Jun2Sep(:,2,i) >= 215);
end

% now cut
sta_hourly_GPS2PERS = permute(precip(vertcat(idx,jdx),:,:), [3 1 2]);

for i = 1:length(sta_hourly_PrPERS_V1(:,1))
    
    resPERS(i,:) = nansum( ( sta_hourly_PrPERS_V1(i,:) - sta_hourly_GPS2PERS(i,:) ).^ 2 );
    RMSE(4,i) = sqrt( resPERS(i,:) / ( Length(i) - 1 ) );

end

meanRMSE = nanmean(RMSE,2);

% RMSE_inv = RMSE';
% bb = find(RMSE_inv == 0);
% RMSE_inv(bb) = NaN;
longitude(15,1) = NaN;  %remove SA08 because it seems to have extreme precip 
latitude(15,1) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalized the RMSE
denom = nanmean(sta_hourly_GPS_V1,2);

for i = 1:length(sta_hourly_GPS_V1(:,1))
%     denomE(i) = max( sta_hourly_PrGPME(i,:) )  - min( sta_hourly_PrGPME(i,:) );
    denomE = nanmean(sta_hourly_PrGPME_V1,2);
    normalized(1,i) = RMSE(1,i)./denom(i); %denomE(i);
    
%     denomF(i) = max( sta_hourly_PrGPMF(i,:) )  - min( sta_hourly_PrGPMF(i,:) );
    denomF = nanmean(sta_hourly_PrGPMF_V1,2);    
    normalized(2,i) = RMSE(2,i)./denom(i); %denomF(i);
    
%     denomC(i) = max( sta_hourly_PrCMOR(i,:) )  - min( sta_hourly_PrCMOR(i,:) );
    denomC = nanmean(sta_hourly_PrCMOR_V1,2);
    normalized(3,i) = RMSE(3,i)./denom(i); %denomC(i);
end 
for i = 1:length(sta_hourly_GPS_V1(:,1))
%     denomP(i) = max( sta_hourly_PrPERS(i,:) ) - min( sta_hourly_PrPERS(i,:) );
    denomP = nanmean(sta_hourly_PrPERS_V1,2);
    normalized(4,i) = RMSE(4,i)./denom(i); %denomP(i);
end

aa = find(normalized == Inf); normalized(aa) = NaN;
meanNRMSE = nanmean(normalized,2);
num_nan = isnan(normalized); cc = find(num_nan(1,:) == 1);
longitude(cc) = NaN; latitude(cc) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
cmap2 = jet(40);

fig = figure(1);

% set(gcf,'Position',[20 20 1500 650])
set(gcf,'Position',[20 20 1300 650])
% ha = tight_subplot(2,4,[.006 .006],[.16 .05],[.03 .02]);
ha = tight_subplot(1,4,[.007 .006],[.1 .08],[.03 .02]);

ct = 0;

for i = 1:4
    
    ct = ct + 1;
    set(fig,'CurrentAxes',ha(ct))
    
    toplot1 = normalized(i,:);    
    
    sc = scatter(longitude,latitude,50,toplot1,'filled');
    sc.MarkerEdgeColor = 'k';
    shading flat; 
%     caxis([0 4])
    caxis([0 20])
    colormap(ha(ct),cmap2)
    
                hold on;
    
    contour(lon,lat,NAMcoast','LevelList',1.5,...
        'LineColor','black','LineWidth',1.25);
    contour(Xlo,Xla,HGT,'LevelList',(500:500:2500),...
        'LineColor','black','LineWidth',1);
    
    xlim([-115 -105])
    ylim([22.5 34])
    
    set(gca,'layer','top')
    
    set(gca,'DataAspectRatio',[5 4.5 1])
    set(gca,'YAxisLocation','right')
    
    set(gca,'xtick',(-114:2:-106))
    set(gca,'ytick',(22:2:32))
    set(gca,'XTickLabelRotation',45, 'TickLabelInterpreter','latex')
    set(gca,'TickLength',[0.03, 0.01])
    set(gca, 'FontSize', 12)
    
    title(list{i},'FontSize',12)
end

colorbar('location','southoutside','position',[.05 0.08 .85 0.015], ...
   'FontSize',10,'FontWeight','bold');

%% Set super-title position, then save it in .png
axes( 'Position', [0, 0.9, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'NRMSE (GPS vs Satellites) based on precip (mm) from 30 June to 12 September', 'FontSize', 12', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top' ) ;