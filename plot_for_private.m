%% Geomap
addpath('/Volumes/mydirve/cloud_annual/m_map');
m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
[CS,CH]=m_etopo2('contourf',[-12000:200:12000],'edgecolor','none');
topo_data=CH.ZData;
lon_topo=linspace(240,300,451);
lat_topo=linspace(30,0,451);
[lon_topo,lat_topo]=meshgrid(lon_topo,lat_topo);
close
figure('pos',[10 10 1500 1500]);
worldmap([0 30],[240 300])
geoshow(lat_topo,lon_topo,topo_data,'DisplayType','texturemap')
hold on
bordersm('countries','linewidth',1,'color','k');
setm(gca,'fontsize',16,'fontweight','bold')
colormap([ m_colmap('blues',80); m_colmap('gland',80)]);
%brighten(.5);
caxis([-8000 8000]);
x=scaleruler('units','km');
setm(handlem('scaleruler1'), ...
    'XLoc',-2e6, ...
    'YLoc',0.5e6, ...
    'TickDir','down', ...
    'MajorTick',0:200:1000,...
    'RulerStyle','patches',...
    'linewidth',2,...
    'fontsize',12)
textm(12,260,'Pacific','fontsize',16,'fontweight','bold');
s=northarrow('latitude',4,'longitude',245,'scaleratio',0.1,'linewidth',3,...
    'FaceColor', [1.000 0.8431 0.0000], ...
      'EdgeColor', [0.0100 0.0100 0.9000]);
 textm(12,260,'Pacific','fontsize',16,'fontweight','bold');
 textm(25.3043,180+180-90.0659,'Gulf of Mexico','fontsize',16,'fontweight','bold');
 textm(14.5401,180+180-74.9676,'Caribbean Sea','fontsize',16,'fontweight','bold');
labelbordersm('Mexico','fontsize',16,'fontweight','bold')
labelbordersm('Guatemala','fontsize',12,'fontweight','bold')
labelbordersm('Belize','fontsize',12,'fontweight','bold');
labelbordersm('Honduras','fontsize',12,'fontweight','bold');
labelbordersm('El Salvador','fontsize',12,'fontweight','bold');
labelbordersm('Nicaragua','fontsize',12,'fontweight','bold');
labelbordersm('Costa Rica','fontsize',12,'fontweight','bold');
labelbordersm('Panama','fontsize',12,'fontweight','bold');
labelbordersm('Cuba','fontsize',16,'fontweight','bold')
labelbordersm('Bahamas','fontsize',14,'fontweight','bold')
textm(20,288,'Dominican Republic','fontsize',14,'fontweight','bold');
labelbordersm('Jamaica','fontsize',10,'fontweight','bold')
labelbordersm('Colombia','fontsize',16,'fontweight','bold')
labelbordersm('Venezuela','fontsize',16,'fontweight','bold')
labelbordersm('Haiti','fontsize',14,'fontweight','bold')
%% Figure 1 mean states (frequency, start peak end dur Imsd)
load('msd_collection_ncep');
load('mean_char_msd');
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

full_freq=NaN(120,60);

for i=1:size(loc_unique,1);
    loc_here=loc_unique(i,:);
    msd_here=msd_collection_ncep(msd_collection_ncep(:,2)==loc_here(1) & msd_collection_ncep(:,3)==loc_here(2),:);
    full_freq(loc_here(1),loc_here(2))=size(msd_here,1)./39;
    
end

figure('pos',[10 10 1500 1500]);
m_proj('miller','lon',[nanmin(double(lon_ncep)) nanmax(double(lon_ncep))],'lat',[nanmin(double(lat_ncep)) nanmax(double(lat_ncep))]);
h=tight_subplot(3,2,[0.025,0.01],[0.03 0.01],[0.05,0.1]);

% freq
axes(h(1));
m_contourf(lon_ncep,lat_ncep,full_freq',0:0.01:1,'linestyle','none');
m_coast();
m_grid('xtick',[],'fontsize',16,'linestyle','none');
colormap(m_colmap('jet',11));
caxis([0 1]);
m_text(243,4,'a) Frequency','fontsize',16,'fontweight','bold')
s=colorbar('fontsize',16);
s.Label.String='/year';

% start
axes(h(2));
m_contourf(lon_ncep,lat_ncep,full_start',140:0.1:190,'linestyle','none');
m_coast();
m_grid('xtick',[],'ytick',[]);
colormap(m_colmap('jet',11));
caxis([140 190]);
m_text(243,4,'b) Onset date','fontsize',16,'fontweight','bold')
s=colorbar('ticks',[datenum(2016,6,1)-datenum(2016,1,1)+1 datenum(2016,7,1)-datenum(2016,1,1)+1],'ticklabels',{'June 1st','July 1st'},'fontsize',16);
s.Label.String='';

% peak
axes(h(3));
m_contourf(lon_ncep,lat_ncep,full_peak',170:0.1:240,'linestyle','none');
m_coast();
m_grid('xtick',[],'fontsize',16,'linestyle','none');
colormap(m_colmap('jet',11));
caxis([190 230]);
m_text(243,4,'c) Peak date','fontsize',16,'fontweight','bold')
s=colorbar('ticks',[datenum(2016,7,15)-datenum(2016,1,1)+1 datenum(2016,8,1)-datenum(2016,1,1)+1],'ticklabels',{'July 15th','August 1st'},'fontsize',16);
s.Label.String='';

% end
axes(h(4));
m_contourf(lon_ncep,lat_ncep,full_end',230:0.1:280,'linestyle','none');
m_coast();
m_grid('xtick',[],'ytick',[]);
colormap(m_colmap('jet',11));
caxis([230 280]);
m_text(243,4,'c) End date','fontsize',16,'fontweight','bold')
s=colorbar('ticks',[datenum(2016,9,1)-datenum(2016,1,1)+1 datenum(2016,10,1)-datenum(2016,1,1)+1],'ticklabels',{'Sep 1st','Oct 1st'},'fontsize',16);
s.Label.String='';

% dur
axes(h(5));
m_contourf(lon_ncep,lat_ncep,full_dur',50:0.1:150,'linestyle','none');
m_coast();
m_grid('fontsize',16,'linestyle','none');
colormap(m_colmap('jet',11));
caxis([50 130]);
m_text(243,4,'d) Duration','fontsize',16,'fontweight','bold')
s=colorbar('fontsize',16);
s.Label.String='days';

% imsd
axes(h(6));
m_contourf(lon_ncep,lat_ncep,imsd_climatology_ncep',0:0.01:0.6,'linestyle','none');
m_coast();
m_grid('ytick',[],'fontsize',16,'linestyle','none');
colormap(m_colmap('jet',11));
caxis([0 0.6]);
m_text(243,4,'e) I_{msd}','fontsize',16,'fontweight','bold')
s=colorbar('fontsize',16);
s.Label.String='';

axes(h(6));
box off

%% Figure 2 annual MSD number
load('msd_num');
figure('pos',[10 10 1000 1000]);
subplot(2,1,1);
plot(1979:2017,msd_num./2511,'linewidth',2);
xlabel('Year','fontsize',16);
ylabel('Frequency','fontsize',16);
set(gca,'fontsize',16);
xlim([1979 2017]);
text(1981,720,'a)','fontsize',16,'fontweight','bold');
subplot(2,1,2);
[acf,lags,bounds] = autocorr(msd_num);
lineHandles = stem(lags,acf,'filled','r-o');
set(lineHandles(1),'MarkerSize',4)
grid('on')
xlabel('Lag (Years)','fontsize',16,'fontweight','bold')
ylabel('Autocorrelation','fontsize',16,'fontweight','bold')
%title('Sample Autocorrelation Function','fontsize',16)
hold('on')
plot(linspace(0,20,1000),ones(1000,1)*bounds(1),'-b');
hold on
plot(linspace(0,20,1000),ones(1000,1)*bounds(2),'-b');
hold on
plot(linspace(0,20,1000),ones(1000,1)*0,'-k');
hold('off')
xlim([0 20]);
set(gca,'fontsize',16);
hold on
text(1,0.8,'b)','fontsize',16,'fontweight','bold');

%% Figure 3 atmospheric and oceanic conditions 
load('mean_start_end_condition');
load('msd_collection_ncep');
load('all_sst_noaa','*start','*end','*peak');
lon_era_low = linspace(nanmin(lon_ncep),nanmax(lon_ncep),121);
lat_era_low = linspace(nanmax(lat_ncep),nanmin(lat_ncep),61);
load('daily_sst_noaa','l*');
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');

figure('pos',[10 10 1500 1500]);

h=tight_subplot(4,3,[0.025,0.01],[0.03 0.03],[0.05,0.05]);

sst_anom_end(isinf(sst_anom_end))=nan;
sst_anom_start(isinf(sst_anom_start))=nan;
sst_anom_peak(isinf(sst_anom_peak))=nan;

color_here =  m_colmap('diverging',11);

color_rain = customcolormap(linspace(0,1,11), {'#523107','#523107','#bf812c','#e2c17e','#f3e9c4','#f6f4f4','#cae9e3','#81cdc1','#379692','#01665e','#003d2e'},11);
color_rain=color_rain(end:-1:1,:);

axes(h(1))
m_contourf(lon_ncep,lat_ncep,precip_anom_start',-4.2:0.01:4.2,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
m_text(243,4,'a) Precipitation','fontsize',16,'fontweight','bold');
colormap(h(1),color_rain);
caxis([-4.2 4.2]);
title('Onset','fontsize',16);
s=colorbar('location','eastoutside','fontsize',16);
s.Label.String='mm/day';

axes(h(4))
m_contourf(lon_era_low,lat_era_low,temp_anom_start',-1:0.001:1,'linestyle','none')
m_coast();
colormap(h(4),color_here);
m_grid('xtick',[],'ytick',[]);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_anom_start(1:4:end,1:4:end))',(v_anom_start(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_text(243,4,'d) 2m Temperature','fontsize',16,'fontweight','bold')
%title('Peak','fontsize',16);
caxis([-0.4 0.4]);   
s=colorbar('location','eastoutside','fontsize',16);
s.Label.String='^{o}C';

axes(h(7))
m_contourf(lon_era_low,lat_era_low,(pres_anom_start)',-0.65:0.001:0.65,'linestyle','none')
colormap(h(7),color_here);
m_grid('xtick',[],'ytick',[]);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_anom_start(1:4:end,1:4:end))',(v_anom_start(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_coast;
m_text(243,4,'g) Pressure','fontsize',16,'fontweight','bold')
%title('End','fontsize',16);
caxis([-0.6 0.6]);   
s=colorbar('location','eastoutside','fontsize',16,'ticks',-0.6:0.2:0.6);
s.Label.String='hpa';

axes(h(10))
m_contourf(lon_sst,lat_sst,(sst_anom_start)',-0.25:0.001:0.25,'linestyle','none')
m_grid('fontsize',16,'linestyle','none');
m_coast;
colormap(h(10),color_here);
m_text(243,4,'j) SST','fontsize',16,'fontweight','bold')
%title('SST','fontsize',16);
caxis([-0.2 0.2]);  
s=colorbar('location','eastoutside','fontsize',16);
s.Label.String='^{o}C';

axes(h(2))
m_contourf(lon_ncep,lat_ncep,(precip_anom_peak)',-4.2:0.01:4.2,'linestyle','none')
m_coast();
colormap(h(2),color_rain);
m_grid('xtick',[],'ytick',[]);
m_text(243,4,'b) Peak','fontsize',16,'fontweight','bold');
title('Peak','fontsize',16);
caxis([-4.2 4.2]);

axes(h(5))
m_contourf(lon_era_low,lat_era_low,(temp_anom_peak)',-1:0.001:1,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_anom_peak(1:4:end,1:4:end))',(v_anom_peak(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
colormap(h(5),color_here);
m_text(243,4,'e)','fontsize',16,'fontweight','bold')
caxis([-0.4 0.4]);   

axes(h(8))
m_contourf(lon_era_low,lat_era_low,(pres_anom_peak)',-0.65:0.001:0.65,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_anom_peak(1:4:end,1:4:end))',(v_anom_peak(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
colormap(h(8),color_here);
m_text(243,4,'h)','fontsize',16,'fontweight','bold')
caxis([-0.6 0.6]);

axes(h(11))
m_contourf(lon_sst,lat_sst,sst_anom_peak',-0.25:0.001:0.25,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
colormap(h(11),color_here);
m_text(243,4,'k)','fontsize',16,'fontweight','bold')
caxis([-0.2 0.2]);

axes(h(3))
m_contourf(lon_ncep,lat_ncep,precip_anom_end',-4.2:0.01:4.2,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
m_text(243,4,'c)','fontsize',16,'fontweight','bold');
colormap(h(3),color_rain);
caxis([-4.2 4.2]);
title('End','fontsize',16);


axes(h(6))
m_contourf(lon_era_low,lat_era_low,temp_anom_end',-1:0.001:1,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_anom_end(1:4:end,1:4:end))',(v_anom_end(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
colormap(h(6),color_here);
m_text(243,4,'f)','fontsize',16,'fontweight','bold')
caxis([-0.4 0.4]);   


axes(h(9))
m_contourf(lon_era_low,lat_era_low,(pres_anom_end)',-0.65:0.001:0.65,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_anom_end(1:4:end,1:4:end))',(v_anom_end(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
colormap(h(9),color_here);
m_text(243,4,'i)','fontsize',16,'fontweight','bold')
caxis([-0.6 0.6]);

axes(h(12))
m_contourf(lon_sst,lat_sst,sst_anom_end',-0.25:0.001:0.25,'linestyle','none')
m_coast();
m_grid('xtick',[],'ytick',[]);
m_coast;
m_text(243,4,'l)','fontsize',16,'fontweight','bold')
colormap(h(12),color_here);
caxis([-0.2 0.2]);



%% Figure 10
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
load('precip_clim');
precip_smooth_clim=smoothdata(precip_clim,3,'gaussian',31);
load('msd_collection_ncep');

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

precip_eof_clim_msd=NaN(366,987);

for i=1:size(loc_unique,1);
    loc_here=loc_unique(i,:);
    precip_eof_clim_msd(:,i)=squeeze(precip_smooth_clim(loc_here(1),loc_here(2),:));
    
end
F=detrend(precip_eof_clim_msd,0);
[EOFs,lambda,contribution,PCs]=EOF_analysis(F');

EOFs(:,1)=EOFs(:,1).*std(PCs(1,:));

PCs_sd=PCs(1,:)./std(PCs(1,:));

close all

save eof_associated PCs_sd EOFs

load('eof_associated');
EOF_1=EOFs(:,1);
load('msd_collection_ncep');

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

EOF_1_plot=NaN(120,60);

for i=1:size(loc_unique,1);
    loc_here=loc_unique(i,:);
    EOF_here=EOF_1(i);
    EOF_1_plot(loc_here(1),loc_here(2))=EOF_here;
end

rainy_period=find(PCs_sd>=0);

fill_x=[rainy_period rainy_period(end:-1:1)];

fill_y=[zeros(length(rainy_period),1);(PCs_sd(rainy_period(end:-1:1)))'];

fill_before_x=[1:(rainy_period(1)-1) (rainy_period(1)-1):-1:1];

fill_before_y=[PCs_sd(1:(rainy_period(1)-1)) zeros(1,length(1:(rainy_period(1)-1)))];

fill_after_x=[(rainy_period(end)+1):366 366:-1:(rainy_period(end)+1)];

fill_after_y=[PCs_sd((rainy_period(end)+1):366) zeros(1,length((rainy_period(end)+1):366))];


subplot(2,1,1);

m_contourf(lon_ncep,lat_ncep,EOF_1_plot',-0.1:0.01:6,'linestyle','none');
m_coast();
m_grid('fontsize',16);
colormap(m_colmap('jet',11));
caxis([0 6]);
title('EOF1 (77.71%)','fontsize',16);
s=colorbar('fontsize',16);
s.Label.String='mm/day';

subplot(2,1,2);

plot(1:366,PCs_sd,'linewidth',2);
hold on
fill(fill_x,fill_y,[1 0.5 0.5],'linestyle','none','FaceAlpha',0.5);
hold on
fill(fill_before_x,fill_before_y,[0.5 0.5 1],'linestyle','none','FaceAlpha',0.5);
hold on
fill(fill_after_x,fill_after_y,[0.5 0.5 1],'linestyle','none','FaceAlpha',0.5);
xlim([1 366]);
set(gca,'xtick',[datenum(2016,1,15):30:datenum(2016,12,31)]-datenum(2016,1,1)+1,'xticklabels',{'J','F',...
    'M','A','M','J','J','A','S','O','N','D'},'fontsize',16);
set(gca,'fontsize',16);
title('PC1','fontsize',16);
text(datenum(2016,7,15)-datenum(2016,1,1)+1,0.5,'Rainy','fontsize',16,'fontweight','bold');

%% Figure 11
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
load('msd_collection_ncep');
load('coef_r_full');
logic_msd=NaN(120,60);

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

for i=1:size(loc_full_data,1);
    logic_msd(loc_full_data(i,1),loc_full_data(i,2))=0;
end

for i=1:size(loc_unique,1);
    logic_msd(loc_unique(i,1),loc_unique(i,2))=1;
end

load('p_logic');

figure('pos',[10 10 1000 1000]);
h=subplot(2,1,2);
m_contourf(lon_ncep,lat_ncep,(coef_plot_full')*1e07,-4.7:0.01:4.7,'linestyle','none');
m_grid('fontsize',16);
m_coast();
colormap(h,m_colmap('diverging',11));
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
title('b) b4 \times 10^{7}','fontsize',16);
s=colorbar('fontsize',16);
s.Label.String='mm/day^{4}';
caxis([-4 4]);

h=subplot(2,1,1);
m_contourf(lon_ncep,lat_ncep,(r_plot_full'),0:0.01:1,'linestyle','none');
m_grid('fontsize',16);
m_coast();
colormap(h,m_colmap('jet',11));
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
title('a) adjusted R^{2}','fontsize',16);
s=colorbar('fontsize',16);
caxis([0.5 1]);

%% Figure 12
load('sign_and_insign.mat')
figure('pos',[10 10 1000 1000]);
h=tight_subplot(3,2,[0.025,0.025],[0.05 0.035],[0.05,0.05]);
axes(h(1));
plot(1:366,precip_n_sign,'linewidth',2);
xlim([1 366]);
set(gca,'xtick',[],'fontsize',16);
set(gca,'fontsize',16);
ylim([0 6]);
text(18,5,'a) Significantly Negative','fontsize',16);

axes(h(2));
plot(1:366,precip_most_n,'linewidth',2);
xlim([1 366]);
set(gca,'xtick',[]);
set(gca,'fontsize',16);
ylim([0 15]);
text(18,13,'b) 92.25^{o}W 17.25^{o}N','fontsize',16);
%ylabel('mm','fontsize',16);

axes(h(3));
plot(1:366,precip_p_sign,'linewidth',2);
xlim([1 366]);
set(gca,'xtick',[],'fontsize',16);
set(gca,'fontsize',16);
ylim([0 6]);
text(18,5,'c) Significantly Positive','fontsize',16);

axes(h(4));
plot(1:366,precip_most_p,'linewidth',2);
xlim([1 366]);
set(gca,'xtick',[datenum(2016,1,15):30:datenum(2016,12,31)]-datenum(2016,1,1)+1,'xticklabels',{'J','F',...
    'M','A','M','J','J','A','S','O','N','D'},'fontsize',16);
set(gca,'fontsize',16);
ylim([0 15]);
text(18,13,'d) 108.75^{o}W 28.25^{o}N','fontsize',16);
%ylabel('mm','fontsize',16);


axes(h(5));
plot(1:366,precip_insign,'linewidth',2);
xlim([1 366]);
set(gca,'xtick',[datenum(2016,1,15):30:datenum(2016,12,31)]-datenum(2016,1,1)+1,'xticklabels',{'J','F',...
    'M','A','M','J','J','A','S','O','N','D'},'fontsize',16);
set(gca,'fontsize',16);
ylim([0 6]);
ylabel('mm','fontsize',16);
text(18,5,'e) Insignificant','fontsize',16);

axes(h(6));
box off
axis off

%% ts try
h1 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 
line(1979:2017,full_start,'color','r','parent',h1);
set(h1, 'ycolor', 'k','ylim',[155 180],'fontsize',16)
set(h1,'ytick',[157 177],'yticklabels',{'Jun 5th','Jun 25th'},'fontsize',16);
box off
ylabel('Onset','fontsize',16);
set(h1, 'yaxislocation', 'left', 'xtick', [])

h2 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 
line(1979:2017,full_end,'color','g','parent',h2);
set(h2, 'ycolor', 'r','ylim',[245 270],'fontsize',16)
set(h2,'ytick',[245 269],'yticklabels',{'Sep 1st','Sep 25th'},'fontsize',16);
box off
ylabel('End','fontsize',16);
set(h2, 'yaxislocation', 'right', 'xtick', [])

h3 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 
line(1979:2017,full_peak,'color','k','parent',h3);
set(h3, 'ycolor', 'b','ylim',[195 225],'fontsize',16);
set(h3,'ytick',[197 214],'yticklabels',{'Jul 15th','Aug 1st'},'fontsize',16);
box off
ylabel('Peak','fontsize',16);
set(h3, 'yaxislocation', 'right', 'xtick', [])

%% see
h1 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 
h2 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 
h3 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 


plot(h1,1:10,1:10,'color','k');
plot(h2,1:10,-(1:10),'color','r');
plot(h3,1:10,[1:5 3:7],'color','b');

set([h1,h2,h3],'Color', 'None')
set([h1 h2 h3], 'YAxisLocation', 'right')

h3.YAxis.TickGapOffset = 50;

h3.YRuler.Axle.VertexData(1,3:4)=12;


%% ts metrics
load('ts_metric');
figure('pos',[10 10 1500 1500]);

hf = axes('position', [0.1 0.1 0.5 0.85]);  
h1=plot(1979:2017,zscore(full_start),'k','linewidth',2);
hold on
h2=plot(1979:2017,zscore(full_end),'r','linewidth',2);
hold on
h3=plot(1979:2017,zscore(full_peak),'b','linewidth',2);
hold on
h4=plot(1979:2017,zscore(full_dur),'g','linewidth',2);
hold on
h5=plot(1979:2017,zscore(full_imsd),'m','linewidth',2);
legend([h1 h2 h3 h4 h5],{'Onset','End','Peak','Dur','I_{msd}'},'location','best','fontsize',16);

set(hf,'fontsize',16,'ytick',[]);

h1 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 
line(1979:2017,NaN(size(full_start)),'color','w','parent',h1);
set(h1, 'ycolor', 'k','ylim',[155 180],'fontsize',16)
set(h1,'ytick',[157 177],'yticklabels',{'Jun 5th','Jun 25th'},'fontsize',16);
box off
ylabel('Onset','fontsize',16);
set(h1, 'yaxislocation', 'left', 'xtick', [])

h2 = axes('position', [0.1 0.1 0.5 0.85],'color','none'); 
line(1979:2017,NaN(size(full_end)),'color','w','parent',h2);
set(h2, 'ycolor', 'r','ylim',[245 270],'fontsize',16)
set(h2,'ytick',[245 269],'yticklabels',{'Sep 1st','Sep 25th'},'fontsize',16);
box off
ylabel('End','fontsize',16);
set(h2, 'yaxislocation', 'right', 'xtick', [])

h3 = axes('position', [0.7 0.1 0.001 0.85],'color','none'); 
line(1979:2017,NaN(size(full_peak)),'color','w','parent',h3);
set(h3, 'ycolor', 'b','ylim',[195 225],'fontsize',16);
set(h3,'ytick',[197 214],'yticklabels',{'Jul 15th','Aug 1st'},'fontsize',16);
box off
ylabel('Peak','fontsize',16);
set(h3, 'yaxislocation', 'right', 'xtick', [])

h4 = axes('position', [0.8 0.1 0.001 0.85],'color','none'); 
line(1979:2017,NaN(size(full_dur)),'color','w','parent',h4);
set(h4, 'ycolor', 'g','ylim',[75 115],'fontsize',16)
box off
ylabel('Duration (days)','fontsize',16);
set(h4, 'yaxislocation', 'right', 'xtick', [])


h5 = axes('position', [0.9 0.1 0.001 0.85],'color','none'); 
line(1979:2017,NaN(size(full_imsd)),'color','w','parent',h5);
set(h5, 'ycolor','m', 'ylim',[0.45 0.65],'fontsize',16)
box off
ylabel('I_{msd}','fontsize',16);
set(h5, 'yaxislocation', 'right', 'xtick', [])

%% plot sst j ja s
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
load('sst_j_ja_s_hm');
figure('pos',[10 10 1500 1500]);

x_1=linspace(258,270,100);
y_1=12*ones(1,100);

x_2=linspace(258,270,100);
y_2=18*ones(1,100);

x_3=258*ones(1,100);
y_3=linspace(12,18,100);

x_4=270*ones(1,100);
y_4=linspace(12,18,100);


% June
subplot(3,2,1);
m_contourf(lon_sst,lat_sst,exp(sst_j)',linspace(nanmin(exp(sst_j(:))),nanmax(exp(sst_j(:))),100),'linestyle','none');
m_coast();
m_grid('xtick',[],'ytick',[]);
colormap(m_colmap('jet',11));
caxis([exp(15) exp(30.5)]);
hold on
m_line(x_1,y_1,'linewidth',2,'color','r');
hold on
m_line(x_2,y_2,'linewidth',2,'color','r');
hold on
m_line(x_3,y_4,'linewidth',2,'color','r');
hold on
m_line(x_4,y_4,'linewidth',2,'color','r');
m_text(243,4,'a) June','fontsize',16,'fontweight','bold')

% July-Aug

subplot(3,2,3);
m_contourf(lon_sst,lat_sst,exp(sst_ja)',linspace(nanmin(exp(sst_ja(:))),nanmax(exp(sst_ja(:))),100),'linestyle','none');
m_coast();
m_grid('xtick',[],'ytick',[]);
colormap(m_colmap('jet',11));
caxis([exp(15) exp(30.5)]);
hold on
m_line(x_1,y_1,'linewidth',2,'color','r');
hold on
m_line(x_2,y_2,'linewidth',2,'color','r');
hold on
m_line(x_3,y_4,'linewidth',2,'color','r');
hold on
m_line(x_4,y_4,'linewidth',2,'color','r');
m_text(243,4,'b) July and August','fontsize',16,'fontweight','bold')

% Sep

subplot(3,2,5);
m_contourf(lon_sst,lat_sst,exp(sst_s)',linspace(nanmin(exp(sst_s(:))),nanmax(exp(sst_s(:))),100),'linestyle','none');
m_coast();
m_grid('fontsize',16,'linestyle','none');
colormap(m_colmap('jet',11));
caxis([exp(15) exp(30.5)]);
hold on
m_line(x_1,y_1,'linewidth',2,'color','r');
hold on
m_line(x_2,y_2,'linewidth',2,'color','r');
hold on
m_line(x_3,y_3,'linewidth',2,'color','r');
hold on
m_line(x_4,y_4,'linewidth',2,'color','r');
m_text(243,4,'c) Sep','fontsize',16,'fontweight','bold')

hp4=get(subplot(3,2,5),'Position');   
s=colorbar('location','southoutside','ticks',linspace(exp(15),exp(30.5),12),'ticklabels',round(log(linspace(exp(15),exp(30.5),12)),1),'fontsize',12,'Position', [hp4(1)-0.03  hp4(2)-0.08  0.4  0.025]);

% hovm



subplot(3,2,[2 4 6]);
contourf(lon_sst(lon_sst<=270 & lon_sst>=258),datenum(2017,5,15):datenum(2017,10,15),sst_hm',linspace(nanmin(sst_hm(:)),nanmax(sst_hm(:)),100),'linestyle','none');
colormap(m_colmap('jet',11));
set(gca,'xtick',[260 264 268],'xticklabels',{'100^{o}W','94^{o}W','90^{o}W'},'fontsize',16);
set(gca,'ytick',[datenum(2017,6,1) datenum(2017,7,1) datenum(2017,8,1) datenum(2017,9,1) datenum(2017,10,1)],'yticklabels',{'Jun 1st','Jul 1st','Aug 1st','Sep 1st','Oct 1st'},'fontsize',16);
hold on
text(259,datenum(2017,5,25),'d) Hovmller (SST)','fontsize',16,'fontweight','bold');
hp4=get(subplot(3,2,[2 4 6]),'Position');   
s=colorbar('location','southoutside','fontsize',12,'Position', [hp4(1)-0.03  hp4(2)-0.08  0.4  0.025]);

%% Figure 5
size_d=0.8;color_d=[0 0 0];
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
load('lon_lat_ncep');
figure('pos',[10 10 1500 1500]);
h=tight_subplot(3,3,[0.025,0.025],[0.05 0.035],[0.05,0.05]);
load('precip_enso_new');
[lat,lon]=meshgrid(lat_ncep,lon_ncep);
color_here= customcolormap(linspace(0,1,11), {'#523107','#523107','#bf812c','#e2c17e','#f3e9c4','#f6f4f4','#cae9e3','#81cdc1','#379692','#01665e','#003d2e'},40);
color_here=color_here(end:-1:1,:);

axes(h(1));
m_contourf(lon_ncep,lat_ncep,precip_jun_en',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_jun_en==1),lat(see_jun_en==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'a) Jun','fontsize',16,'fontweight','bold');
title('El Nino','fontsize',16,'fontweight','bold');

axes(h(2));
m_contourf(lon_ncep,lat_ncep,precip_jun_la',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_jun_la==1),lat(see_jun_la==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'d) Jun','fontsize',16,'fontweight','bold');
title('La Nina','fontsize',16,'fontweight','bold');

axes(h(3));
m_contourf(lon_ncep,lat_ncep,precip_jun_ne',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_jun_ne==1),lat(see_jun_ne==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'g) Jun','fontsize',16,'fontweight','bold');
title('Neutral','fontsize',16,'fontweight','bold');

axes(h(4));
m_contourf(lon_ncep,lat_ncep,precip_ja_en',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_ja_en==1),lat(see_ja_en==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'b) Jul-Aug','fontsize',16,'fontweight','bold');

axes(h(5));
m_contourf(lon_ncep,lat_ncep,precip_ja_la',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_ja_la==1),lat(see_ja_la==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'e) Jul-Aug','fontsize',16,'fontweight','bold');

axes(h(6));
m_contourf(lon_ncep,lat_ncep,precip_ja_ne',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_ja_ne==1),lat(see_ja_ne==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'h) Jul-Aug','fontsize',16,'fontweight','bold');

axes(h(7));
m_contourf(lon_ncep,lat_ncep,precip_sep_en',-4:0.01:4,'linestyle','none');
m_grid('fontsize',16,'linestyle','none');
m_coast();
hold 
m_scatter(lon(see_sep_en==1),lat(see_sep_en==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'c) Sep','fontsize',16,'fontweight','bold');

axes(h(8));
m_contourf(lon_ncep,lat_ncep,precip_sep_la',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_sep_la==1),lat(see_sep_la==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'f) Sep','fontsize',16,'fontweight','bold');

axes(h(9));
m_contourf(lon_ncep,lat_ncep,precip_sep_ne',-4:0.01:4,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
hold 
m_scatter(lon(see_sep_ne==1),lat(see_sep_ne==1),size_d,color_d);
colormap(color_here)
caxis([-3 3]);
m_text(243,4,'i) Sep','fontsize',16,'fontweight','bold');

hp4=get(h(9),'Position');   
s=colorbar('fontsize',12,'Position', [hp4(1)+hp4(3)+0.01  hp4(2)-0.015  0.01  0.95],...
    'ticks',-3:0.3:3);
s.Label.String='mm/day';

%% Figure 6
size_d=0.5;color_d=[0 0 0];
load('lon_lat_ncep');
load('enso_patterns_new');
load('see');
color_here=m_colmap('diverging',24);
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
lon_era_low = linspace(nanmin(lon_ncep),nanmax(lon_ncep),121);
lat_era_low = linspace(nanmax(lat_ncep),nanmin(lat_ncep),61);
[lat,lon]=meshgrid(lat_era_low,lon_era_low);
figure('pos',[10 10 1500 1500]);
m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
h=tight_subplot(3,3,[0.025,0.025],[0.035 0.035],[0.03,0.08]);
axes(h(1));
m_contourf(lon_era_low,lat_era_low,pres_start_en',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_start_en==1),lat(see_start_en==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_start_en(1:4:end,1:4:end))',(v_start_en(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'fontsize',12,'linestyle','none');
title('El Nino','fontsize',16,'fontweight','bold');
caxis([-1.2 1.2]);
m_text(243,4,'a) Onset','fontsize',16,'fontweight','bold')
hold on
m_quiver(248,32,1,0,2,'color','k','autoscale','off')

axes(h(2));
m_contourf(lon_era_low,lat_era_low,pres_start_la',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_start_la==1),lat(see_start_la==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_start_la(1:4:end,1:4:end))',(v_start_la(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
title('La Nina','fontsize',16,'fontweight','bold');
caxis([-1.2 1.2]);
m_text(243,4,'b)','fontsize',16,'fontweight','bold')

axes(h(3));
m_contourf(lon_era_low,lat_era_low,pres_start_ne',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_start_ne==1),lat(see_start_ne==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_start_ne(1:4:end,1:4:end))',(v_start_ne(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
title('Neutral','fontsize',16,'fontweight','bold');
caxis([-1.2 1.2]);
m_text(243,4,'c)','fontsize',16,'fontweight','bold')

axes(h(4));
m_contourf(lon_era_low,lat_era_low,pres_peak_en',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_peak_en==1),lat(see_peak_en==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_peak_en(1:4:end,1:4:end))',(v_peak_en(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'fontsize',12,'linestyle','none');
caxis([-1.2 1.2]);
m_text(243,4,'d) Peak','fontsize',16,'fontweight','bold')

axes(h(5));
m_contourf(lon_era_low,lat_era_low,pres_peak_la',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_peak_la==1),lat(see_peak_la==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_peak_la(1:4:end,1:4:end))',(v_peak_la(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-1.2 1.2]);
m_text(243,4,'e)','fontsize',16,'fontweight','bold')

axes(h(6));
m_contourf(lon_era_low,lat_era_low,pres_peak_ne',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_peak_ne==1),lat(see_peak_ne==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_peak_ne(1:4:end,1:4:end))',(v_peak_ne(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-1.2 1.2]);
m_text(243,4,'f)','fontsize',16,'fontweight','bold')

axes(h(7));
m_contourf(lon_era_low,lat_era_low,pres_end_en',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_end_en==1),lat(see_end_en==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_end_en(1:4:end,1:4:end))',(v_end_en(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('fontsize',12,'linestyle','none');
caxis([-1.2 1.2]);
m_text(243,4,'g) End','fontsize',16,'fontweight','bold')

axes(h(8));
m_contourf(lon_era_low,lat_era_low,pres_end_la',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_end_la==1),lat(see_end_la==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_end_la(1:4:end,1:4:end))',(v_end_la(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('ytick',[],'fontsize',12,'linestyle','none');
caxis([-1.2 1.2]);
m_text(243,4,'k)','fontsize',16,'fontweight','bold')

axes(h(9));
m_contourf(lon_era_low,lat_era_low,pres_end_ne',-2:0.01:2,'linestyle','none');
m_coast();
hold on
m_scatter(lon(see_end_ne==1),lat(see_end_ne==1),size_d,color_d);
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_end_ne(1:4:end,1:4:end))',(v_end_ne(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('ytick',[],'fontsize',12,'linestyle','none');
caxis([-1.2 1.2]);
m_text(243,4,'l)','fontsize',16,'fontweight','bold')
hp4=get(h(9),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.025  0.9],'fontsize',14,'ticks',-1.2:0.2:1.2);
s.Label.String='hpa';

%% MSD phases patterns
load('/Volumes/mydrive/msd_data/msd_clim_era_low.mat','lon_era_low','lat_era_low');
load('msd_phase_patterns');
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
color_here=m_colmap('diverging',256);
[lat,lon]=meshgrid(lat_era_low,lon_era_low);
size_d=0.5;color_d=[0 0 0];

figure('pos',[10 10 1500 1500]);
m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
h=tight_subplot(2,2,[0.025,0.025],[0.035 0.035],[0.03,0.08]);
axes(h(1));
m_contourf(lon_era_low,lat_era_low,pres_msd_1',-2:0.01:2,'linestyle','none');
m_coast();
hold on
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_msd_1(1:4:end,1:4:end))',(v_msd_1(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-0.3 0.3]);
m_text(243,4,'a) P1','fontsize',16,'fontweight','bold')

axes(h(2));
m_contourf(lon_era_low,lat_era_low,pres_msd_2',-2:0.01:2,'linestyle','none');
m_coast();
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_msd_2(1:4:end,1:4:end))',(v_msd_2(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-0.3 0.3]);
m_text(243,4,'b) P2','fontsize',16,'fontweight','bold')

axes(h(3));
m_contourf(lon_era_low,lat_era_low,pres_msd_3',-2:0.01:2,'linestyle','none');
m_coast();
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_msd_3(1:4:end,1:4:end))',(v_msd_3(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('fontsize',16,'linestyle','none');
caxis([-0.3 0.3]);
m_text(243,4,'c) P3','fontsize',16,'fontweight','bold')

axes(h(4));
m_contourf(lon_era_low,lat_era_low,pres_msd_4',-2:0.01:2,'linestyle','none');
m_coast();
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_msd_4(1:4:end,1:4:end))',(v_msd_4(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-0.3 0.3]);
m_text(243,4,'d) P4','fontsize',16,'fontweight','bold')

hp4=get(h(4),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.025  0.9],'fontsize',14);
s.Label.String='hpa';

%% Figure 9
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
lon_era_low = linspace(nanmin(lon_ncep),nanmax(lon_ncep),121);
lat_era_low = linspace(nanmax(lat_ncep),nanmin(lat_ncep),61);
load('mjo_patterns_msd');
figure('pos',[10 10 1500 1500]);
m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
h=tight_subplot(2,2,[0.025,0.025],[0.035 0.035],[0.03,0.08]);
axes(h(1));
m_contourf(lon_era_low,lat_era_low,pres_mjo_81',-2:0.01:2,'linestyle','none');
m_coast();
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_mjo_81(1:4:end,1:4:end))',(v_mjo_81(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-0.6 0.6]);
m_text(243,4,'a) Phase 8-1','fontsize',16,'fontweight','bold')

axes(h(2));
m_contourf(lon_era_low,lat_era_low,pres_mjo_23',-2:0.01:2,'linestyle','none');
m_coast();
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_mjo_23(1:4:end,1:4:end))',(v_mjo_23(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-0.6 0.6]);
m_text(243,4,'b) Phase 2-3','fontsize',16,'fontweight','bold')

axes(h(3));
m_contourf(lon_era_low,lat_era_low,pres_mjo_45',-2:0.01:2,'linestyle','none');
m_coast();
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_mjo_45(1:4:end,1:4:end))',(v_mjo_45(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('fontsize',16,'linestyle','none');
caxis([-0.6 0.6]);
m_text(243,4,'c) Phase 4-5','fontsize',16,'fontweight','bold')

axes(h(4));
m_contourf(lon_era_low,lat_era_low,pres_mjo_67',-2:0.01:2,'linestyle','none');
m_coast();
colormap(color_here);
hold on
[lon,lat]=meshgrid(lon_era_low,lat_era_low);
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),(u_mjo_67(1:4:end,1:4:end))',(v_mjo_67(1:4:end,1:4:end))',2,'color','k','autoscale','off');
hold on
m_grid('xtick',[],'ytick',[]);
caxis([-0.6 0.6]);
m_text(243,4,'d) Phase 6-7','fontsize',16,'fontweight','bold')

hp4=get(h(4),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.025  0.9],'fontsize',14,'ticks',-0.6:0.2:0.6);
s.Label.String='hpa';

%% proportion
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
addpath('/Volumes/mydirve/cloud_annual/m_map');
load('prop_msd_mjo.mat');
load('msd_collection_ncep');
figure('pos',[10 10 1500 1500]);
m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
h=tight_subplot(4,4,[0.025,0.025],[0.035 0.035],[0.03,0.08]);

axes(h(1));
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,1))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'fontsize',12,'fontweight','bold','linestyle','none');
caxis([0 0.6]);
m_text(243,4,'a) P1','fontsize',16,'fontweight','bold')
title('Phase 8-1','fontsize',16);

axes(h(2));
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,2))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'b)','fontsize',16,'fontweight','bold')
title('Phase 2-3','fontsize',16);

axes(h(3));
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,3))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'c)','fontsize',16,'fontweight','bold')
title('Phase 4-5','fontsize',16);

axes(h(4));
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,4))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'d)','fontsize',16,'fontweight','bold')
title('Phase 6-7','fontsize',16);

axes(h(5));
m_contourf(lon_ncep,lat_ncep,(prop_2(:,:,1))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'fontsize',12,'fontweight','bold','linestyle','none');
caxis([0 0.6]);
m_text(243,4,'e) P2','fontsize',16,'fontweight','bold')

axes(h(6));
m_contourf(lon_ncep,lat_ncep,(prop_2(:,:,2))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'f)','fontsize',16,'fontweight','bold')

axes(h(7));
m_contourf(lon_ncep,lat_ncep,(prop_2(:,:,3))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'g)','fontsize',16,'fontweight','bold')

axes(h(8));
m_contourf(lon_ncep,lat_ncep,(prop_2(:,:,4))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'h)','fontsize',16,'fontweight','bold')

axes(h(9));
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,1))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'fontsize',12,'fontweight','bold','linestyle','none');
caxis([0 0.6]);
m_text(243,4,'i) P3','fontsize',16,'fontweight','bold')

axes(h(10));
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,2))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'j)','fontsize',16,'fontweight','bold')

axes(h(11));
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,3))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'k)','fontsize',16,'fontweight','bold')

axes(h(12));
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,4))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('xtick',[],'ytick',[]);
caxis([0 0.6]);
m_text(243,4,'l)','fontsize',16,'fontweight','bold')

axes(h(13));
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,1))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('fontsize',12,'linestyle','none');
caxis([0 0.6]);
m_text(243,4,'m) P4','fontsize',16,'fontweight','bold')

axes(h(14));
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,2))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('ytick',[],'fontsize',12,'fontweight','bold','linestyle','none');
caxis([0 0.6]);
m_text(243,4,'n)','fontsize',16,'fontweight','bold')

axes(h(15));
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,3))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('ytick',[],'fontsize',12,'fontweight','bold','linestyle','none');
caxis([0 0.6]);
m_text(243,4,'o)','fontsize',16,'fontweight','bold')

axes(h(16));
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,4))',0:0.01:1,'linestyle','none');
m_coast();
colormap(m_colmap('jet',11));
hold on
m_grid('ytick',[],'fontsize',12,'fontweight','bold','linestyle','none');
caxis([0 0.6]);
m_text(243,4,'p)','fontsize',16,'fontweight','bold')

hp4=get(h(16),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.025  0.9],'fontsize',14);
s.Label.String='';
%%
load('s_p');
figure('pos',[10 10 1500 1500]);
m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
h=tight_subplot(2,2,[0.025,0.025],[0.035 0.035],[0.03,0.05]);

axes(h(1));
m_contourf(lon_ncep,lat_ncep,(s_p2(:,:,1))',0:0.01:35,'linestyle','none');
m_coast();
m_grid('xtick',[],'ytick',[],'linestyle','none');
colormap(jet);
caxis([0 20]);
m_text(243,4,'a) End El','fontsize',16,'fontweight','bold')

axes(h(2));
m_contourf(lon_ncep,lat_ncep,(s_p2(:,:,2))',0:0.01:35,'linestyle','none');
m_coast();
m_grid('xtick',[],'ytick',[],'linestyle','none');
colormap(jet);
caxis([0 20]);
m_text(243,4,'b) End no - El','fontsize',16,'fontweight','bold')

hp4=get(h(2),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.01  0.46],'fontsize',14);
s.Label.String='mm/day';

axes(h(3));
m_contourf(lon_ncep,lat_ncep,(s_pmin(:,:,1))',0:0.01:35,'linestyle','none');
m_coast();
m_grid('fontsize',16,'fontweight','bold','linestyle','none');
colormap(jet);
caxis([0 10]);
m_text(243,4,'c) Peak La','fontsize',16,'fontweight','bold')

axes(h(4));
m_contourf(lon_ncep,lat_ncep,(s_pmin(:,:,2))',0:0.01:35,'linestyle','none');
m_coast();
m_grid('xtick',[],'ytick',[],'linestyle','none');
colormap(jet);
caxis([0 10]);
m_text(243,4,'c) Peak no- La','fontsize',16,'fontweight','bold')


hp4=get(h(4),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.01  0.46],'fontsize',14);
s.Label.String='mm/day';


%% Figure 13
addpath('/Volumes/mydirve/cloud_annual/m_map');
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
load('full_model');
load('msd_collection_ncep');
logic_msd=NaN(120,60);
color_here=m_colmap('diverging',11);

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

for i=1:size(loc_full_data,1);
    logic_msd(loc_full_data(i,1),loc_full_data(i,2))=0;
end

for i=1:size(loc_unique,1);
    logic_msd(loc_unique(i,1),loc_unique(i,2))=1;
end

r_mean=nanmean(R_adj_mjo);
coef_4_mean=nanmean(coef_4);

r_plot_full=NaN(120,60);
coef_4_mean_plot=NaN(120,60);

for i=1:size(loc_full_data,1);
    r_plot_full(loc_full_data(i,1),loc_full_data(i,2))=r_mean(i);
    coef_4_mean_plot(loc_full_data(i,1),loc_full_data(i,2))=coef_4_mean(i);
    %p_logic_full(loc_full_data(i,1),loc_full_data(i,2))=nansum(p_4(i)<=0.05);
end


figure('pos',[10 10 1000 1000]);
h=subplot(2,1,2);
m_contourf(lon_ncep,lat_ncep,(coef_4_mean_plot')*1e07,-4.7:0.01:4.7,'linestyle','none');
m_grid('fontsize',16);
m_coast();
colormap(h,color_here);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
title('b4 \times 10^{7}','fontsize',16);
s=colorbar('fontsize',16);
s.Label.String='mm/day^{4}';
caxis([-4 4]);

h=subplot(2,1,1);
m_contourf(lon_ncep,lat_ncep,(r_plot_full'),0:0.01:1,'linestyle','none');
m_grid('fontsize',16);
m_coast();
colormap(h,m_colmap('jet',11));
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
title('adjusted R^{2}','fontsize',16);
colorbar('fontsize',16);
caxis([0.4 0.8]);

%% mjo R_adj
addpath('/Volumes/mydirve/cloud_annual/m_map');
addpath('/Volumes/mydirve/cloud_annual/tight_subplot');
load('annual_model_used');
load('annual_poly_model','loc_full_data');
load('msd_collection_ncep');
logic_msd=NaN(120,60);
color_here=m_colmap('diverging',11);

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

for i=1:size(loc_full_data,1);
    logic_msd(loc_full_data(i,1),loc_full_data(i,2))=0;
end

for i=1:size(loc_unique,1);
    logic_msd(loc_unique(i,1),loc_unique(i,2))=1;
end

r_mean=nanmean(R_adj_mjo);

coef_rmm1_mean=nanmean(coef_rmm1);

coef_rmm2_mean=nanmean(coef_rmm2);

coef_4_mean=nanmean(coef_4);

r_plot_full=NaN(120,60);

coef_rmm1_mean_plot=NaN(120,60);

coef_rmm2_mean_plot=NaN(120,60);

coef_4_mean_plot=NaN(120,60);

for i=1:size(loc_full_data,1);
    r_plot_full(loc_full_data(i,1),loc_full_data(i,2))=r_mean(i);
    coef_rmm1_mean_plot(loc_full_data(i,1),loc_full_data(i,2))=coef_rmm1_mean(i);
    coef_rmm2_mean_plot(loc_full_data(i,1),loc_full_data(i,2))=coef_rmm2_mean(i);
    %p_logic_full(loc_full_data(i,1),loc_full_data(i,2))=nansum(p_4(i)<=0.05);
    coef_4_mean_plot(loc_full_data(i,1),loc_full_data(i,2))=coef_4_mean(i);
end

figure('pos',[10 10 1500 1500]);
h=tight_subplot(2,2,[0.025,0.025],[0.035 0.035],[0.03,0.03]);
axes(h(1));
m_contourf(lon_ncep,lat_ncep,(r_plot_full'),0:0.01:1,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(1),m_colmap('jet',11));
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
m_text(243,4,'a) adjusted R^{2}','fontsize',16,'fontweight','bold')
s=colorbar('fontsize',16);
caxis([0.4 0.8]);

axes(h(2));
m_contourf(lon_ncep,lat_ncep,(coef_4_mean_plot')*1e07,-4.7:0.01:4.7,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(2),color_here);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
m_text(243,4,'b) b4 \times 10^{7}','fontsize',16,'fontweight','bold')
s=colorbar('fontsize',16);
s.Label.String='10^{7} mm/day^{4}';
caxis([-4 4]);

axes(h(3));
m_contourf(lon_ncep,lat_ncep,(coef_rmm1_mean_plot'),-2:0.01:2,'linestyle','none');
m_grid('fontsize',16,'linestyle','none');
m_coast();
colormap(h(3),color_here);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
m_text(243,4,'c) a1','fontsize',16,'fontweight','bold')
s=colorbar('fontsize',16);
s.Label.String='mm/rmm';
caxis([-1.5 1.5]);

axes(h(4));
m_contourf(lon_ncep,lat_ncep,(coef_rmm2_mean_plot'),-2:0.01:2,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(4),color_here);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
hold on
[lon,lat]=meshgrid(lon_ncep,lat_ncep);
m_scatter(lon((logic_msd==1)'),lat((logic_msd==1)'),1.75,'k');
m_text(243,4,'d) a2','fontsize',16,'fontweight','bold')
s=colorbar('fontsize',16);
s.Label.String='mm/rmm';
caxis([-1.5 1.5]);

%% plot for p
figure('pos',[10 10 1500 1500]);
h=tight_subplot(3,3,[0.025,0.025],[0.035 0.035],[0.03,0.05]);
axes(h(1));
m_contourf(lon_ncep,lat_ncep,(s_p1(:,:,1))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(1),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'a) P1','fontsize',16,'fontweight','bold')
caxis([0 20]);
title('El Nino','fontsize',16);

axes(h(2));
m_contourf(lon_ncep,lat_ncep,(s_p1(:,:,2))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(2),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'b)','fontsize',16,'fontweight','bold')
caxis([0 20]);
title('La Nina','fontsize',16);

axes(h(3));
m_contourf(lon_ncep,lat_ncep,(s_p1(:,:,3))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(3),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'c)','fontsize',16,'fontweight','bold')
caxis([0 20]);
title('Neutral','fontsize',16);

axes(h(4));
m_contourf(lon_ncep,lat_ncep,(s_p2(:,:,1))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(4),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'d) P2','fontsize',16,'fontweight','bold')
caxis([0 20]);

axes(h(5));
m_contourf(lon_ncep,lat_ncep,(s_p2(:,:,2))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(5),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'e)','fontsize',16,'fontweight','bold')
caxis([0 20]);

axes(h(6));
m_contourf(lon_ncep,lat_ncep,(s_p2(:,:,3))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(6),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'f)','fontsize',16,'fontweight','bold')
caxis([0 20]);

hp4=get(h(6),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.01  0.57],'fontsize',14);
s.Label.String='mm/day';

axes(h(7));
m_contourf(lon_ncep,lat_ncep,(s_pmin(:,:,1))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(7),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'g) P_{min}','fontsize',16,'fontweight','bold')
m_grid('fontsize',16,'linestyle','none');
caxis([0 10]);

axes(h(8));
m_contourf(lon_ncep,lat_ncep,(s_pmin(:,:,2))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(8),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'h)','fontsize',16,'fontweight','bold')
caxis([0 10]);

axes(h(9));
m_contourf(lon_ncep,lat_ncep,(s_pmin(:,:,3))',0:0.1:30,'linestyle','none');
m_grid('xtick',[],'ytick',[]);
m_coast();
colormap(h(9),jet);
% hold on
% m_contour(lon_ncep,lat_ncep,logic_msd','color','k');
m_text(243,4,'i)','fontsize',16,'fontweight','bold')
caxis([0 10]);

hp4=get(h(9),'Position');   
s=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.01  0.01  0.27],'fontsize',14);
s.Label.String='mm/day';



















    




    














    
    
    