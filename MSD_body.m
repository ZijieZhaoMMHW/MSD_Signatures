%% just try!!!
cd('/Volumes/mydrive/msd_data/new_rainfall');
p=double(ncread('rainfall_low_1993_2017.nc','time'));
precip=double(ncread('rainfall_low_1993_2017.nc','tp'));

time_used_large=datevec(datenum(1900,1,1,0,0,0)+double(p)./24-0.5);
time_used_large=time_used_large(:,1:3);

time_unique_large=unique(time_used_large,'rows');

time_used=time_unique_large;

precip_daily=NaN(121,61,9131);

precip=double(precip);

precip(precip==-32767)=nan;

for i=1:size(time_used,1);
    index_here=find(time_used_large(:,1)==time_used(i,1) & time_used_large(:,2)==time_used(i,2) & time_used_large(:,3)==time_used(i,3));
    precip_here=precip(:,:,index_here);
    precip_here=nansum(precip_here,3);
    precip_daily(:,:,i)=precip_here;
end

precip_clim=[];

time_unique=unique(time_used(:,2:3),'rows');

for i=1:size(time_unique,1);
    month_here=time_unique(i,1);
    
    day_here=time_unique(i,2);
    
    index_used=find(time_used(:,2)==month_here & time_used(:,3)==day_here);
    
    precip_here=precip_daily(:,:,index_used);
    
    precip_here=nanmean(precip_here,3);
    
    precip_clim=cat(3,precip_clim,precip_here);
end

precip_clim_smooth=smoothdata(precip_clim,3,'gaussian',31)*1000;

imsd_climatology=NaN(121,61);

for i=1:size(precip_clim,1);
    for j=1:size(precip_clim,2);
        precip_here=squeeze(precip_clim_smooth(i,j,:));
        if nansum(isnan(precip_here))~=length(precip_here);
            potential_p1=[];
            potential_p2=[];
            potential_loc1=[];
            potential_loc2=[];
            
            [p,loc]=findpeaks(precip_here);
            [p,order]=sort(p,'descend');
            loc=loc(order);
            
            p_here=p;
            loc_here=loc;
            
            for m=1:(length(p_here)-1);
                for n=m+1:length(p_here);
                    p_1_2=p_here([m n]);
                    loc_1_2=loc_here([m n]);
                    
                    [loc_f,one_or_two]=nanmin(loc_1_2);
                    [loc_b,two_or_one]=nanmax(loc_1_2);
                    
                    p_f=p_1_2(one_or_two);
                    p_b=p_1_2(two_or_one);
                    
                    precip_f=precip_here(1:loc_f);
                    precip_b=precip_here(loc_b:length(precip_here));
                    
                    mdl_f=fitlm(1:loc_f,precip_f);
                    mdl_b=fitlm(loc_b:length(precip_here),precip_b);
                    
                    coef_f=mdl_f.Coefficients.Estimate(2);
                    pv_f=mdl_f.Coefficients.pValue(2);
                    coef_b=mdl_b.Coefficients.Estimate(2);
                    pv_b=mdl_b.Coefficients.pValue(2);
                    
                    precip_f_b=precip_here(loc_f:loc_b);
                    period_f_b=loc_f:loc_b;
                    [p_s,l_s]=nanmin(precip_f_b);
                    
                    loc_s=period_f_b(l_s);
                    
                    mdl_h_f=fitlm(loc_f:loc_s,precip_here(loc_f:loc_s));
                    mdl_h_b=fitlm(loc_s:loc_b,precip_here(loc_s:loc_b));
                    
                    coef_h_f=mdl_h_f.Coefficients.Estimate(2);
                    coef_h_b=mdl_h_b.Coefficients.Estimate(2);
                    pv_h_f=mdl_h_f.Coefficients.pValue(2);
                    pv_h_b=mdl_h_b.Coefficients.pValue(2);
                    
                    if (loc_b-loc_f+1)>=30 && (loc_b-loc_f+1)<=124 && coef_f>0 && pv_f<0.05 && coef_b<0 && pv_b<0.05 &&...
                            coef_h_f<0 && pv_h_f<0.05 && coef_h_b>0 && pv_h_b<0.05
                        potential_p1=[potential_p1;p_f];
                        potential_p2=[potential_p2;p_b];
                        potential_loc1=[potential_loc1;loc_f];
                        potential_loc2=[potential_loc2;loc_b];
                    end
                end
            end
            
            if ~isempty(potential_p1)
                sum_p1_p2=potential_p1+potential_p2;
                p1_final=potential_p1(sum_p1_p2==nanmax(sum_p1_p2));
                p2_final=potential_p2(sum_p1_p2==nanmin(sum_p1_p2));
                loc1_final=potential_loc1(sum_p1_p2==nanmin(sum_p1_p2));
                loc2_final=potential_loc2(sum_p1_p2==nanmin(sum_p1_p2));
                pmin=nanmin(precip_here(loc1_final:loc2_final));
                
                pmax=nanmax([p1_final p2_final]);
                pmin=nanmin(precip_here(loc1_final:loc2_final));
                
                imsd_here=(pmax-pmin)./pmax;
                
                imsd_climatology(i,j)=imsd_here;
            end
            
        end
    end
end

%% detect 1979 - 2017
cd('/Volumes/mydrive/msd_data/rainfall');

lon=double(ncread('precip.1999.nc','lon'));
lat=double(ncread('precip.1999.nc','lat'));

lon_used=find(lon>=240 & lon<=300);
lat_used=find(lat<=30 & lat>=0);

precip=[];

for i=1979:2017
    file_here=['precip.' num2str(i) '.nc'];
    precip_here=ncread(file_here,'precip',[lon_used(1) lat_used(1) 1],[length(lon_used),length(lat_used),inf]);
    precip=cat(3,precip,precip_here);
end

precip=double(precip);

precip(abs(precip)>1000)=nan;

time_used=datevec(datenum(1979,1,1):datenum(2017,12,31));

time_used=time_used(:,1:3);

precip_clim=[];

time_unique=unique(time_used(:,2:3),'rows');

for i=1:size(time_unique,1);
    month_here=time_unique(i,1);
    day_here=time_unique(i,2);
    
    index_used=find(time_used(:,2)==month_here & time_used(:,3)==day_here);
    
    precip_here=precip(:,:,index_used);
    
    precip_here=nanmean(precip_here,3);
    
    precip_clim=cat(3,precip_clim,precip_here);
end

addpath('/Users/MAC/Desktop/project/����')

binary_msd=NaN(120,60);

imsd_climatology=NaN(120,60);

dif_clim=NaN(120,60);

ind1_cli=NaN(120,60);
ind2_cli=NaN(120,60);

for i=1:size(precip_clim,1)
    for j=1:size(precip_clim,2)
        precip_here=squeeze(precip_clim(i,j,:));
        if nansum(isnan(precip_here))~=length(precip_here)
            clim_here=smoothdata(precip_here,1,'gaussian',31);
            %clim_here=wrunavg(precip_here,31,3,1);
            %clim_here=precip_here;
            
            period_pmax1=datenum(0,5,15):datenum(0,7,15);
            
            period_pmax2=datenum(0,8,15):datenum(0,10,15);
            
            rain_p1=clim_here(period_pmax1);
            rain_p2=clim_here(period_pmax2);
            
            [pmax1,loc1]=nanmax(rain_p1);
            [pmax2,loc2]=nanmax(rain_p2);
            
            ind1=period_pmax1(loc1);
            ind2=period_pmax2(loc2);
            
            
            mdl_start=fitlm(1:ind1,clim_here(1:ind1));
            mdl_end=fitlm(ind2:length(clim_here),clim_here(ind2:end));
            trend_start=mdl_start.Coefficients.Estimate(2);
            p_start=mdl_start.Coefficients.pValue(2);
            trend_end=mdl_end.Coefficients.Estimate(2);
            p_end=mdl_end.Coefficients.pValue(2);

            
            if ind1==period_pmax1(end) || ind2==period_pmax2(1) ... %%|| nanmax(clim_here(ind1+15:ind2-15))>=pmax1 || nanmax(clim_here(ind1+15:ind2-15))>=pmax2 ...
                    || (~ismember(nanmax(clim_here),[pmax1,pmax2])) || trend_start<=0 || trend_end>=0 || p_start>0.05 || p_end>0.05

                ind1=NaN;
                ind2=NaN;
            else
                pmax=nanmax([pmax1 pmax2]);
                pmin=nanmean(clim_here(ind1:ind2));
                
                
                imsd=(pmax-pmin)./pmax;
                binary_msd(i,j)=1;
                imsd_climatology(i,j)=imsd;
                dif_clim(i,j)=pmax-pmin;
                ind1_cli(i,j)=ind1;
                ind2_cli(i,j)=ind2;
            end
            
        end
    end
end

imsd_climatology_ncep = imsd_climatology;

lon_ncep=lon(lon_used);
lat_ncep=lat(lat_used);

dif_clim_ncep=dif_clim;

[x,y]=find(binary_msd==1);

msd_collection=[];

time_here=datenum(time_used);

precip_smooth=smoothdata(precip,3,'gaussian',31);

for i=1:length(x)
    precip_here=squeeze(precip_smooth(x(i),y(i),:));
    if nansum(isnan(precip_here))~=length(precip_here)
        %precip_here=smoothdata(precip_here,1,'gaussian',31);
        %clim_here=wrunavg(precip_here,31,3,1);
        %clim_here=precip_here;
        
        year_used=unique(time_used(:,1));
        
        for j=1:length(year_used)
            year_here=year_used(j);
            period_pmax1=datenum(year_here,5,15):datenum(year_here,7,15);
            index_pmax1=find(ismember(time_here,period_pmax1));
            
            period_pmax2=datenum(year_here,8,15):datenum(year_here,10,15);
            index_pmax2=find(ismember(time_here,period_pmax2));
            
            rain_p1=precip_here(index_pmax1);
            rain_p2=precip_here(index_pmax2);
            
            [pmax1,loc1]=nanmax(rain_p1);
            [pmax2,loc2]=nanmax(rain_p2);
            
            ind1=index_pmax1(loc1);
            ind2=index_pmax2(loc2);
            
            precip_year=precip_here((datenum(year_here,1,1)-datenum(1979,1,1)+1):(datenum(year_here,12,31)-datenum(1979,1,1)+1));
            
            mdl_start=fitlm(1:(time_here(ind1)-datenum(year_here,1,1)+1),precip_year(1:time_here(ind1)-datenum(year_here,1,1)+1));
            mdl_end=fitlm((time_here(ind2)-datenum(year_here,1,1)+1):length(precip_year),precip_year(time_here(ind2)-datenum(year_here,1,1)+1:end));
            trend_start=mdl_start.Coefficients.Estimate(2);
            p_start=mdl_start.Coefficients.pValue(2);
            trend_end=mdl_end.Coefficients.Estimate(2);
            p_end=mdl_end.Coefficients.pValue(2);

            
            if ind1==index_pmax1(end) || ind2==index_pmax2(1) ||  ...
                    (~ismember(nanmax(precip_here((datenum(year_here,1,1)-datenum(1979,1,1)+1):(datenum(year_here,12,31)-datenum(1979,1,1)+1))),[pmax1,pmax2])) || trend_start<=0 || trend_end>=0 || p_start>0.05 || p_end>0.05

                ind1=NaN;
                ind2=NaN;
            else
                pmax=nanmax([pmax1 pmax2]);
                pmin=nanmean(precip_here(ind1:ind2));
                [t_pmin,loc_tmin]=nanmin(precip_here(ind1:ind2));
                ind_full=ind1:ind2;
                
                imsd=(pmax-pmin)./pmax;
                imsd_here=[year_here x(i) y(i) time_here(ind1) time_here(ind1)-datenum(year_here,1,1)+1 time_here(ind2) time_here(ind2)-datenum(year_here,1,1)+1 time_here(ind_full(loc_tmin)) time_here(ind_full(loc_tmin))-datenum(year_here,1,1)+1 pmax pmin imsd];
            
            msd_collection=[msd_collection;imsd_here];
                
            end
            
            
        end
    end
end

precip_ncep=precip;
msd_collection_ncep=msd_collection;
lon_ncep=lon(lon_used);
lat_ncep=lat(lat_used);
dif_clim_ncep=dif_clim;
                
save msd_collection_ncep msd_collection_ncep lat_ncep lon_ncep imsd_climatology_ncep

%% detect MSD states

load('msd_collection_ncep');

loc_unique=unique(msd_collection(:,2:3),'rows');

full_start=NaN(120,60);
full_end=NaN(120,60);
full_peak=NaN(120,60);
full_dur=NaN(120,60);

for i=1:size(loc_unique,1);
    loc_here=loc_unique(i,:);
    msd_here=msd_collection_ncep(msd_collection(:,2)==loc_here(1) & msd_collection(:,3)==loc_here(2),:);
    start_here=msd_here(:,5);
    end_here=msd_here(:,7);
    peak_here=msd_here(:,9);
    dur_here=end_here-start_here+1;
    
    full_start(loc_here(1),loc_here(2))=nanmean(start_here);
    full_end(loc_here(1),loc_here(2))=nanmean(end_here);
    full_peak(loc_here(1),loc_here(2))=nanmean(peak_here);
    full_dur(loc_here(1),loc_here(2))=nanmean(dur_here);
end

save mean_char_msd full_start full_end full_peak full_dur imsd_climatology_ncep
save daily_cpc_precip precip

%% detect MSD trends (Imsd)
imsd_spatial=NaN(120,60,39);

load('msd_collection_ncep');
loc= msd_collection_ncep(:,2:3);
loc_unique=unique(loc,'rows');

for i=1:size(loc_unique,1)
    loc_here=loc_un(i,:);
    msd_here=msd_collection_ncep(msd_collection_ncep(:,2)==loc_here(1) & msd_collection_ncep(:,3)==loc_here(2),:);
    
    time_here=1979:2017;
    
    imsd_here=zeros(39,1);
    
    for j=1:size(msd_here,1);
        year_here=msd_here(j,1);
        imsd_here(year_here-1979+1)=msd_here(j,end);
    end
    
    imsd_spatial(loc_here(1),loc_here(2),:)=imsd_here;
end

trend_spatial=NaN(120,60);
p_spatial=NaN(120,60);

for i=1:size(loc_unique,1);
    loc_here=loc_un(i,:);
    ts_here=squeeze(imsd_spatial(loc_here(1),loc_here(2),:));
    mdl_here=fitlm(1:39,ts_here);
    t_here=mdl_here.Coefficients.Estimate(2);
    p_here=mdl_here.Coefficients.pValue(2);
    trend_spatial(loc_here(1),loc_here(2))=t_here;
    p_spatial(loc_here(1),loc_here(2))=p_here;
end

% seems no trend......

% what about annual number and enso

year_unique=1979:2017;

msd_num=[];

for i=1979:2017
    msd_num=[msd_num;size(msd_collection_ncep(msd_collection_ncep(:,1)==i,:),1)];
end

full_index=[];

for i=1979:2017
    full_index=[full_index;datenum(i,7,1)];
end

plot(full_index,msd_num);
xlim([datenum(1979,1,1) datenum(2017,12,31)]);
ylim([350 750]);

hold on

% en=[1980 1983 1987 1988 1992 1995 1998 2003 2005 2007 2010 2015 2016];
% la=[1984 1985 1989 1996 1999 2000 2001 2006 2008 2009 2011 2012 2017];
en=[1980 1983 1987 1988 1992 1995 1998 2003 2005 2007 2010 2015 2016]-1;
la=[1984 1985 1989 1996 1999 2000 2001 2006 2008 2009 2011 2012 2017]-1;

msd_prop_e=NaN(120,60);
msd_prop_l=NaN(120,60);
msd_prop_n=NaN(120,60);

for i=1:size(loc_unique,1)
    loc_here=loc_unique(i,:);
    msd_here=msd_collection_ncep(msd_collection_ncep(:,2)==loc_here(1) & msd_collection_ncep(:,3)==loc_here(2),:);
    %     en_prop=nansum(ismember(msd_here(:,1),en))./size(msd_here,1);
    %     la_prop=nansum(ismember(msd_here(:,1),la))./size(msd_here,1);
    %     n_prop=1-en_prop-la_prop;
    %
    %     msd_prop_e(loc_here(1),loc_here(2))=en_prop;
    %     msd_prop_l(loc_here(1),loc_here(2))=la_prop;
    %     msd_prop_n(loc_here(1),loc_here(2))=n_prop;
    
    msd_prop_e(loc_here(1),loc_here(2))=nanmean(msd_here(ismember(msd_here(:,1),en),7)-msd_here(ismember(msd_here(:,1),en),5));
    msd_prop_l(loc_here(1),loc_here(2))=nanmean(msd_here(ismember(msd_here(:,1),la),7)-msd_here(ismember(msd_here(:,1),la),5));
    msd_prop_n(loc_here(1),loc_here(2))=nanmean(msd_here(~ismember(msd_here(:,1),en) & ~ismember(msd_here(:,1),la),7)-msd_here(~ismember(msd_here(:,1),en) & ~ismember(msd_here(:,1),la),5));
    
end

% subplot(3,1,1);
%  m_contourf(lon_ncep,lat_ncep,(msd_prop_e)',50:0.1:150,'linestyle','none');
%  colormap(jet);
%  m_grid;
%  m_coast();
%  caxis([50 130]);
%  subplot(3,1,2);
%  m_contourf(lon_ncep,lat_ncep,(msd_prop_l)',50:0.1:150,'linestyle','none');
%  colormap(jet);
%  m_grid;
%  m_coast();
%   caxis([50 130]);
%  subplot(3,1,3);
%  m_contourf(lon_ncep,lat_ncep,(msd_prop_n)',50:0.1:150,'linestyle','none');
%  colormap(jet);
%  m_grid;
%  m_coast();
%   caxis([50 130]);

en_event_ts=[];
la_event_ts=[];

for i=1:length(en);
    en_event_ts=[en_event_ts;[datenum(en(i),1,1) datenum(en(i),12,31)]];
end

for i=1:length(la);
    la_event_ts=[la_event_ts;[datenum(la(i),1,1) datenum(la(i),12,31)]];
end

for j=1:size(en_event_ts,1)
    start_ind=en_event_ts(j,1);
    end_ind=en_event_ts(j,2);
    
    x=[start_ind end_ind end_ind start_ind];
    y=[-1000 -1000 1000 1000];
    hold on
    fill(x,y,[1 0.5 0.5],'linestyle','none','FaceAlpha',0.5);
end

for j=1:size(la_event_ts,1)
    start_ind=la_event_ts(j,1);
    end_ind=la_event_ts(j,2);
    
    x=[start_ind end_ind end_ind start_ind];
    y=[-1000 -1000 1000 1000];
    hold on
    fill(x,y,[0.5 0.5 1],'linestyle','none','FaceAlpha',0.5);
end

%% MSD performance in EOF

% full rainfall
load('precip_clim');
precip_clim_smooth = smoothdata(precip_clim,3,'gaussian',31);

precip_2d=NaN(366,7200);
for i=1:size(precip_clim_smooth,3);
    data_here=squeeze(precip_clim_smooth(:,:,i));
    data_here=data_here(:);
    precip_2d(i,:)=data_here;
end

[no_need,mu,sd]=zscore(precip_2d);


land_log=nansum(isnan(no_need))==size(no_need,1);
land_loc=find(nansum(isnan(no_need))==size(no_need,1));

precip_2d_no_nan=precip_2d(:,~land_log);
sd=sd(~land_log);
mu=mu(~land_log);

F=detrend(precip_2d_no_nan,0);
[EOFs,lambda,contribution,PCs]=EOF_analysis(F');

% rainfall over MSD area

precip_2d=NaN(366,987);

for i=1:size(loc_unique,1)
    precip_here=precip_clim_smooth(loc_unique(i,1),loc_unique(i,2),:);
    precip_here=precip_here(:);
    precip_2d(:,i)=precip_here;
    
    
end

F=detrend(precip_2d,0);
[EOFs,lambda,contribution,PCs]=EOF_analysis(F');

EOF_full_1=NaN(120,60);
EOF_1=EOFs(:,1);

for i=1:size(loc_unique);
    loc_here=loc_unique(i,:);
    EOF_full_1(loc_here(1),loc_here(2))=EOF_1(i);
end

%% mean states during MSD

% precip
load('daily_sst_noaa');
load('daily_cpc_precip');
load('msd_collection_ncep');
time_used=datevec(datenum(1979,1,1):datenum(2017,12,31));

m_d_unique=unique(time_used(:,2:3),'rows');

temp=double(ncread('temp_low_1979_2017.nc','t2m'));
temp(temp==-32767)=nan;
temp=temp-273.15;
u_wind=double(ncread('u_low_1979_2017.nc','u10'));
u_wind(u_wind==-32767)=nan;
v_wind=double(ncread('v_low_1979_2017.nc','v10'));
v_wind(v_wind==-32767)=nan;
pres=double(ncread('pres_low_1979_2017.nc','sp'));
pres(pres==-32767)=nan;
% sst=double(ncread('sst_low_1979_2017.nc','sst'));
% sst(sst==-32767)=nan;
% sst=sst-273.15;
% sst(sst<5)=nan;

precip_anom=NaN(size(precip));
temp_anom=NaN(size(temp));
u_anom=NaN(size(u_wind));
v_anom=NaN(size(v_wind));
sst_anom=NaN(size(sst));
pres_anom=NaN(size(pres));

for i=1:size(m_d_unique,1);
    m_d_here=m_d_unique(i,:);
    index_here=find(time_used(:,2)==m_d_here(1) & time_used(:,3)==m_d_here(2));
    %precip_anom(:,:,index_here)=precip(:,:,index_here)-nanmean(precip(:,:,index_here),3);
    %temp_anom(:,:,index_here)=temp(:,:,index_here)-nanmean(temp(:,:,index_here),3);
    %u_anom(:,:,index_here)=u_wind(:,:,index_here)-nanmean(u_wind(:,:,index_here),3);
    %v_anom(:,:,index_here)=v_wind(:,:,index_here)-nanmean(v_wind(:,:,index_here),3);
    sst_anom(:,:,index_here)=sst(:,:,index_here)-nanmean(sst(:,:,index_here),3);
    %pres_anom(:,:,index_here)=pres(:,:,index_here)-nanmean(pres(:,:,index_here),3);
end

start_full=[];
end_full=[];
peak_full=[];

for i=1:size(msd_collection_ncep,1);
    start_full=[start_full (msd_collection_ncep(i,4))-datenum(1979,1,1)+1];
    end_full=[end_full (msd_collection_ncep(i,6))-datenum(1979,1,1)+1];
    peak_full=[peak_full (msd_collection_ncep(i,8))-datenum(1979,1,1)+1];
end
 
precip_anom_start = nanmean(precip_anom(:,:,start_full),3);
precip_anom_end = nanmean(precip_anom(:,:,end_full),3);
precip_anom_peak = nanmean(precip_anom(:,:,peak_full),3);

temp_anom_start = nanmean(temp_anom(:,:,start_full),3);
temp_anom_end = nanmean(temp_anom(:,:,end_full),3);
temp_anom_peak = nanmean(temp_anom(:,:,peak_full),3);

u_anom_start = nanmean(u_anom(:,:,start_full),3);
u_anom_end = nanmean(u_anom(:,:,end_full),3);
u_anom_peak = nanmean(u_anom(:,:,peak_full),3);

v_anom_start = nanmean(v_anom(:,:,start_full),3);
v_anom_end = nanmean(v_anom(:,:,end_full),3);
v_anom_peak = nanmean(v_anom(:,:,peak_full),3);

sst_anom_start = nanmean(sst_anom(:,:,start_full),3);
sst_anom_end = nanmean(sst_anom(:,:,end_full),3);
sst_anom_peak = nanmean(sst_anom(:,:,peak_full),3);

pres_anom_start = nanmean(pres_anom(:,:,start_full),3);
pres_anom_end = nanmean(pres_anom(:,:,end_full),3);
pres_anom_peak = nanmean(pres_anom(:,:,peak_full),3);



% % enso
% e_year=[1980 1983 1987 1988 1992 1995 1998 2003 2005 2007 2010 2015 2016];
% l_year=[1984 1985 1989 1996 1999 2000 2001 2006 2008 2009 2011 2012 2017];
% 
% start_full_e=[];
% start_full_l=[];
% start_full_n=[];
% 
% for i=1:size(msd_collection_ncep,1);
%     msd_here=msd_collection_ncep(i,:);
%     if ismember(msd_here(1),e_year)
%         start_full_e=[start_full_e;msd_here(4)-datenum(1979,1,1)+1];
%     elseif ismember(msd_here(1),l_year)
%         start_full_l=[start_full_l;msd_here(4)-datenum(1979,1,1)+1];
%     else
%         start_full_n=[start_full_n;msd_here(4)-datenum(1979,1,1)+1];
%     end
% end
% 
% precip_anom_mean_e = nanmean(precip_anom(:,:,start_full_e),3);
% precip_anom_mean_l = nanmean(precip_anom(:,:,start_full_l),3);
% precip_anom_mean_n = nanmean(precip_anom(:,:,start_full_n),3);
% 
% 
% imsd_e=nanmean(msd_collection_ncep(ismember(msd_collection_ncep(:,1),e_year),12));
% imsd_l=nanmean(msd_collection_ncep(ismember(msd_collection_ncep(:,1),l_year),12));

%% annual MSD numbers
msd_num=[];
for i=1979:2017
    msd_here=msd_collection_ncep(msd_collection_ncep(:,1)==i,:);
    msd_num=[msd_num;size(msd_here,1)];
end

sst_start = nanmean(sst(:,:,start_full),3);
sst_end = nanmean(sst(:,:,end_full),3);
sst_peak = nanmean(sst(:,:,peak_full),3);

m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
figure
m_contourf(lon_era_low,lat_era_low,sst_june',25:0.1:30,'linestyle','none');
colormap(jet);
m_coast();
m_grid();

figure
m_contourf(lon_era_low,lat_era_low,sst_julaug',25:0.1:30,'linestyle','none');
colormap(jet);
m_coast();
m_grid();

figure
m_contourf(lon_era_low,lat_era_low,sst_sep',25:0.1:30,'linestyle','none');
colormap(jet);
m_coast();
m_grid();

date_full=datevec(datenum(1979,1,1):datenum(2017,12,31));

sst_june=nanmean(sst(:,:,date_full(:,2)==6),3);
sst_julaug=nanmean(sst(:,:,date_full(:,2)==7 | date_full(:,2)==8),3);
sst_sep=nanmean(sst(:,:,date_full(:,2)==9),3);



%% try enso stuff again
oni=csvread('ENSO_year_index.csv');
oni=oni(oni(:,1)>=1979 & oni(:,1)<=2017,:);

oni_mean=[];

for i=1:39;
    oni_here=oni(i,2:13);
    oni_mean=[oni_mean;nanmean(oni_here(9:12))];
end
full_start=[];
full_peak=[];
full_end=[];
full_dur=[];
full_imsd=[];
for i=1979:2017;
    start_here=msd_collection_ncep(msd_collection_ncep(:,1)==i,5);
    start_here=nanmean(start_here);
    full_start=[full_start;start_here];
    
    end_here=msd_collection_ncep(msd_collection_ncep(:,1)==i,7);
    end_here=nanmean(end_here);
    full_end=[full_end;end_here];
    
    peak_here=msd_collection_ncep(msd_collection_ncep(:,1)==i,9);
    peak_here=nanmean(peak_here);
    full_peak=[full_peak;peak_here];
    
    dur_here=msd_collection_ncep(msd_collection_ncep(:,1)==i,7)-msd_collection_ncep(msd_collection_ncep(:,1)==i,5)+1;
    dur_here=nanmean(dur_here);
    full_dur=[full_dur;dur_here];
    
    imsd_here=msd_collection_ncep(msd_collection_ncep(:,1)==i,end);
    imsd_here=nanmean(imsd_here);
    full_imsd=[full_imsd;imsd_here];
end

oni_msd=[];

for i=1:size(msd_collection_ncep,1);
    msd_here=msd_collection_ncep(i,:);
    oni_here=oni_mean(msd_here(1)-1979+1);
    oni_msd=[oni_msd;oni_here];
end

tbl_for_cor=table(full_start,full_end,full_peak,full_dur,full_imsd,oni_mean,'variablenames',{'Onset','End','Peak','Duration','Imsd','Oni'});




corrplot(tbl_for_cor,'testR','on');
title('Correlation plot','fontsize',16);

%% sst? may june-july sep
load('all_sst_noaa','sst');
time_used=datevec(datenum(1979,1,1):datenum(2017,12,31));
sst_j=double(nanmean(sst(:,:,time_used(:,2)==6),3));
sst_ja=double(nanmean(sst(:,:,time_used(:,2)==7 | time_used(:,2)==8),3));
sst_s=double(nanmean(sst(:,:,time_used(:,2)==9),3));

% hovmoller diagram
date_used=datevec(datenum(1979,5,15):datenum(1979,10,15));
date_used=date_used(:,2:3);
sst_hm=NaN(240,120,size(date_used,1));

for i=1:size(date_used,1);
    time_here=find(time_used(:,2)==date_used(i,1) & time_used(:,3)==date_used(i,2));
    sst_here=nanmean(sst(:,:,time_here),3);
    sst_hm(:,:,i)=sst_here;
end

sst_hm=sst_hm(lon_sst<=270 & lon_sst>=258,lat_sst>=12 & lat_sst<=18,:);
sst_hm=squeeze(nanmean(sst_hm,2));

%% precip in El/La/neutral 

en=[1980 1983 1987 1988 1992 1995 1998 2003 2005 2007 2010 2015 2016]-1;
la=[1984 1985 1989 1996 1999 2000 2001 2006 2008 2009 2011 2012 2017 2018]-1;
year_used=1979:2017;
ne=year_used(~ismember(year_used,en) & ~ismember(year_used,la));

date_used=datevec(datenum(1979,1,1):datenum(2017,12,31));

load('daily_anom','precip_anom');

precip_sep_en=nanmean(precip_anom(:,:,date_used(:,2)==9 & ismember(date_used(:,1),en)),3);
length_here=nansum(date_used(:,2)==9 & ismember(date_used(:,1),en));
data_here=mat2cell(precip_anom(:,:,date_used(:,2)==9 & ismember(date_used(:,1),en)),ones(120,1),ones(60,1),length_here);
see_sep_en=cellfun(@ttest,data_here);

precip_sep_la=nanmean(precip_anom(:,:,date_used(:,2)==9 & ismember(date_used(:,1),la)),3);
length_here=nansum(date_used(:,2)==9 & ismember(date_used(:,1),la));
data_here=mat2cell(precip_anom(:,:,date_used(:,2)==9 & ismember(date_used(:,1),la)),ones(120,1),ones(60,1),length_here);
see_sep_la=cellfun(@ttest,data_here);

precip_sep_ne=nanmean(precip_anom(:,:,date_used(:,2)==9 & ismember(date_used(:,1),ne)),3);
length_here=nansum(date_used(:,2)==9 & ismember(date_used(:,1),ne));
data_here=mat2cell(precip_anom(:,:,date_used(:,2)==9 & ismember(date_used(:,1),ne)),ones(120,1),ones(60,1),length_here);
see_sep_ne=cellfun(@ttest,data_here);

precip_jun_en=nanmean(precip_anom(:,:,date_used(:,2)==6 & ismember(date_used(:,1),en)),3);
length_here=nansum(date_used(:,2)==6 & ismember(date_used(:,1),en));
data_here=mat2cell(precip_anom(:,:,date_used(:,2)==6 & ismember(date_used(:,1),en)),ones(120,1),ones(60,1),length_here);
see_jun_en=cellfun(@ttest,data_here);

precip_jun_la=nanmean(precip_anom(:,:,date_used(:,2)==6 & ismember(date_used(:,1),la)),3);
length_here=nansum(date_used(:,2)==6 & ismember(date_used(:,1),la));
data_here=mat2cell(precip_anom(:,:,date_used(:,2)==6 & ismember(date_used(:,1),la)),ones(120,1),ones(60,1),length_here);
see_jun_la=cellfun(@ttest,data_here);

precip_jun_ne=nanmean(precip_anom(:,:,date_used(:,2)==6 & ismember(date_used(:,1),ne)),3);
length_here=nansum(date_used(:,2)==6 & ismember(date_used(:,1),ne));
data_here=mat2cell(precip_anom(:,:,date_used(:,2)==6 & ismember(date_used(:,1),ne)),ones(120,1),ones(60,1),length_here);
see_jun_ne=cellfun(@ttest,data_here);

precip_ja_en=nanmean(precip_anom(:,:,(date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),en)),3);
length_here=nansum((date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),en));
data_here=mat2cell(precip_anom(:,:,(date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),en)),ones(120,1),ones(60,1),length_here);
see_ja_en=cellfun(@ttest,data_here);

precip_ja_la=nanmean(precip_anom(:,:,(date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),la)),3);
length_here=nansum((date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),la));
data_here=mat2cell(precip_anom(:,:,(date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),la)),ones(120,1),ones(60,1),length_here);
see_ja_la=cellfun(@ttest,data_here);

precip_ja_ne=nanmean(precip_anom(:,:,(date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),ne)),3);
length_here=nansum((date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),ne));
data_here=mat2cell(precip_anom(:,:,(date_used(:,2)==7 | date_used(:,2)==8) & ismember(date_used(:,1),ne)),ones(120,1),ones(60,1),length_here);
see_ja_ne=cellfun(@ttest,data_here);

save precip_enso_new *_ne *_la *_en


figure();
m_contourf(lon_ncep,lat_ncep,precip_sep_en',-4:0.01:4,'linestyle','none');
colormap(jet);
colormap(color_here);
m_coast();
caxis([-3 3]);

figure();
m_contourf(lon_ncep,lat_ncep,precip_sep_la',-4:0.01:4,'linestyle','none');
colormap(jet);
colormap(color_here);
m_coast();
caxis([-3 3]);

figure();
m_contourf(lon_ncep,lat_ncep,precip_sep_ne',-4:0.01:4,'linestyle','none');
m_coast();
colormap(jet);
colormap(color_here);
caxis([-3 3]);

%% precip west???

cd('/Volumes/mydrive/msd_data/rainfall');

lon=double(ncread('precip.1999.nc','lon'));
lat=double(ncread('precip.1999.nc','lat'));

lon_used=find(lon>=100 & lon<=120);
lat_used=find(lat<=30 & lat>=0);

precip=[];

for i=1979:2017
    file_here=['precip.' num2str(i) '.nc'];
    precip_here=ncread(file_here,'precip',[lon_used(1) lat_used(1) 1],[length(lon_used),length(lat_used),inf]);
    precip=cat(3,precip,precip_here);
end

precip=double(precip);

precip(abs(precip)>1000)=nan;

time_used=datevec(datenum(1979,1,1):datenum(2017,12,31));

time_used=time_used(:,1:3);

y_m=unique(time_used(:,2:3),'rows');

precip_west_anom=NaN(size(precip));

for i=1:size(y_m,1);
    index_here=find(time_used(:,2)==y_m(i,1) & time_used(:,3)==y_m(i,2));
    precip_west_anom(:,:,index_here)=precip(:,:,index_here)-nanmean(precip(:,:,index_here),3);
end

precip_west_start=NaN(40,60);
precip_west_peak=NaN(40,60);
precip_west_end=NaN(40,60);

start_index=msd_collection_ncep(:,4)-datenum(1979,1,1)+1;
peak_index=msd_collection_ncep(:,8)-datenum(1979,1,1)+1;
end_index=msd_collection_ncep(:,6)-datenum(1979,1,1)+1;

precip_west_start=nanmean(precip_west_anom(:,:,start_index),3);
precip_west_peak=nanmean(precip_west_anom(:,:,peak_index),3);
precip_west_end=nanmean(precip_west_anom(:,:,end_index),3);

figure
m_contourf(lon(lon_used),lat(lat_used),precip_west_start',-10:0.01:10,'linestyle','none');
colormap(color_here);
caxis([-4 4]);
m_coast();
m_grid;

figure
m_contourf(lon(lon_used),lat(lat_used),precip_west_peak',-10:0.01:10,'linestyle','none');
colormap(color_here);
caxis([-4 4]);
m_coast();
m_grid;


figure
m_contourf(lon(lon_used),lat(lat_used),precip_west_end',-10:0.01:10,'linestyle','none');
colormap(color_here);
caxis([-4 4]);
m_coast();
m_grid;

%% hovo equator
%u
addpath('/Volumes/mydrive/phd_project/era_high_scale');
u10_full=ncread('u10_msd_hovo_1993_2017.nc','u10');
u10_full(abs(u10_full)>10000)=nan;
t=ncread('u10_msd_hovo_1993_2017.nc','time');

t=double(double(t)./24);
t=datevec(t+datenum(1900,1,1)-3./24);

unique_t=unique(t(:,1:3),'rows');
u10_pacific_daily=NaN(121,31,size(unique_t,1));

for i=1:size(unique_t,1);
    index_here=find(t(:,1)==unique_t(i,1) & t(:,2)==unique_t(i,2) & t(:,3)==unique_t(i,3));
    u10_pacific_daily(:,:,i)=nanmean(u10_full(:,:,index_here),3);
end

y_m=unique(unique_t(:,2:3),'rows');

u10_pacific_anom=NaN(size(u10_pacific_daily));

for i=1:size(y_m,1);
    index_here=find(unique_t(:,2)==y_m(i,1) & unique_t(:,3)==y_m(i,2));
    u10_pacific_anom(:,:,index_here)=u10_pacific_daily(:,:,index_here)-nanmean(u10_pacific_daily(:,:,index_here),3);
end

%v
addpath('/Volumes/mydrive/phd_project/era_high_scale');
v10_full=ncread('v10_msd_hovo_1993_2017.nc','v10');
v10_full(abs(v10_full)>10000)=nan;
t=ncread('v10_msd_hovo_1993_2017.nc','time');

t=double(double(t)./24);
t=datevec(t+datenum(1900,1,1)-3./24);

unique_t=unique(t(:,1:3),'rows');
v10_pacific_daily=NaN(121,31,size(unique_t,1));

for i=1:size(unique_t,1);
    index_here=find(t(:,1)==unique_t(i,1) & t(:,2)==unique_t(i,2) & t(:,3)==unique_t(i,3));
    v10_pacific_daily(:,:,i)=nanmean(v10_full(:,:,index_here),3);
end

y_m=unique(unique_t(:,2:3),'rows');

v10_pacific_anom=NaN(size(v10_pacific_daily));

for i=1:size(y_m,1);
    index_here=find(unique_t(:,2)==y_m(i,1) & unique_t(:,3)==y_m(i,2));
    v10_pacific_anom(:,:,index_here)=v10_pacific_daily(:,:,index_here)-nanmean(v10_pacific_daily(:,:,index_here),3);
end


start_full=msd_collection_ncep(:,4)-datenum(1979,1,1)+1;
peak_full=msd_collection_ncep(:,8)-datenum(1979,1,1)+1;
end_full=msd_collection_ncep(:,6)-datenum(1979,1,1)+1;

u10_pacific_start=nanmean(u10_pacific_anom(:,:,start_full),3);
v10_pacific_start=nanmean(v10_pacific_anom(:,:,start_full),3);

u10_pacific_peak=nanmean(u10_pacific_anom(:,:,peak_full),3);
v10_pacific_peak=nanmean(v10_pacific_anom(:,:,peak_full),3);

u10_pacific_end=nanmean(u10_pacific_anom(:,:,end_full),3);
v10_pacific_end=nanmean(v10_pacific_anom(:,:,end_full),3);


% u850
addpath('/Volumes/mydrive/phd_project/era_high_scale');
u850_full=ncread('u850_msd_hovo_1979_2017.nc','u');
u850_full(abs(u850_full)>10000)=nan;
t=ncread('u850_msd_hovo_1979_2017.nc','time');

t=double(double(t)./24);
t=datevec(t+datenum(1900,1,1));

unique_t=unique(t(:,1:3),'rows');
u850_pacific_daily=NaN(121,31,size(unique_t,1));

for i=1:size(unique_t,1);
    index_here=find(t(:,1)==unique_t(i,1) & t(:,2)==unique_t(i,2) & t(:,3)==unique_t(i,3));
    u850_pacific_daily(:,:,i)=nanmean(u850_full(:,:,index_here),3);
end

y_m=unique(unique_t(:,2:3),'rows');

u850_pacific_anom=NaN(size(u850_pacific_daily));

for i=1:size(y_m,1);
    index_here=find(unique_t(:,2)==y_m(i,1) & unique_t(:,3)==y_m(i,2));
    u850_pacific_anom(:,:,index_here)=u850_pacific_daily(:,:,index_here)-nanmean(u850_pacific_daily(:,:,index_here),3);
end

% v850
addpath('/Volumes/mydrive/phd_project/era_high_scale');
v850_full=ncread('v850_msd_hovo_1979_2017.nc','v');
v850_full(abs(v850_full)>10000)=nan;
t=ncread('v850_msd_hovo_1979_2017.nc','time');

t=double(double(t)./24);
t=datevec(t+datenum(1900,1,1));

unique_t=unique(t(:,1:3),'rows');
v850_pacific_daily=NaN(121,31,size(unique_t,1));

for i=1:size(unique_t,1);
    index_here=find(t(:,1)==unique_t(i,1) & t(:,2)==unique_t(i,2) & t(:,3)==unique_t(i,3));
    v850_pacific_daily(:,:,i)=nanmean(v850_full(:,:,index_here),3);
end

y_m=unique(unique_t(:,2:3),'rows');

v850_pacific_anom=NaN(size(v850_pacific_daily));

for i=1:size(y_m,1);
    index_here=find(unique_t(:,2)==y_m(i,1) & unique_t(:,3)==y_m(i,2));
    v850_pacific_anom(:,:,index_here)=v850_pacific_daily(:,:,index_here)-nanmean(v850_pacific_daily(:,:,index_here),3);
end

lon=double(ncread('v10_msd_hovo_1993_2017.nc','longitude'));
lat=double(ncread('v10_msd_hovo_1993_2017.nc','latitude'));
[lon_full,lat_full]=meshgrid(lon,lat);

start_full=msd_collection_ncep(:,4)-datenum(1979,1,1)+1;
peak_full=msd_collection_ncep(:,8)-datenum(1979,1,1)+1;
end_full=msd_collection_ncep(:,6)-datenum(1979,1,1)+1;

u850_pacific_start=nanmean(u850_pacific_daily(:,:,start_full),3);
v850_pacific_start=nanmean(v850_pacific_daily(:,:,start_full),3);

u850_pacific_peak=nanmean(u850_pacific_daily(:,:,peak_full),3);
v850_pacific_peak=nanmean(v850_pacific_daily(:,:,peak_full),3);

u850_pacific_end=nanmean(u850_pacific_daily(:,:,end_full),3);
v850_pacific_end=nanmean(v850_pacific_daily(:,:,end_full),3);



figure
m_quiver(lon_full(1:5:end,1:5:end),lat_full(1:5:end,1:5:end),(u850_pacific_start(1:5:end,1:5:end))',(v850_pacific_start(1:5:end,1:5:end))',...
    2,'autoscale','off');
m_coast();
m_grid;

figure
m_quiver(lon_full(1:5:end,1:5:end),lat_full(1:5:end,1:5:end),(u850_pacific_peak(1:5:end,1:5:end))',(v850_pacific_peak(1:5:end,1:5:end))',...
    2,'autoscale','off');
m_coast();
m_grid;

figure
m_quiver(lon_full(1:5:end,1:5:end),lat_full(1:5:end,1:5:end),(u850_pacific_end(1:5:end,1:5:end))',(v850_pacific_end(1:5:end,1:5:end))',...
    2,'autoscale','off');
m_coast();
m_grid;

%% ENSO patterns
en=[1980 1983 1987 1988 1992 1995 1998 2003 2005 2007 2010 2015 2016]-1;
la=[1984 1985 1989 1996 1999 2000 2001 2006 2008 2009 2011 2012 2017 2018]-1;
year_used=1979:2017;
ne=year_used(~ismember(year_used,en) & ~ismember(year_used,la));

date_used=datevec(datenum(1979,1,1):datenum(2017,12,31));
load('msd_collection_ncep');
year_full=msd_collection_ncep(:,1);
start_full=msd_collection_ncep(:,4)-datenum(1979,1,1)+1;
end_full=msd_collection_ncep(:,6)-datenum(1979,1,1)+1;
peak_full=msd_collection_ncep(:,8)-datenum(1979,1,1)+1;

start_en=start_full(ismember(year_full,en));
end_en=end_full(ismember(year_full,en));
peak_en=peak_full(ismember(year_full,en));

start_la=start_full(ismember(year_full,la));
end_la=end_full(ismember(year_full,la));
peak_la=peak_full(ismember(year_full,la));

start_ne=start_full(ismember(year_full,ne));
end_ne=end_full(ismember(year_full,ne));
peak_ne=peak_full(ismember(year_full,ne));

load('daily_anom','pres_anom','u_anom','v_anom');

% start
pres_start_en=nanmean(pres_anom(:,:,start_en),3);
data_here=mat2cell(pres_anom(:,:,start_en),ones(121,1),ones(61,1),length(start_en));
see_start_en=cellfun(@ttest,data_here);
u_start_en=nanmean(u_anom(:,:,start_en),3);
v_start_en=nanmean(v_anom(:,:,start_en),3);

pres_start_la=nanmean(pres_anom(:,:,start_la),3);
data_here=mat2cell(pres_anom(:,:,start_la),ones(121,1),ones(61,1),length(start_la));
see_start_la=cellfun(@ttest,data_here);
u_start_la=nanmean(u_anom(:,:,start_la),3);
v_start_la=nanmean(v_anom(:,:,start_la),3);

pres_start_ne=nanmean(pres_anom(:,:,start_ne),3);
data_here=mat2cell(pres_anom(:,:,start_ne),ones(121,1),ones(61,1),length(start_ne));
see_start_ne=cellfun(@ttest,data_here);
u_start_ne=nanmean(u_anom(:,:,start_ne),3);
v_start_ne=nanmean(v_anom(:,:,start_ne),3);

% peak
pres_peak_en=nanmean(pres_anom(:,:,peak_en),3);
data_here=mat2cell(pres_anom(:,:,peak_en),ones(121,1),ones(61,1),length(peak_en));
see_peak_en=cellfun(@ttest,data_here);
u_peak_en=nanmean(u_anom(:,:,peak_en),3);
v_peak_en=nanmean(v_anom(:,:,peak_en),3);

pres_peak_la=nanmean(pres_anom(:,:,peak_la),3);
data_here=mat2cell(pres_anom(:,:,peak_la),ones(121,1),ones(61,1),length(peak_la));
see_peak_la=cellfun(@ttest,data_here);
u_peak_la=nanmean(u_anom(:,:,peak_la),3);
v_peak_la=nanmean(v_anom(:,:,peak_la),3);

pres_peak_ne=nanmean(pres_anom(:,:,peak_ne),3);
data_here=mat2cell(pres_anom(:,:,peak_ne),ones(121,1),ones(61,1),length(peak_ne));
see_peak_ne=cellfun(@ttest,data_here);
u_peak_ne=nanmean(u_anom(:,:,peak_ne),3);
v_peak_ne=nanmean(v_anom(:,:,peak_ne),3);

% end
pres_end_en=nanmean(pres_anom(:,:,end_en),3);
data_here=mat2cell(pres_anom(:,:,end_en),ones(121,1),ones(61,1),length(end_en));
see_end_en=cellfun(@ttest,data_here);
u_end_en=nanmean(u_anom(:,:,end_en),3);
v_end_en=nanmean(v_anom(:,:,end_en),3);

pres_end_la=nanmean(pres_anom(:,:,end_la),3);
data_here=mat2cell(pres_anom(:,:,end_la),ones(121,1),ones(61,1),length(end_la));
see_end_la=cellfun(@ttest,data_here);
u_end_la=nanmean(u_anom(:,:,end_la),3);
v_end_la=nanmean(v_anom(:,:,end_la),3);

pres_end_ne=nanmean(pres_anom(:,:,end_ne),3);
data_here=mat2cell(pres_anom(:,:,end_ne),ones(121,1),ones(61,1),length(end_ne));
see_end_ne=cellfun(@ttest,data_here);
u_end_ne=nanmean(u_anom(:,:,end_ne),3);
v_end_ne=nanmean(v_anom(:,:,end_ne),3);

save enso_patterns_new *_end_* *_peak_* *_start_*

%% MJO phase?
% classification for MSD
load('msd_collection_ncep.mat')
start_end_phase_4=NaN(21872,4,2);

for i=1:size(msd_collection_ncep,1);
    msd_here=msd_collection_ncep(i,:);
    start_here=msd_here(4);
    end_here=msd_here(6);
    peak_here=msd_here(8);
    
    start_end_phase_4(i,1,1)=start_here;
    start_end_phase_4(i,1,2)=round(quantile([start_here peak_here],0.3));
    start_end_phase_4(i,2,1)=round(quantile([start_here peak_here],0.3))+1;
    start_end_phase_4(i,2,2)=peak_here;
    start_end_phase_4(i,3,1)=peak_here+1;
    start_end_phase_4(i,3,2)=round(quantile([peak_here+1 end_here],0.7));
    start_end_phase_4(i,4,1)=round(quantile([peak_here+1 end_here],0.7))+1;
    start_end_phase_4(i,4,2)=end_here;
end

% phase index

phase_1_index=[];
phase_2_index=[];
phase_3_index=[];
phase_4_index=[];

for i=1:size(msd_collection_ncep,1);
    phase_1_index=[phase_1_index start_end_phase_4(i,1,1):start_end_phase_4(i,1,2)];
    phase_2_index=[phase_2_index start_end_phase_4(i,2,1):start_end_phase_4(i,2,2)];
    phase_3_index=[phase_3_index start_end_phase_4(i,3,1):start_end_phase_4(i,3,2)];
    phase_4_index=[phase_4_index start_end_phase_4(i,4,1):start_end_phase_4(i,4,2)];
end

% patterns in different phases
load('daily_anom','pres_anom','u_anom','v_anom');
event_index_1=zeros(14245,1);
event_index_2=zeros(14245,1);
event_index_3=zeros(14245,1);
event_index_4=zeros(14245,1);

for i=1:14245;
    event_index_1(i)=nansum((phase_1_index-datenum(1979,1,1)+1)==i)./length(phase_1_index);
    event_index_2(i)=nansum((phase_2_index-datenum(1979,1,1)+1)==i)./length(phase_2_index);
    event_index_3(i)=nansum((phase_3_index-datenum(1979,1,1)+1)==i)./length(phase_3_index);
    event_index_4(i)=nansum((phase_4_index-datenum(1979,1,1)+1)==i)./length(phase_4_index);
end

pres_msd_1=nansum(pres_anom.*repmat(reshape(event_index_1,1,1,14245),121,61,1),3);
u_msd_1=nansum(u_anom.*repmat(reshape(event_index_1,1,1,14245),121,61,1),3);
v_msd_1=nansum(v_anom.*repmat(reshape(event_index_1,1,1,14245),121,61,1),3);

pres_msd_2=nansum(pres_anom.*repmat(reshape(event_index_2,1,1,14245),121,61,1),3);
u_msd_2=nansum(u_anom.*repmat(reshape(event_index_2,1,1,14245),121,61,1),3);
v_msd_2=nansum(v_anom.*repmat(reshape(event_index_2,1,1,14245),121,61,1),3);

pres_msd_3=nansum(pres_anom.*repmat(reshape(event_index_3,1,1,14245),121,61,1),3);
u_msd_3=nansum(u_anom.*repmat(reshape(event_index_3,1,1,14245),121,61,1),3);
v_msd_3=nansum(v_anom.*repmat(reshape(event_index_3,1,1,14245),121,61,1),3);

pres_msd_4=nansum(pres_anom.*repmat(reshape(event_index_4,1,1,14245),121,61,1),3);
u_msd_4=nansum(u_anom.*repmat(reshape(event_index_4,1,1,14245),121,61,1),3);
v_msd_4=nansum(v_anom.*repmat(reshape(event_index_4,1,1,14245),121,61,1),3);

% mjo phase patterns
addpath('/Volumes/mydrive/msd_data/mjo_data');
load('daily_anom','pres_anom','u_anom','v_anom');
date_used=datenum(1974,6,1):datenum(2018,8,30);
rmm1=double(ncread('rmm1.nc','RMM1'));
rmm2=double(ncread('rmm2.nc','RMM2'));
amp=double(ncread('amp.nc','amplitude'));
p=double(ncread('phase.nc','phase'));

rmm1=rmm1((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);
rmm2=rmm2((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);
amp=amp((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);
p=p((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);
    
date_used=datenum(1979,1,1):datenum(2017,12,31);
date_used=datevec(date_used);

index_8_1 = find(ismember(date_used(:,2),6:9) & (ismember(p,[1 8])));
index_2_3 = find(ismember(date_used(:,2),6:9) & (ismember(p,[2 3])));
index_4_5 = find(ismember(date_used(:,2),6:9) & (ismember(p,[4 5])));
index_6_7 = find(ismember(date_used(:,2),6:9) & (ismember(p,[6 7])));

pres_mjo_81=nanmean(pres_anom(:,:,index_8_1),3);
u_mjo_81=nanmean(u_anom(:,:,index_8_1),3);
v_mjo_81=nanmean(v_anom(:,:,index_8_1),3);

pres_mjo_23=nanmean(pres_anom(:,:,index_2_3),3);
u_mjo_23=nanmean(u_anom(:,:,index_2_3),3);
v_mjo_23=nanmean(v_anom(:,:,index_2_3),3);

pres_mjo_45=nanmean(pres_anom(:,:,index_4_5),3);
u_mjo_45=nanmean(u_anom(:,:,index_4_5),3);
v_mjo_45=nanmean(v_anom(:,:,index_4_5),3);

pres_mjo_67=nanmean(pres_anom(:,:,index_6_7),3);
u_mjo_67=nanmean(u_anom(:,:,index_6_7),3);
v_mjo_67=nanmean(v_anom(:,:,index_6_7),3);


%map proportion
prop_1=NaN(120,60,4);
prop_2=NaN(120,60,4);
prop_3=NaN(120,60,4);
prop_4=NaN(120,60,4);
date_used=datenum(1979,1,1):datenum(2017,12,31);
loc_unique=unique(msd_collection_ncep(:,2:3),'rows');
for i=1:size(loc_unique,1);
    loc_here=loc_unique(i,:);
    index_here = msd_collection_ncep(:,2)==loc_here(1) & msd_collection_ncep(:,3)==loc_here(2);
    
    start_end_here=squeeze(start_end_phase_4(index_here,:,:));
    
    phase_1=[];
    phase_2=[];
    phase_3=[];
    phase_4=[];
    
    for j=1:size(start_end_here,1);
        
        phase_1_here=squeeze(start_end_here(j,1,:));
        phase_1_here=phase_1_here(1):phase_1_here(2);
        phase_1=[phase_1 phase_1_here];
        
        phase_2_here=squeeze(start_end_here(j,2,:));
        phase_2_here=phase_2_here(1):phase_2_here(2);
        phase_2=[phase_2 phase_2_here];
        
        phase_3_here=squeeze(start_end_here(j,3,:));
        phase_3_here=phase_3_here(1):phase_3_here(2);
        phase_3=[phase_3 phase_3_here];
        
        phase_4_here=squeeze(start_end_here(j,4,:));
        phase_4_here=phase_4_here(1):phase_4_here(2);
        phase_4=[phase_4 phase_4_here];
        
    end
    
    prop_1(loc_here(1),loc_here(2),1)=nansum(ismember(phase_1,date_used(p==1 | p==8)))./length(phase_1);
    prop_1(loc_here(1),loc_here(2),2)=nansum(ismember(phase_1,date_used(p==2 | p==3)))./length(phase_1);
    prop_1(loc_here(1),loc_here(2),3)=nansum(ismember(phase_1,date_used(p==4 | p==5)))./length(phase_1);
    prop_1(loc_here(1),loc_here(2),4)=nansum(ismember(phase_1,date_used(p==6 | p==7)))./length(phase_1);
    
    prop_2(loc_here(1),loc_here(2),1)=nansum(ismember(phase_2,date_used(p==1 | p==8)))./length(phase_2);
    prop_2(loc_here(1),loc_here(2),2)=nansum(ismember(phase_2,date_used(p==2 | p==3)))./length(phase_2);
    prop_2(loc_here(1),loc_here(2),3)=nansum(ismember(phase_2,date_used(p==4 | p==5)))./length(phase_2);
    prop_2(loc_here(1),loc_here(2),4)=nansum(ismember(phase_2,date_used(p==6 | p==7)))./length(phase_2);
    
    prop_3(loc_here(1),loc_here(2),1)=nansum(ismember(phase_3,date_used(p==1 | p==8)))./length(phase_3);
    prop_3(loc_here(1),loc_here(2),2)=nansum(ismember(phase_3,date_used(p==2 | p==3)))./length(phase_3);
    prop_3(loc_here(1),loc_here(2),3)=nansum(ismember(phase_3,date_used(p==4 | p==5)))./length(phase_3);
    prop_3(loc_here(1),loc_here(2),4)=nansum(ismember(phase_3,date_used(p==6 | p==7)))./length(phase_3);
    
    prop_4(loc_here(1),loc_here(2),1)=nansum(ismember(phase_4,date_used(p==1 | p==8)))./length(phase_4);
    prop_4(loc_here(1),loc_here(2),2)=nansum(ismember(phase_4,date_used(p==2 | p==3)))./length(phase_4);
    prop_4(loc_here(1),loc_here(2),3)=nansum(ismember(phase_4,date_used(p==4 | p==5)))./length(phase_4);
    prop_4(loc_here(1),loc_here(2),4)=nansum(ismember(phase_4,date_used(p==6 | p==7)))./length(phase_4);
    
end
    
    
    
figure
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,1))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,2))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,3))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_1(:,:,4))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar


figure
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,1))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,2))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,3))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_3(:,:,4))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar


figure
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,1))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,2))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,3))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

figure
m_contourf(lon_ncep,lat_ncep,(prop_4(:,:,4))',0:0.01:1,'linestyle','none');
colormap(jet);
caxis([0 0.6]);
colorbar

% draw example?
load('precip_clim');
example_ts=smoothdata(squeeze(precip_clim(45,18,:)),1,'gaussian',31);
period_1=150:200;
period_2=250:300;

t1=period_1(find(example_ts(period_1)==nanmax(example_ts(period_1))));
t2=period_2(find(example_ts(period_2)==nanmax(example_ts(period_2))));

period_p=t1:t2;
tmin=period_p(find(example_ts(t1:t2)==nanmin(example_ts(t1:t2))));
figure('pos',[10 10 1000 1000]);
plot(1:366,example_ts,'linewidth',2);
phase_1=t1:round(quantile([t1 tmin],0.35));
phase_2=round(quantile([t1 tmin],0.35)):tmin;
phase_3=tmin:round(quantile([tmin t2],0.65));
phase_4=round(quantile([tmin t2],0.65)):t2;
hold on
x1_phase1=phase_1;
y1_phase1=zeros(1,length(x1_phase1));
x2_phase1=phase_1(end:-1:1);
y2_phase1=example_ts(x2_phase1);
h1=fill([x1_phase1(:);x2_phase1(:)],[y1_phase1(:);y2_phase1(:)],'red','linestyle','none');
hold on
x1_phase1=phase_2;
y1_phase1=zeros(1,length(x1_phase1));
x2_phase1=phase_2(end:-1:1);
y2_phase1=example_ts(x2_phase1);
h2=fill([x1_phase1(:);x2_phase1(:)],[y1_phase1(:);y2_phase1(:)],'yellow','linestyle','none');
hold on
x1_phase1=phase_3;
y1_phase1=zeros(1,length(x1_phase1));
x2_phase1=phase_3(end:-1:1);
y2_phase1=example_ts(x2_phase1);
h3=fill([x1_phase1(:);x2_phase1(:)],[y1_phase1(:);y2_phase1(:)],'green','linestyle','none');
hold on
x1_phase1=phase_4;
y1_phase1=zeros(1,length(x1_phase1));
x2_phase1=phase_4(end:-1:1);
y2_phase1=example_ts(x2_phase1);
h4=fill([x1_phase1(:);x2_phase1(:)],[y1_phase1(:);y2_phase1(:)],'blue','linestyle','none');
hold on
plot(t1*ones(1000,1),(linspace(0,12,1000))','k--','linewidth',2);
hold on
plot(t2*ones(1000,1),(linspace(0,12,1000))','k--','linewidth',2);
ylim([0 12]);
xlim([150 290]);
set(gca,'xtick',[datenum(2016,6,15) datenum(2016,7,15) datenum(2016,8,15) datenum(2016,9,15)]-datenum(2016,1,1)+1,...
    'xticklabels',{'Jun15th','Jul15th','Aug15th','Sep15th'});
legend([h1 h2 h3 h4],{'P1','P2','P3','P4'},'fontsize',16,'fontweight','bold');
ylabel('mm');
set(gca,'fontsize',16);

%% plot en la normal

en=[1980 1983 1987 1988 1992 1995 1998 2003 2005 2007 2010 2015 2016]-1;
la=[1984 1985 1989 1996 1999 2000 2001 2006 2008 2009 2011 2012 2017 2018]-1;
year_used=1979:2017;
ne=year_used(~ismember(year_used,en) & ~ismember(year_used,la));

precip_sm=smoothdata(precip,3,'gaussian',31);

msd_en=msd_collection_ncep(ismember(msd_collection_ncep(:,1),en),:);
msd_la=msd_collection_ncep(ismember(msd_collection_ncep(:,1),la),:);
msd_ne=msd_collection_ncep(ismember(msd_collection_ncep(:,1),ne),:);

precip_msd_en=NaN(size(msd_en,1),366);
precip_msd_la=NaN(size(msd_la,1),366);
precip_msd_ne=NaN(size(msd_ne,1),366);

for i=1:size(msd_en,1);
    msd_here=msd_en(i,:);
    precip_here=squeeze(precip_sm(msd_here(2),msd_here(3),(datenum(msd_here(1),1,1):datenum(msd_here(1),12,31))-datenum(1979,1,1)+1));
    precip_should=NaN(366,1);
    
    if length(precip_here)==366
        precip_should=precip_here;
    else
        precip_should([1:59 61:366])=precip_here;
        precip_should(60)=nanmean(precip_should([59 61]));
    end
    
    precip_msd_en(i,:)=precip_should;
end

for i=1:size(msd_la,1);
    msd_here=msd_la(i,:);
    precip_here=squeeze(precip_sm(msd_here(2),msd_here(3),(datenum(msd_here(1),1,1):datenum(msd_here(1),12,31))-datenum(1979,1,1)+1));
    precip_should=NaN(366,1);
    
    if length(precip_here)==366
        precip_should=precip_here;
    else
        precip_should([1:59 61:366])=precip_here;
        precip_should(60)=nanmean(precip_should([59 61]));
    end
    
    precip_msd_la(i,:)=precip_should;
end

for i=1:size(msd_ne,1);
    msd_here=msd_ne(i,:);
    precip_here=squeeze(precip_sm(msd_here(2),msd_here(3),(datenum(msd_here(1),1,1):datenum(msd_here(1),12,31))-datenum(1979,1,1)+1));
    precip_should=NaN(366,1);
    
    if length(precip_here)==366
        precip_should=precip_here;
    else
        precip_should([1:59 61:366])=precip_here;
        precip_should(60)=nanmean(precip_should([59 61]));
    end
    
    precip_msd_ne(i,:)=precip_should;
end

figure
plot(nanmean(precip_msd_en))

figure
plot(nanmean(precip_msd_la))

figure
plot(nanmean(precip_msd_ne))

%% p1 p2 pmin
load('MSD');
MSD=MSD{:,:};
pmin=MSD(:,11);
p2=MSD(:,13);
year_msd=MSD(:,1);
    
en=[1980 1983 1987 1988 1992 1995 1998 2003 2005 2007 2010 2015 2016]-1;
la=[1984 1985 1989 1996 1999 2000 2001 2006 2008 2009 2011 2012 2017 2018]-1;
year_used=1979:2017;
no_en=year_used(~ismember(year_used,en));
no_la=year_used(~ismember(year_used,la));

loc_unique=unique(MSD(:,2:3),'rows');

s_pmin=NaN(120,60,2);
s_p2=NaN(120,60,2);

for i=1:size(loc_unique,1);
    s_pmin(loc_unique(i,1),loc_unique(i,2),1)=nanmean(pmin(MSD(:,2)==loc_unique(i,1) & MSD(:,3)==loc_unique(i,2) & ismember(year_msd,la)));
    s_pmin(loc_unique(i,1),loc_unique(i,2),2)=nanmean(pmin(MSD(:,2)==loc_unique(i,1) & MSD(:,3)==loc_unique(i,2) & ismember(year_msd,no_la)));
    
    s_p2(loc_unique(i,1),loc_unique(i,2),1)=nanmean(p2(MSD(:,2)==loc_unique(i,1) & MSD(:,3)==loc_unique(i,2) & ismember(year_msd,en)));
    s_p2(loc_unique(i,1),loc_unique(i,2),2)=nanmean(p2(MSD(:,2)==loc_unique(i,1) & MSD(:,3)==loc_unique(i,2) & ismember(year_msd,no_en)));
end
    



















    










    
    
    
    
                        
                        
                    
                    
                    
            
            
                    
                    
                    
            
            
            
            
        
        