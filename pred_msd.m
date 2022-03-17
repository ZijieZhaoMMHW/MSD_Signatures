%%

F=detrend(precip_eof,0);
[EOFs,lambda,contribution,PCs]=EOF_analysis(F');
 PCs_sd=PCs(1,:)./std(PCs(1,:));
 
 rainy_season_se=[];
 
 for i=1979:2017
     pc_here=PCs_sd((datenum(i,1,1):datenum(i,12,31))-datenum(1979,1,1)+1);
     r_p=find(pc_here>0);
     period_here=datenum(i,1,1):datenum(i,12,31);
     rainy_season_se=[rainy_season_se;[period_here(r_p(1)) period_here(r_p(end))]];
 end
 
loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

std_explained=NaN(39,987);

for i=339:987
    precip_here=precip_eof(:,i);
    
    for j=1979:2017
        start_here=rainy_season_se(j-1979+1,1)-datenum(1979,1,1)+1;
        end_here=rainy_season_se(j-1979+1,2)-datenum(1979,1,1)+1;
        rainy_here=precip_here(start_here:end_here);
        length_here=length(rainy_here);
        
        equ_here=['a+b*sin(4*pi*x./' num2str(length_here) '+c)'];
        
        model_s_c=fit((start_here:end_here)',rainy_here,equ_here);
        model_fitted=model_s_c((start_here:end_here)');
        
%         plot((start_here:end_here)',rainy_here);
%         hold on
%         plot((start_here:end_here)',model_fitted);
        
        std_explained(j-1979+1,i)=nanstd(model_fitted)./nanstd(rainy_here);
    end
end

msd_existence=NaN(39,987);

for i=1:987;
    loc_here=loc_unique(i,:);
    msd_here=msd_collection_ncep(msd_collection_ncep(:,2)==loc_here(1) & msd_collection_ncep(:,3)==loc_here(2),:);
    
    msd_e=zeros(39,1);
    
    msd_year_here=msd_here(:,1);
    
    msd_e(ismember(1979:2017,msd_year_here))=1;
    
    msd_existence(:,i)=msd_e;
end

%% load ENSO

oni=csvread('ENSO_year_index.csv');
oni=oni((1978:2017)-oni(1,1)+1,:);

oni_ts=oni(:,2:end);
oni_ts=oni_ts(:);

year_used=repmat(1978:2017,1,12);
month_used=repmat(1:12,1,40);

y_m_used=[sort(year_used(:)) month_used(:)];

oni_mean=[];

for i=1979:2017
    index_here=find(y_m_used(:,1)==i & y_m_used(:,2)==2);
    oni_here=oni_ts((index_here-5):1:index_here);
    oni_here=nanmean(oni_here);
    
    oni_mean=[oni_mean;oni_here];
end

plot(std_explained(:,2)./std(std_explained(:,2)))
hold on
plot(oni_mean./std(oni_mean))

see_std=std_explained(:,2)./std(std_explained(:,2));
see_oni=oni_mean./std(oni_mean);

corr([see_oni(2:end) see_std(1:(end-1))])

full_start=[];
full_end=[];

for i=1979:2017
    msd_here=msd_collection_ncep(msd_collection_ncep(:,1)==i,:);
    start_here=nanmean(msd_here(:,5));
    end_here=nanmean(msd_here(:,7));
    
    full_start=[full_start;start_here];
    full_end=[full_end;end_here];
end

csvwrite('std_explained.csv',std_explained);

%% try the full climatology
load('precip_clim');
load('msd_collection_ncep');

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

precip_clim_smooth=smoothdata(precip_clim,3,'gaussian',31);

precip_clim_eof=NaN(366,987);

for i=1:size(loc_unique,1);
    loc_here=loc_unique(i,:);
    precip_here=squeeze(precip_clim_smooth(loc_here(1),loc_here(2),:));
    precip_clim_eof(:,i)=precip_here;
end

addpath('/Applications/MATLAB_R2017a.app/toolbox/EOF_analysis')

F=detrend(precip_clim_eof,0);
[EOFs,lambda,contribution,PCs]=EOF_analysis(F');
 PCs_sd=PCs(1,:)./std(PCs(1,:));
 
 rainy_period=find(PCs_sd>=0);
 start_end=[rainy_period(1) rainy_period(end)];
 
 std_explained_eof=NaN(987,1);
 R_adj=NaN(987,1);
 
 for i=1:987;
     precip_here=precip_clim_eof(:,i);
     rainy_here=precip_here(start_end(1):start_end(2));
     
     length_here=length(rainy_here);
     
     tbl=table((start_end(1):start_end(end))',rainy_here,'variablenames',{'Time','Rain'});
     
     mdl=fitlm(tbl,'Rain~Time+Time^2+Time^3');
     
     model_fitted=mdl.Fitted;
     
%      %equ_here=['a+b*sin(4*pi*x./' num2str(length_here) '+c)'];
%      
%      equ_here=['a+b*x+c*x.^2+d*x.^3+e*x.^4'];
%      
%      model_s_c=fit((start_end(1):start_end(end))',rainy_here,equ_here);
%      model_fitted=model_s_c((start_end(1):start_end(end))');
     
     std_explained_eof(i)=nanstd(model_fitted)./nanstd(rainy_here);
     
     R_adj(i)=mdl.Rsquared.Adjusted;
     
 end
 
 r_plot=NaN(120,60);
 
 predict_plot=NaN(120,60);
 
 for i=1:987;
     loc_here=loc_unique(i,:);
     std_here=std_explained_eof(i);
     r_here=R_adj(i);
     predict_plot(loc_here(1),loc_here(2))=std_here;
     r_plot(loc_here(1),loc_here(2))=r_here;
 end
 
 %% eof_determination_rainy
 load('eof_associated');
EOF_1=EOFs(:,1);
load('msd_collection_ncep');

loc_unique=unique(msd_collection_ncep(:,2:3),'rows');

EOF_1_plot=NaN(120,60);

for i=1:size(loc_unique,1);
    loc_here=loc_unique(i,:);
    EOF_here=EOF_1(i);
    EOF_1(loc_here(1),loc_here(2))=EOF_here;
end

rainy_period=find(PCs_sd>=0);

%% modelling for full region
load('precip_clim');
loc_full_data=[];

precip_clim=smoothdata(precip_clim,3,'gaussian',31);

for i=1:120;
    for j=1:60
        precip_here=squeeze(precip_clim(i,j,:));
        if ~(sum(isnan(precip_here))>0)
            loc_full_data=[loc_full_data;[i j]];
        end
    end
end

precip_2d=NaN(366,2118);

for i=1:size(loc_full_data,1);
    precip_2d(:,i)=squeeze(precip_clim(loc_full_data(i,1),loc_full_data(i,2),:));
end

R_adj=NaN(2118,1);
coef_4=NaN(2118,1);

rainy_period=[138 301];

for i=1:size(precip_2d,2);
     precip_here=precip_2d(:,i);
     rainy_here=precip_here(rainy_period(1):rainy_period(end));
     
     length_here=length(rainy_here);
     
     tbl=table((rainy_period(1):rainy_period(end))',rainy_here,'variablenames',{'Time','Rain'});
     
     mdl=fitlm(tbl,'Rain~Time+Time^2+Time^3+Time^4');
     
     model_fitted=mdl.Fitted;
     
%      %equ_here=['a+b*sin(4*pi*x./' num2str(length_here) '+c)'];
%      
%      equ_here=['a+b*x+c*x.^2+d*x.^3+e*x.^4'];
%      
%      model_s_c=fit((start_end(1):start_end(end))',rainy_here,equ_here);
%      model_fitted=model_s_c((start_end(1):start_end(end))');
     
     R_adj(i)=mdl.Rsquared.Adjusted;
     coef_4(i)=mdl.Coefficients.Estimate(end);
     p_4(i)=mdl.Coefficients.pValue(end);
 end

r_plot_full=NaN(120,60);
coef_plot_full=NaN(120,60);
p_logic_full=NaN(120,60);

for i=1:size(loc_full_data,1);
    r_plot_full(loc_full_data(i,1),loc_full_data(i,2))=R_adj(i);
    coef_plot_full(loc_full_data(i,1),loc_full_data(i,2))=coef_4(i);
    p_logic_full(loc_full_data(i,1),loc_full_data(i,2))=nansum(p_4(i)<=0.05);
end

%% plot for three conditions

% coef>0 sign
load('coef_r_full');
load('p_logic');
load('precip_clim');

loc_positive_sign=[];

for i=1:120;
    for j=1:60;
        p_here=p_logic_full(i,j);
        coef_here=coef_plot_full(i,j);
        
        if p_here==1 && coef_here>0
            loc_positive_sign=[loc_positive_sign;[i j]];
        end
    end
end

precip_p_sign=NaN(366,size(loc_positive_sign,1));

for i=1:size(loc_positive_sign,1);
    precip_p_sign(:,i)=squeeze(precip_clim(loc_positive_sign(i,1),loc_positive_sign(i,2),:));
end

precip_p_sign=nanmean(precip_p_sign,2);

precip_p_sign=smoothdata(precip_p_sign,1,'gaussian',31);




loc_negative_sign=[];

for i=1:120;
    for j=1:60;
        p_here=p_logic_full(i,j);
        coef_here=coef_plot_full(i,j);
        
        if p_here==1 && coef_here<0
            loc_negative_sign=[loc_negative_sign;[i j]];
        end
    end
end

precip_n_sign=NaN(366,size(loc_negative_sign,1));

for i=1:size(loc_negative_sign,1);
    precip_n_sign(:,i)=squeeze(precip_clim(loc_negative_sign(i,1),loc_negative_sign(i,2),:));
end

precip_n_sign=nanmean(precip_n_sign,2);

precip_n_sign=smoothdata(precip_n_sign,1,'gaussian',31);

loc_insign=[];

for i=1:120;
    for j=1:60;
        p_here=p_logic_full(i,j);
        coef_here=coef_plot_full(i,j);
        
        if p_here==0
            loc_insign=[loc_insign;[i j]];
        end
    end
end

precip_insign=NaN(366,size(loc_insign,1));

for i=1:size(loc_insign,1);
    precip_insign(:,i)=squeeze(precip_clim(loc_insign(i,1),loc_insign(i,2),:));
end

precip_insign=nanmean(precip_insign,2);

[x,y]=find(coef_plot_full==nanmin(coef_plot_full(:)));
lon_n=lon_ncep(x);
lat_n=lat_ncep(y);
precip_most_n=smoothdata(squeeze(precip_clim(x,y,:)),1,'gaussian',31);

[x,y]=find(coef_plot_full==nanmax(coef_plot_full(:)));
lon_p=lon_ncep(x);
lat_p=lat_ncep(y);
precip_most_p=smoothdata(squeeze(precip_clim(x,y,:)),1,'gaussian',31);

%% modelling for each period
load('daily_cpc_precip');
loc_full_data=[];

precip_sm=smoothdata(precip,3,'gaussian',31);

for i=1:120;
    for j=1:60
        precip_here=squeeze(precip_sm(i,j,:));
        if ~(sum(isnan(precip_here))>0)
            loc_full_data=[loc_full_data;[i j]];
        end
    end
end

precip_2d=NaN(14245,2115);

for i=1:size(loc_full_data,1);
    precip_2d(:,i)=squeeze(precip_sm(loc_full_data(i,1),loc_full_data(i,2),:));
end

R_adj=NaN(39,2115);
coef_4=NaN(39,2115);

rainy_period=[138 301];

for i=1:size(precip_2d,2);
    tic
    for j=1979:2017
    
     precip_here=squeeze(precip_2d((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1,i));
     rainy_here=precip_here(rainy_period(1):rainy_period(end));
     
     length_here=length(rainy_here);
     
     tbl=table((rainy_period(1):rainy_period(end))',rainy_here,'variablenames',{'Time','Rain'});
     
     mdl=fitlm(tbl,'Rain~Time+Time^2+Time^3+Time^4');
     
     model_fitted=mdl.Fitted;
     
%      %equ_here=['a+b*sin(4*pi*x./' num2str(length_here) '+c)'];
%      
%      equ_here=['a+b*x+c*x.^2+d*x.^3+e*x.^4'];
%      
%      model_s_c=fit((start_end(1):start_end(end))',rainy_here,equ_here);
%      model_fitted=model_s_c((start_end(1):start_end(end))');
     
     R_adj(j-1979+1,i)=mdl.Rsquared.Adjusted;
     coef_4(j-1979+1,i)=mdl.Coefficients.Estimate(end);
     p_4(j-1979+1,i)=mdl.Coefficients.pValue(end);
    end
    toc
end
 
r_mean=nanmean(R_adj);
coef_mean=nanmean(coef_4);

r_plot_full=NaN(120,60);
coef_plot_full=NaN(120,60);

for i=1:size(loc_full_data,1);
    r_plot_full(loc_full_data(i,1),loc_full_data(i,2))=r_mean(i);
    coef_plot_full(loc_full_data(i,1),loc_full_data(i,2))=coef_mean(i);
    %p_logic_full(loc_full_data(i,1),loc_full_data(i,2))=nansum(p_4(i)<=0.05);
end

r_new=NaN(2115,1);

for i=1:size(coef_4,2);
    coef_here=coef_4(:,i);
    oni=tbl_for_cor.Oni;
    mdl=fitlm([(1979:2017)' oni],coef_here);
    r_new(i)=mdl.Rsquared.Adjusted;
end

r_new_plot=NaN(120,60);

for i=1:size(loc_full_data,1);
    r_new_plot(loc_full_data(i,1),loc_full_data(i,2))=r_new(i);
end

% what about mjo and enso
load('daily_cpc_precip');
loc_full_data=[];

precip_sm=smoothdata(precip,3,'gaussian',31);

for i=1:120;
    for j=1:60
        precip_here=squeeze(precip_sm(i,j,:));
        if ~(sum(isnan(precip_here))>0)
            loc_full_data=[loc_full_data;[i j]];
        end
    end
end

precip_2d=NaN(14245,2115);

for i=1:size(loc_full_data,1);
    precip_2d(:,i)=squeeze(precip_sm(loc_full_data(i,1),loc_full_data(i,2),:));
end


addpath('/Volumes/mydrive/msd_data/mjo_data');
date_used=datenum(1974,6,1):datenum(2018,8,30);
rmm1=double(ncread('rmm1.nc','RMM1'));
rmm2=double(ncread('rmm2.nc','RMM2'));
rmm1=rmm1((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);
rmm2=rmm2((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);

oni_index=csvread('ENSO_year_index.csv');

R_adj_mjo=NaN(39,2115);

coef_rmm1=NaN(39,2115);
coef_rmm2=NaN(39,2115);
coef_4=NaN(39,2115);
p_coef=NaN(39,2115);
p_rmm1=NaN(39,2115);
p_rmm2=NaN(39,2115);
rainy_period=[138 301];

for i=1:size(precip_2d,2);
    tic
    for j=1979:2017
    
     precip_here=squeeze(precip_2d((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1,i));
     rmm1_here=squeeze(rmm1((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1));
     rmm2_here=squeeze(rmm2((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1));
     
     
     
     rainy_here=precip_here(rainy_period(1):rainy_period(end));
     
     length_here=length(rainy_here);
     
     tbl=table((rainy_period(1):rainy_period(end))',rainy_here,rmm1_here(rainy_period(1):rainy_period(end)),rmm2_here(rainy_period(1):rainy_period(end)),'variablenames',{'Time','Rain','RMM1','RMM2'});
     
     mdl=fitlm(tbl,'Rain~Time+Time^2+Time^3+Time^4+RMM1+RMM2');
     
     model_fitted=mdl.Fitted;
     
%      %equ_here=['a+b*sin(4*pi*x./' num2str(length_here) '+c)'];
%      
%      equ_here=['a+b*x+c*x.^2+d*x.^3+e*x.^4'];
%      
%      model_s_c=fit((start_end(1):start_end(end))',rainy_here,equ_here);
%      model_fitted=model_s_c((start_end(1):start_end(end))');
     
     R_adj_mjo(j-1979+1,i)=mdl.Rsquared.Adjusted;
     coef_rmm1(j-1979+1,i)=mdl.Coefficients.Estimate(3);
     coef_rmm2(j-1979+1,i)=mdl.Coefficients.Estimate(4);
     coef_4(j-1979+1,i)=mdl.Coefficients.Estimate(end);
     p_coef(j-1979+1,i)=mdl.Coefficients.pValue(end);
     p_rmm1(j-1979+1,i)=mdl.Coefficients.pValue(3);
     p_rmm2(j-1979+1,i)=mdl.Coefficients.pValue(4);
    end
    toc
end

oni_index=csvread('ENSO_year_index.csv');
year_index=oni_index(:,1);
oni_index=oni_index(:,2:end);

R_adj_mjo_enso=NaN(39,2115);

coef_rmm1=NaN(39,2115);
coef_rmm2=NaN(39,2115);
coef_enso=NaN(39,2115);

rainy_period=[138 301];

for i=1:size(precip_2d,2);
    tic
    for j=1979:2017
    
     precip_here=squeeze(precip_2d((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1,i));
     rmm1_here=squeeze(rmm1((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1));
     rmm2_here=squeeze(rmm2((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1));
     
     enso_here=[repmat(oni_index(year_index==j,5),15,1);repmat(oni_index(year_index==j,6),30,1);...
         repmat(oni_index(year_index==j,7),31,1);repmat(oni_index(year_index==j,8),31,1);...
         repmat(oni_index(year_index==j,9),30,1);repmat(oni_index(year_index==j,10),27,1)];
     
     
     rainy_here=precip_here(rainy_period(1):rainy_period(end));
     
     length_here=length(rainy_here);
     
     tbl=table((rainy_period(1):rainy_period(end))',rainy_here,rmm1_here(rainy_period(1):rainy_period(end)),rmm2_here(rainy_period(1):rainy_period(end)),enso_here,'variablenames',{'Time','Rain','RMM1','RMM2','ENSO'});
     
     mdl=fitlm(tbl,'Rain~Time+Time^2+Time^3+Time^4+RMM1+RMM2+ENSO');
     
     model_fitted=mdl.Fitted;
     
%      %equ_here=['a+b*sin(4*pi*x./' num2str(length_here) '+c)'];
%      
%      equ_here=['a+b*x+c*x.^2+d*x.^3+e*x.^4'];
%      
%      model_s_c=fit((start_end(1):start_end(end))',rainy_here,equ_here);
%      model_fitted=model_s_c((start_end(1):start_end(end))');
     
     R_adj_mjo_enso(j-1979+1,i)=mdl.Rsquared.Adjusted;
     coef_rmm1(j-1979+1,i)=mdl.Coefficients.Estimate(end-2);
     coef_rmm2(j-1979+1,i)=mdl.Coefficients.Estimate(end-1);
     coef_enso(j-1979+1,i)=mdl.Coefficients.Estimate(end);
    end
    toc
end

r_mean=nanmean(R_adj_mjo);

coef_rmm1_mean=nanmean(coef_rmm1);

coef_rmm2_mean=nanmean(coef_rmm2);

r_plot_full=NaN(120,60);

coef_rmm1_mean_plot=NaN(120,60);

coef_rmm2_mean_plot=NaN(120,60);

for i=1:size(loc_full_data,1);
    r_plot_full(loc_full_data(i,1),loc_full_data(i,2))=r_mean(i);
    coef_rmm1_mean_plot(loc_full_data(i,1),loc_full_data(i,2))=coef_rmm1_mean(i);
    coef_rmm2_mean_plot(loc_full_data(i,1),loc_full_data(i,2))=coef_rmm2_mean(i);
    %p_logic_full(loc_full_data(i,1),loc_full_data(i,2))=nansum(p_4(i)<=0.05);
end

figure
m_contourf(lon,lat,(coef_rmm1_mean_plot)',-2:0.01:2,'linestyle','none');
m_coast();
m_grid;
colormap(color_here);
caxis([-1.5 1.5]);

figure
m_contourf(lon,lat,(coef_rmm2_mean_plot)',-2:0.01:2,'linestyle','none');
m_coast();
m_grid;
colormap(color_here);
caxis([-1.5 1.5]);

%% climatological mjo

addpath('/Volumes/mydrive/msd_data/mjo_data');
date_used=datenum(1974,6,1):datenum(2018,8,30);
rmm1=double(ncread('rmm1.nc','RMM1'));
rmm2=double(ncread('rmm2.nc','RMM2'));
rmm1=rmm1((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);
rmm2=rmm2((datenum(1979,1,1):datenum(2017,12,31))-datenum(1974,6,1)+1);
    
rainy_period=[138 301];
rmm1_full=NaN(39,164);
rmm2_full=NaN(39,164);

for j=1979:2017
    rmm1_here=squeeze(rmm1((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1));
    rmm2_here=squeeze(rmm2((datenum(j,1,1):datenum(j,12,31))-datenum(1979,1,1)+1));
    
    rmm1_full(j-1979+1,:)=rmm1_here(rainy_period(1):rainy_period(2));
    rmm2_full(j-1979+1,:)=rmm2_here(rainy_period(1):rainy_period(2));
end

rmm1_clim=nanmean(rmm1_full);
rmm2_clim=nanmean(rmm2_full);

plot(rmm1_clim,rmm2_clim,'.-');

figure('pos',[10 10 1500 1500]);
m_proj('miller','lon',[180+60 180+120],'lat',[0 30]);
h=tight_subplot(3,3,[0.025,0.025],[0.035 0.035],[0.03,0.05]);


    
    




 

    


        
    

        
        
     





