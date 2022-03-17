%% mjo annual model

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

%% annual model

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


R_adj_mjo=NaN(39,2115);

coef_4=NaN(39,2115);
p_coef=NaN(39,2115);
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
     
     R_adj_mjo(j-1979+1,i)=mdl.Rsquared.Adjusted;
     coef_4(j-1979+1,i)=mdl.Coefficients.Estimate(end);
     p_coef(j-1979+1,i)=mdl.Coefficients.pValue(end);
    end
    toc
end

%% climatology
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
p_4=NaN(2118,1);
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