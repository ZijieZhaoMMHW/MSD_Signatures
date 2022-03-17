figure('pos',[10 10 1500 1500]);
subplot(1,2,1);
lin1x=linspace(0,0,1000);lin1y=linspace(1,4,1000);
lin2x=linspace(sqrt(2)/2,4,1000);lin2y=linspace(sqrt(2)/2,4,1000);
lin3x=linspace(1,4,1000);lin3y=linspace(0,0,1000);
lin4x=linspace(sqrt(2)/2,4,1000);lin4y=linspace(-4,-sqrt(2)/2,1000);
lin5x=linspace(0,0,1000);lin5y=linspace(-4,-1,1000);
lin6x=linspace(-4,-sqrt(2)/2,1000);lin6y=linspace(-4,-sqrt(2)/2,1000);
lin7x=linspace(-4,-1,1000);lin7y=linspace(0,0,1000);
lin8x=linspace(-4,-sqrt(2)/2,1000);lin8y=linspace(sqrt(2)/2,4,1000);

t=0:pi/100:2*pi;
x=sin(t);
y=cos(t);

plot(x,y,'b','linewidth',2);
hold on
plot([lin1x(:) lin2x(:) lin3x(:) lin4x(:) lin5x(:) lin6x(:) lin7x(:) flipud(lin8x(:))],...
    [lin1y(:) lin2y(:) lin3y(:) flipud(lin4y(:)) lin5y(:) lin6y(:) lin7y(:) lin8y(:)],'b','linewidth',2);

xlabel('RMM1','fontweight','bold','fontsize',16);
ylabel('RMM2','fontweight','bold','fontsize',16);

title('Phase defined by RMM1/RMM2','fontsize',16,'fontweight','bold');
set(gca,'fontsize',16);
set(gca,'fontsize',16);

subplot(1,2,2);
load('precip_clim');
example_ts=smoothdata(squeeze(precip_clim(45,18,:)),1,'gaussian',31);
period_1=150:200;
period_2=250:300;

t1=period_1(find(example_ts(period_1)==nanmax(example_ts(period_1))));
t2=period_2(find(example_ts(period_2)==nanmax(example_ts(period_2))));

period_p=t1:t2;
tmin=period_p(find(example_ts(t1:t2)==nanmin(example_ts(t1:t2))));
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
