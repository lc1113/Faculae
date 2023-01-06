% patch analysis for NIR-GECO. this script is for finding AP peaks from DAQ
% data and calculate the Ca sensitivity of GECI.
% reference: stimulated AP multicycles
clear all; clc;
%% load and select
% Movie loading path
dir_process = 'D:\lc_data\20221018_LIN_CAV24\A7\cell2\';
dire = -1;  % dire = 1 means positive GECI; dire = -1 means negative GECI
peak_thres = 20;    % mV

pathname = uigetdir(dir_process,'Select patch-imaging subfolder');
listing = dir(pathname);
for i =1:size(listing,1)
    if listing(i).isdir == 1 && ~strcmp(listing(i).name,'.')  && ~strcmp(listing(i).name,'..') && exist([pathname '\' listing(i).name '\matlab variables.mat'])
        load ([pathname '\' listing(i).name '\matlab variables.mat']);
        break
    end
end
if exist([dir_process '\patch param.txt'])
    b = importdata([dir_process '\patch param.txt']);
    datab = b.data;Cm = datab(1,5);Rm = datab(1,4);Ra = datab(1,3);%pF 

else
    Rm = 'N.A.';Cm = 'N.A.';Ra = 'N.A.';
end
% constants
% cycles = 80;

if exist([pathname '\movie_info.txt'])
    c = importfile([pathname '\movie_info.txt']);
elseif exist([pathname '\movie.txt'])
    c = importfile([pathname '\movie.txt']);
end
movname = '\movie.bin';
ncol = c.DO(find(strncmp(c.laser,'nrow',4) == 1));         % x invert
nrow = c.DO(find(strncmp(c.laser,'ncol',4) == 1));         % y invert
camera_bias = c.DO(find(strncmp(c.laser,'Binning',7) == 1)).^2*100;  % background due to camera bias (100 for bin 1x1)
dt_mov = c.DO(find(strncmp(c.laser,'Exposure',8) == 1));    % exposure time in millisecond (484 Hz)
Fs = samprate;
DAQname = '\movie_DAQ.txt';
dnsamp = Fs/(1000/dt_mov);        % downsampling rate = DAQ rate/camera rate
dnsamp = round(dnsamp);

% load DAQ data
% load DAQ data
tmp = importdata([pathname DAQname]);   % import data
data = tmp.data;                    % get array
Vm = data(:,2)*100;                 % Vm in millivolt, column vector
dt_daq = dt_mov/dnsamp;             % DAQ dt in millisecond
t_daq = [0:length(Vm)-1]*dt_daq/10^3;       % DAQ time axis in second
a=importdata([pathname '\movie_DAQ.txt']);
data=a.data;
AI_scaled=data(:,1);
AI_10Vm=data(:,2)*100;
time=(1:length(AI_scaled)')'./Fs;
figure;
set(gcf,'outerposition',get(0,'screensize'));
plot(time,AI_scaled,time,AI_10Vm);
legend('AI\_scaled','Vm (mV)','Location','Northeast');
hold on
xlim=[0,max(time)];
ylim=[-70,-60];
xL=xlim;yL=ylim;
set(gca,'xtick',[0:5:max(time)])
% plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
box off
axis([xL yL])
axis tight
saveas(gca,[pathname '\0 waveform of AI.fig']);
saveas(gca,[pathname '\0 waveform of AI.png']);

%% loading the video movie
% load movie
fname = [pathname movname];
[mov, nframe] = readBinMov(fname, nrow, ncol);
mov = single(mov);img = mean(mov, 3);
% [~,intens_raw] = clicky_v2(mov, dt_mov, c.DO(find(strncmp(c.laser,'Binning',7) == 1)), img,camera_bias,max(img,[],'all'), 'Please select interested regions');
[~, intens_raw] = clicky(mov, img, 'select only 1 ROI, right click when done');
intens_rembkg = intens_raw(:,1)-intens_raw(:,2);
% select ROI for analysis
background = mean(intens_raw(:,size(intens_raw,2)));

saveas(gca,[pathname '\1 clicky analysis.fig']);
saveas(gca,[pathname '\1 clicky analysis.png']);
len = size(intens_raw,1);
t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second

%% dump kernel by left click on axes
h = figure('Name','dump kernel by left click on axes');
set(gcf,'Position',get(0,'ScreenSize'));
ax = {};dumpIndex = ones(1,cycles);
headPts_fluo = size(headPts,2)/dnsamp;    tailPts_fluo = size(tailPts,2)/dnsamp;
bleachPts =  2*samprate/dnsamp;    cyclePts = (hiPts+lowPts)/dnsamp;
for i = 1:cycles
    kernel_Vm(1:hiPts+lowPts,i) = AI_10Vm(size(headPts,2)+(i-1)*(hiPts+lowPts)-round(2*samprate)+1:size(headPts,2)+i*(hiPts+lowPts)-round(2*samprate));
    kernel_Fluo_dump(1:cyclePts,i) = intens_raw(headPts_fluo+(i-1)*cyclePts-round(bleachPts)+1:headPts_fluo+i*cyclePts-round(bleachPts));
    ax{2*i-1} = subplot(cycles,2,2*i-1);    plot([0:hiPts+lowPts-1],kernel_Vm(:,i));
    ax{2*i} = subplot(cycles,2,2*i);    plot([0:cyclePts-1],kernel_Fluo_dump(:,i));
    if i<cycles
        set(ax{2*i-1},'xticklabel',[]);    set(ax{2*i},'xticklabel',[]);
    end
    
end
hold on
% set([h gcf],'hittest','off')           % turn off hittest 
set(h.Children,'buttondownfcn',{@buttondownfcn,cycles,ax,pathname});         % assign function to gca
saveas(h,[pathname '\2 dumped kernel.fig']);
saveas(h,[pathname '\2 dumped kernel.png']);
%% photobleaching correction
% exponential fitting on remained kernels
intens_rembkg_norm = intens_rembkg./mean(intens_rembkg(1:round(headPts_fluo/100)));
plot(t_mov,intens_rembkg_norm');
F = @(x,xdata) x(1).*exp(-x(2).*xdata) + x(3).*exp(-x(4).*xdata);
x0 = [0.9 0.004 0.1 0.0001];
period_dataPts = [ones(1,round(bleachPts)), zeros(1,cyclePts-round(bleachPts))];
period_dataPts_false = [zeros(1,round(bleachPts)), zeros(1,cyclePts-round(bleachPts))];
% Cycle each period to give steps
dataPts = [ones(1,round(headPts_fluo)-round(bleachPts))];
for i = 1:cycles
if logical(dumpIndex(i))
    dataPts = [dataPts period_dataPts];
else
    dataPts = [dataPts period_dataPts_false];
end
end
dataPts = [dataPts zeros(1,round(bleachPts)) ones(1,round(tailPts_fluo))];
xdata = t_mov(logical(dataPts)); ydata = double(intens_rembkg_norm(logical(dataPts))');
x = lsqcurvefit(F,x0,xdata,ydata)
hold on
plot(t_mov,F(x,t_mov));
hold off
saveas(gca,[pathname '\3_1 photobleaching correction.fig']);
saveas(gca,[pathname '\3_1 photobleaching correction.png']);
close(gcf);
intens_corr = intens_rembkg_norm./F(x,t_mov)';
plot(t_mov,intens_rembkg_norm./F(x,t_mov)');
saveas(gca,[pathname '\3_2 trace after photobleaching correction.fig']);
saveas(gca,[pathname '\3_2 trace after photobleaching correction.png']);
%% calculate Ca deltaF/F versus AP number
% generate kernel
figure()
tkernel_Vm = repmat([0:hiPts+lowPts-1]*dt_daq/10^3,cycles,1)';
tkernel_Fluo = repmat([0:cyclePts-1]*dt_mov/10^3,cycles,1)';
locs_Vpeak = {};
Num_Vpeak = [];
for i = 1:cycles
kernel_Fluo(1:cyclePts,i) = intens_corr(headPts_fluo+(i-1)*cyclePts-round(bleachPts)+1:headPts_fluo+i*cyclePts-round(bleachPts));
[~,locs] = findpeaks(kernel_Vm(:,i),'MinPeakHeight',peak_thres);
Num_Vpeak=[Num_Vpeak size(locs,1)];
locs_Vpeak{i} = locs;

end
kernel_Fluo_smo = smoothdata(kernel_Fluo,1,'sgolay',30);
if dire == 1
    FluoPeak = max(kernel_Fluo_smo,[],1);
elseif dire == -1
    FluoPeak = min(kernel_Fluo_smo,[],1);
end
FluoSteady = mean(kernel_Fluo_smo(1:round(bleachPts),:),1);
SensiPeak = (FluoPeak-FluoSteady)./FluoSteady;
FluoPeak_half = (FluoPeak+FluoSteady)./2;
kernel_Fluo_norm = kernel_Fluo./FluoSteady;
kernel_Fluo_smo_norm = kernel_Fluo_smo./FluoSteady;
SNR = abs(FluoPeak-FluoSteady)./std(kernel_Fluo(round(bleachPts)-round(bleachPts*0.25):round(bleachPts),:),1);
% h = figure('Name','dump kernel by left click on axes');
% set(gcf,'Position',get(0,'ScreenSize'));
thalf_rise = [];
thalf_decay = [];
for i = 1:cycles
    if dumpIndex(i) == 1
        subplot(cycles,1,i);
        plot(tkernel_Fluo(:,i),kernel_Fluo_norm(:,i),tkernel_Fluo(:,i),kernel_Fluo_smo_norm(:,i));
        hold on
        Peak_Pt = find(kernel_Fluo_smo(:,i)==FluoPeak(i));
        if dire == -1
            thalf_rise_Pt = min(find(kernel_Fluo_smo(round(bleachPts)+1:Peak_Pt(1),i) <= FluoPeak_half(i)));
            thalf_decay_Pt = min(find(kernel_Fluo_smo(Peak_Pt(1):end,i) >= FluoPeak_half(i)));
        else
            thalf_rise_Pt = min(find(kernel_Fluo_smo(round(bleachPts)+1:Peak_Pt(1),i) >= FluoPeak_half(i)));
            thalf_decay_Pt = min(find(kernel_Fluo_smo(Peak_Pt(1):end,i) <= FluoPeak_half(i)));
        end
        thalf_rise = [thalf_rise thalf_rise_Pt*dt_mov];
        thalf_decay = [thalf_decay thalf_decay_Pt*dt_mov];
        plot(tkernel_Fluo(round(bleachPts)+thalf_rise_Pt,i),kernel_Fluo_smo_norm(round(bleachPts)+thalf_rise_Pt,i),'r.');
        plot(tkernel_Fluo(Peak_Pt(1)+thalf_decay_Pt-1,i),kernel_Fluo_smo_norm(Peak_Pt(1)+thalf_decay_Pt,i),'r.');
        hold off
    else
        thalf_rise = [thalf_rise NaN];
        thalf_decay = [thalf_decay NaN];
    end
end
saveas(gca,[pathname '\4 norm kernel with thalf.fig']);
saveas(gca,[pathname '\4 norm kernel with thalf.png']);
figure()
plot(tkernel_Fluo(:,logical(dumpIndex)),kernel_Fluo_norm(:,logical(dumpIndex)));
saveas(gca,[pathname '\5-1 norm kernel stack.fig']);
saveas(gca,[pathname '\5-1 norm kernel stack.png']);
close(gcf)
figure()
plot(tkernel_Fluo(:,logical(dumpIndex)),kernel_Fluo_smo_norm(:,logical(dumpIndex)));
saveas(gca,[pathname '\5-1 norm smooth kernel stack.fig']);
saveas(gca,[pathname '\5-1 norm smooth kernel stack.png']);
close(gcf)
%%
Num_Vpeak_save = Num_Vpeak(logical(dumpIndex));
SensiPeak_save = SensiPeak(logical(dumpIndex));
thalf_rise_save = thalf_rise(logical(dumpIndex));
thalf_decay_save = thalf_decay(logical(dumpIndex));
SNR_save = SNR(logical(dumpIndex));

save([pathname '\analysis.mat'],'AI_10Vm','dumpIndex','FluoPeak','FluoSteady','intens_rembkg','intens_corr','kernel_Fluo','kernel_Fluo_norm','kernel_Fluo_smo','kernel_Fluo_smo_norm','kernel_Vm','tkernel_Fluo','tkernel_Vm');
xlswrite([pathname '\analysis.xlsx'],{'AP No.','Peak response','t_half_rise','t_half_decay','SNR','Ra','Rm','Cm'},'Raw','A1');
xlswrite([pathname '\analysis.xlsx'],Num_Vpeak_save','Raw','A2');
xlswrite([pathname '\analysis.xlsx'],SensiPeak_save','Raw','B2');
xlswrite([pathname '\analysis.xlsx'],thalf_rise_save','Raw','C2');
xlswrite([pathname '\analysis.xlsx'],thalf_decay_save','Raw','D2');
xlswrite([pathname '\analysis.xlsx'],SNR_save','Raw','E2');
xlswrite([pathname '\analysis.xlsx'],Ra,'Raw','F2');
xlswrite([pathname '\analysis.xlsx'],Rm,'Raw','G2');
xlswrite([pathname '\analysis.xlsx'],Cm,'Raw','H2');

Num_Vpeak_uni = unique(Num_Vpeak(logical(dumpIndex)));
for i = 1:size(Num_Vpeak_uni,2)
    lo = (Num_Vpeak_save == Num_Vpeak_uni(i));
    SensiPeak_uni(i) = mean(SensiPeak_save(lo));
    thalf_rise_uni(i) = mean(thalf_rise_save(lo));
    thalf_decay_uni(i) = mean(thalf_decay_save(lo));
    SNR_uni(i) = mean(SNR_save(lo));
end
xlswrite([pathname '\analysis.xlsx'],{'AP No.','Peak response','t_half_rise','t_half_decay','SNR'},'Average','A1');
xlswrite([pathname '\analysis.xlsx'],Num_Vpeak_uni','Average','A2');
xlswrite([pathname '\analysis.xlsx'],SensiPeak_uni','Average','B2');
xlswrite([pathname '\analysis.xlsx'],thalf_rise_uni','Average','C2');
xlswrite([pathname '\analysis.xlsx'],thalf_decay_uni','Average','D2');
xlswrite([pathname '\analysis.xlsx'],SNR_uni','Average','E2');
%% functions
function buttondownfcn(hobj,~,cycles,ax,pathname)
dumpIndex = evalin('base','dumpIndex');
for i = 1:cycles
    if isequal(ax{2*i-1}.Position,get(hobj,'Position')) | isequal(ax{2*i}.Position,get(hobj,'Position'))
        dumpIndex(i) = 0;
        title(ax{2*i-1},'dumped kernel');
        title(ax{2*i},'dumped kernel');
    end
end
% set(hobj.Parent,'userdata',dumpIndex);      % add this line insdie the function
disp(dumpIndex);
assignin('base','dumpIndex',dumpIndex);
saveas(hobj.Parent,[pathname '\2 dumped kernel.fig']);
saveas(hobj.Parent,[pathname '\2 dumped kernel.png']);
end