% this program is for library screens with stardist/cellpose process(dual channels), a typical experiment is
% selecting NIR FP mutants with higher brightness in 293T cells
% (LLP cell line express ,miRFP680-3xNLS-P2A-PAmCherry-NLS-EGFP-3xNLS).
% For screens that need two channels, which
% including one marker channel(SDC-GFP in this case) and one library
% channel (SDC-Cy5 in this case), please set the marker channel as the
% first channel in MDA list while acquiring, otherwise this program will
% use wrong channel for calculating.
clear all;clc;
% [WARNING] the varient 'background' should be assigned as zero if flat field
% correction is proceed.

% This script is written by Chang Lin
%% MM FIJI Startup
simulate = 0;   % 0 for screening in A317, 1 for simulation on PC.
EdgeArea = 40;   % minimum area size (um2) 200 or 300 for cell membrane, 50 for nuclei
EdgeDis = 8;    % minimum edge distance (um) 10 for cell membrane, 8 for nuclei  
import org.micromanager.internal.MMStudio;
import mmcorej.*;
import org.micromanager.api.*;
import ij.*;
gui = MMStudio(false);  % start up MM2.0 gui within MATLAB
mmc = gui.getCore;  % get the MM2.0 CMMCore
acq = gui.getAcquisitionEngine; % org.micromanager.acquisition.internal.AcquisitionWrapperEngine
MDA = gui.getAcquisitionManager();
IJ = ij.IJ;
slm = gui.live();   % org.micromanager.internal.SnapLiveManager
Miji;   % start up FIJI within MATLAB
cd('C:\Users\ZouOptics\Desktop\MM2.0\code\matlab');
mmc.setConfig("System","Startup");
mmc.setConfig("System","Startup");  % this line repeat is necessary for 561/594 startup
%% path selection
dir_flat = "C:\Users\ZouOptics\Desktop\MM2.0\flatfield\20220407";    % get direction of flat field on Z1
if simulate == 1
    dir_flat = "E:\pku\research_group\topic_screening_platform\flatfield\20220407";    % get direction of flat field on PC
end
Channel_userdif = str2num(cell2mat(inputdlg({'Channel Num'},'input channel number',[1 35])))
switch(Channel_userdif)
    case{2}
        [file_flatMarker,path_flat] = uigetfile(dir_flat + "\*.*",'select bg-substrated Marker Channel flat field image');
        file_flatMarker = string(file_flatMarker);
    case{1}
        file_flatMarker = [];
end
[file_flatLib,path_flat] = uigetfile(dir_flat + "\*.*",'select bg-substrated Library Channel flat field image');
file_flatLib = string(file_flatLib);
file_flat = [file_flatMarker file_flatLib];
[file_background1,path_flat] = uigetfile(dir_flat + "\*.*",'select marker channel black image');
[file_background2,path_flat] = uigetfile(dir_flat + "\*.*",'select library channel black image');
file_background = [string(file_background1) string(file_background2)];
% [file_ilastikModel,path_ilastikModel] = uigetfile(dir_ilastikModel + "\*.*",'select ilastik model');
path_flat = replace(path_flat,'\','/'); % for FIJI reading

%% spot ROI test and setting
% after acquisition, set up optical path and find the spot location,
% select the spot and run this section
IJ.runMacro("roiManager('reset');");
IJ.runMacro('roiManager("Add");');
MIJ.run("Set Measurements...", "area mean centroid stack redirect=None decimal=3");
IJ.runMacro("roiManager('Associate', 'true');");
IJ.runMacro('roiManager("Measure");');
spotROI = MIJ.getResultsTable();
IJ.runMacro("IJ.deleteRows(0,1);");
Xspot = spotROI(1,3);
Yspot = spotROI(1,4);
%% after configure, copy the dir_process to below
dir_process = "D:\data\20220811_lin_attb06\attb06-2"; % get direction to process
dir = replace(dir_process,'\','/');
%% saving spot ROI
IJ.runMacro("roiManager('Select', 0);");
IJ.runMacro("roiManager('rename', 'laser spot');");
IJ.runMacro("roiManager('save selected', '" + dir +"/laser spot.zip"+ "');");
IJ.runMacro("roiManager('Deselect');");
IJ.runMacro("roiManager('Delete');");
clear spotROI;
%% please create grid and start the screening
% for cellpose screening, use 28x28 grid
% for stardist screening, use 30x30 grid
% keep the MDA window open!!!
path_raw = dir_process + "\rawdata";
mkdir(path_raw);
Channel_snap = strings(Channel_userdif,1);
for i = 1:Channel_userdif
list = {'SDC-LED','SDC-GFP','SDC-mCherry','SDC-Cy5'};
[indx,tf] = listdlg('PromptString',['Select channel_' num2str(i) ':'],...
                           'SelectionMode','single',...
                           'ListString',list);
Channel_snap(i) = list(indx);
end
pl_snap = gui.getPositionList(); % class org.micromanager.PositionList
Num_posi = pl_snap.getNumberOfPositions();
Num_grid = sqrt(Num_posi);
msp_snap = pl_snap.getPositions();
position_slice = zeros(Num_posi,3);
position_label = string(zeros(Num_posi,1));
XYStage_Label = msp_snap(1).getDefaultXYStage();
PFSStage_Label = msp_snap(1).getDefaultZStage();
PhysicalSize = double(mmc.getPixelSizeUm());
pl_snap.save(strcat(dir,'/screens.pos'));
for i = 1:Num_posi
    position_slice(i,1) = msp_snap(i).getX;
    position_slice(i,2) = msp_snap(i).getY;
    position_slice(i,3) = msp_snap(i).getZ;
    position_label(i,1) = msp_snap(i).getLabel();
end
for i = 1:Num_grid:Num_posi
    pl_snap.setPositions(msp_snap(i:i+Num_grid-1));
    for j = 1:Channel_userdif
        MDA.loadAcquisition('D:\Softwares\MM2.0\acqsetting\AcqSettings_'+Channel_snap(j)+'_multiP.txt');
        MDA.runAcquisition(Channel_snap(j),path_raw);
        IJ.runMacro('close();');
    end
end
%% process multi-channel single-time-point ome-xml tiff
[file_ome,path_ome,indx] = uigetfile(path_raw + "\*.*",'select marker channel image');
for j =1:Channel_userdif
    for i = 1:Num_grid
        file_omeopen = strcat(Channel_snap(j),'_',num2str(i),file_ome(end-28:end-11),num2str(i-1,"%03d"),file_ome(end-7:end));
        path_omeopen = strcat(path_raw,'\',Channel_snap(j),'_',num2str(i),'\');
        MIJ.run("Bio-Formats", "open=[" + strcat(path_omeopen,file_omeopen) + "] color_mode=Default concatenate_series open_all_series split_channels view=Hyperstack stack_order=XYCZT"); % open "before" ome-tiff file
        IJ.runMacro("rename('"+Channel_snap(j)+"_"+num2str(i)+"');");
    end
    Conc = 'image'+string([1:Num_grid]')+'=['+Channel_snap(j)+'_'+string([1:Num_grid]')+']';
    MIJ.run("Concatenate...",strcat("title=",Channel_snap(j)," ",join(Conc,1)));
end
serialNum = Num_posi;
reader=bfGetReader([path_ome,file_ome]);
omeMeta = reader.getMetadataStore();
PixelsizeC =  omeMeta.getPixelsSizeC(0);
PixelsizeX = omeMeta.getPixelsSizeX(0);
PixelsizeY = omeMeta.getPixelsSizeY(0);
PhysicalSize = double(omeMeta.getPixelsPhysicalSizeX(0).value());
sizeC = double(PixelsizeC.getNumberValue());
sizeX = double(PixelsizeX.getNumberValue());
sizeY = double(PixelsizeY.getNumberValue());
clear indx PixelsizeC PixelsizeX PixelsizeY;

% flat field correction and generation of image sequence for each channel
Mean_flat = zeros(1,size(file_flat,2));
MIJ.run("Set Measurements...", "area mean centroid stack redirect=None decimal=3");
for i = 1:size(file_flat,2)
    IJ.runMacro("open('"+path_flat+file_flat(i)+"');");
    MIJ.run("Select All");
    MIJ.run("Measure");
    flat_data = MIJ.getResultsTable();
    Mean_flat(i) = flat_data(2);
    IJ.runMacro("IJ.deleteRows(0,1);");
    IJ.runMacro("open('"+path_flat+file_background(i)+"');");
end
window_selected = MIJ.getListImages;
window_selected = string(window_selected);
window_selected(find(strncmp(window_selected,'Preview',7) == 1))= [];
% Channel = strings(sizeC,1);
% for no= 0:sizeC-1 
%     Channel(no+1,1) = omeMeta.getChannelName(0,no);
% end
Channel = Channel_snap;
% temp_Channel = Channel_snap(:,ones(1,2))';
% Channel = reshape(temp_Channel',numel(temp_Channel),1);
% Channel = strcat(Channel,["-before" "-before" "-after" "-after"]');
% process "before" and "after" image
for no= 1:size(Channel,1)
    MIJ.run("Calculator Plus", "i1=["+window_selected(no)+"] i2=["+file_background(2-mod(no,2))+"] operation=[Subtract: i2 = (i1-i2) x k1 + k2] k1=1 k2=0 create");
    MIJ.selectWindow(window_selected(no));
    IJ.runMacro('close();');
    MIJ.selectWindow("Result");
    MIJ.run("Properties...", "channels=1 slices=1 frames="+num2str(serialNum)+" unit=micron pixel_width="+num2str(PhysicalSize)+" pixel_height="+num2str(PhysicalSize)+" voxel_depth=1.0000000 frame=[0.00 sec]");
    IJ.runMacro("rename('"+Channel(no)+"');");
    MIJ.run("Calculator Plus", "i1=["+Channel(no)+"] i2=["+file_flat(2-mod(no,2))+"] operation=[Divide: i2 = (i1/i2) x k1 + k2] k1="+Mean_flat(2-mod(no,2))+" k2=0 create");
    MIJ.selectWindow(Channel(no));
    IJ.runMacro('close();');
    MIJ.selectWindow("Result");
    MIJ.run("Properties...", "channels=1 slices=1 frames="+num2str(serialNum)+" unit=micron pixel_width="+num2str(PhysicalSize)+" pixel_height="+num2str(PhysicalSize)+" voxel_depth=1.0000000 frame=[0.00 sec]");
    IJ.runMacro("rename('"+Channel(no)+"');");
    path_save = dir + "/channel_" + Channel(no,1);
    mkdir(path_save);
    MIJ.run("Image Sequence... ", "format=TIFF name=stack_ save=[" + path_save + "]");
end;
for i = 1:size(file_flat,2)
    MIJ.selectWindow(file_flat(i));
    IJ.runMacro('close();');
    MIJ.selectWindow(file_background(i));
    IJ.runMacro('close();');
end
%% run stardist in a new FIJI window (for stardist screening)
% import post marker channel bg-substracted image sequence by bioformat, be
% sure that the t axis is set as the third axis. run stardist using default parameters and save the ROIs. 
% [[NOTE]] you can reduce "persentile high" to 99 to eliminate the overexpressed cells

%% run cellpose in anaconda powershell prompt (for cellpose screening)
% conda activate cellpose
[file_cellpose,path_cellpose,indx] = uigetfile(dir_process + "\*.*",'select processed post-marker images for segmentation');
system(['C:\Users\ZouOptics\anaconda3\envs\cellpose\python -m cellpose'...
    ' --dir ' [path_cellpose]...
    ' --pretrained_model cyto --chan 0 --use_gpu --fast_mode --save_png --no_npy']);
% open the processed image in a new FIJI by bioformat, then use
% "Plugins-BIOP-Image Analysis-ROIs-Label Image to ROIs" to get the ROI.
% saving the ROIs as .zip
%  --diameter 30
%% process stardist/cellpose returned result for confocal
[file_pro,path_pro,indx_pro] = uigetfile(dir_process + "\*.*",'select processed RoiSet');
MIJ.run("ROI Manager...");
IJ.runMacro("roiManager('reset');");
IJ.runMacro("roiManager('Open','"+replace(path_pro,'\','/') + file_pro+"');");
ResultsTable = cell(size(Channel,1),1);
IJ.runMacro('roiManager("count");');
MIJ.selectWindow("Log");
roiManagerSize = str2num(MIJ.getLog());
MIJ.run("Close" );
% IJ.runMacro("roiManager('save', '" + dir +"/libraryROI.zip"+ "');");
for no= 0:size(Channel,1)-1
    MIJ.selectWindow(Channel(no+1));
    %IJ.runMacro('roiManager("translate", 1, 0);');
    IJ.runMacro('roiManager("Measure");');
    ResultsTable(no+1,1) = {MIJ.getResultsTable()};
    IJ.runMacro("IJ.deleteRows(0, "+roiManagerSize+");");
end
MeanIntMarker = ResultsTable{1,1}(:,2);
MeanIntLib = ResultsTable{2,1}(:,2);
Area = ResultsTable{1,1}(:,1);
Slice = ResultsTable{1,1}(:,5);
Xrev = ResultsTable{1,1}(:,3);
Yrev = ResultsTable{1,1}(:,4);
Xabs = zeros(roiManagerSize,1);
Yabs = zeros(roiManagerSize,1);
Zabs = zeros(roiManagerSize,1);
camera = string(omeMeta.getDetectorID(0,0));
if camera == "camera_confocal"
    theta = atan(16/(2454-512));    %camera confocal# num calculated from file G:\lc data\20210227_LC_dish\dish test 2\CELL_8\\CELL_7_MMStack_1-Pos000_000.ome.tif
else
    theta = 0;  %camera widefield unknown
end
Xcor = Xrev-Xspot;
Ycor = Yrev-Yspot;
Xcor_rev = Xcor*cos(theta)+Ycor*sin(theta);
Ycor_rev = Ycor*cos(theta)-Xcor*sin(theta);
Xabs = position_slice(Slice,1)+Xcor_rev;
Yabs = position_slice(Slice,2)+Ycor_rev;
Zabs = position_slice(Slice,3);
PAstate = ones(roiManagerSize,1); % 1 refers to unlabeled, 2 refers to labeled, 3 refers to discard
Index = [1:roiManagerSize]';
AnalysisTable = table(Index,Area,Slice,Xrev,Yrev,Xcor_rev,Ycor_rev,Zabs,MeanIntMarker,MeanIntLib,PAstate);
clear Index Area Slice Xrev Yrev Xcor_rev Ycor_rev Zabs MeanIntMarker MeanIntLib PAstate ResultsTable Xabs Yabs Xcor Ycor;
rows = (AnalysisTable.Area<EdgeArea | AnalysisTable.Xrev<EdgeDis | AnalysisTable.Xrev>sizeX*PhysicalSize-EdgeDis | AnalysisTable.Yrev<EdgeDis | AnalysisTable.Yrev>sizeX*PhysicalSize-EdgeDis);
AnalysisTable = AnalysisTable(~rows,:);
%% find ROIgate (brightness)
figure();hold;
Xscatter = AnalysisTable.MeanIntMarker;
Yscatter = AnalysisTable.MeanIntLib;
scatter(Xscatter,Yscatter,'Marker','.');
% dscatter(Xscatter,Yscatter,'Marker','.','log',1);
xlabel('Marker brightness');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% xlim([1 100000]);
% ylim([1 100000]);
title('select a gate');
ylabel('Library brightness');
gate = drawpolygon;
set(gate,'userdata',[Xscatter Yscatter]);
position_gate = customWait(gate);
in = inpolygon(log(Xscatter),log(Yscatter),log(position_gate(:,1)),log(position_gate(:,2)));
close(gcf);
figure();
plot(Xscatter(in),Yscatter(in),'r.',Xscatter(~in),Yscatter(~in),'b.');hold on
drawpolygon('Position',position_gate);
xlabel('Marker brightness');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% xlim([1 100000]);
% ylim([1 100000]);
ylabel('Library brightness');
title(gca,['Selected cell num:',num2str(sum(in))]);
path_Analysis = dir_process + "/Analysis";
mkdir(path_Analysis);
saveas(gca,path_Analysis+'/0 gating result.fig');
saveas(gca,path_Analysis+'/0 gating result.png');
PATable = AnalysisTable(in,:);
clear in;
%% generate position list
pl = gui.getPositionList(); % class org.micromanager.PositionList
pl.clearAllPositions();
XYStage_Label = mmc.getXYStageDevice();
PFSStage_Label = 'TIPFSOffset'; % for stage1 in A317 , use 'TIPFSOffset',for demo use 'Z'
if simulate == 1
    PFSStage_Label = 'Z';
end
msp = cell(size(PATable,1),1);
for i = 1:size(PATable,1)
PosiX = position_slice(PATable.Slice(i),1)+PATable.Xcor_rev(i);
PosiY = position_slice(PATable.Slice(i),2)+PATable.Ycor_rev(i);
PosiZ = position_slice(PATable.Slice(i),3);
msp(i,1) = org.micromanager.MultiStagePosition(XYStage_Label,PosiX,PosiY,PFSStage_Label,PosiZ);
msp{i,1}.setLabel(['cell_',num2str(PATable.Index(i))]);
pl.addPosition(msp{i,1});
end
%% save analysis data
save(path_Analysis+'/LibAnalysis.mat','AnalysisTable','PATable');
xlswrite(path_Analysis+'/LibAnalysis.xlsx',[AnalysisTable.Properties.VariableNames;table2cell(AnalysisTable)],'AnalysisTable');
xlswrite(path_Analysis+'/LibAnalysis.xlsx',[PATable.Properties.VariableNames;table2cell(PATable)],'PATable');
%%  run PA_manually
run PA_manually_FP_screen;
%% functions
function pos = customWait(hROI)
title('adjust your gate');
% set(gca, 'YScale', 'linear');
% set(gca, 'XScale', 'linear');
% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);
l = addlistener(hROI,'ROIMoved',@movedCallback);
% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;

end

function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end

function movedCallback(src,evt)
coordinate = get(src,'userdata');
in = inpolygon(log(coordinate(:,1)),log(coordinate(:,2)),log(src.Position(:,1)),log(src.Position(:,2)));
title(gca,['Selected cell num:',num2str(sum(in))]);
end