function GrainAnalyzer_Run(fileDir,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GrainAnalyzer_Run
% This program is used to extract grain number and average length, width, area and area/minBoundRectangle from a batch of images.
% Usage:
%   GrainAnalyzer_Run(fileDir,save_fig,imageRotationAngle)
%      fileDir: full path to images storage Dir, e.g., 'D:/test'
%      save_fig: whether to save the result images (saved in the fileDir folder) for review, =1 for saving; =0 for not saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars -except fileDir save_fig imageRotationAngle
close all;
fclose all;
global hue_low;
global hue_high;
global sat_low;
global edge_plot;
global medfilt_P;
global area_low;
global area_high;
global len_low;
global len_high;
global wid_low;
global wid_high;
global fig_handle;
global enlarge_status;
global pixel_scale;
global erode_depth;
global bw_threshold;
global val_high;
global val_low;
global fid;
global fid2;

global xyList2;
global plotPointHandle;
global plotPointHandle2;
global plotLineHandle;
global currPt;
xyList2=zeros(0,2);
plotPointHandle={};
plotPointHandle2={};
plotLineHandle={};
currPt=zeros(0,2);

hue_low=0;
hue_high=1;
sat_low=0;
medfilt_P=0;
edge_plot=[];
enlarge_status=0;
pixel_scale=0;
erode_depth=0;
bw_threshold=0.5;
val_low=0;
val_high=1;

switch nargin
    case 2
        [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant('grain',fileDir,save_fig);
    case 1
        [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant('grain',fileDir);
    case 0
        [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant('grain');
    otherwise
        error('Too many input parameters are given!');
end

if strcmp(writeMode,'w')
    fprintf(fid, 'fileName\tGrainArea\tGrainLength\tGrainWidth\tAreaRatio\n');
    fprintf(fid2, 'fileName\tGrainNum\tGrainLen_mean\tGrainWid_mean\tGrainArea_mean\tAreaRatio_mean\n');
end
fprintf(fid, '%s\n',strcat('*********',datestr(now,'HH:MM:SS mm/dd/yyyy')));
fprintf(fid2, '%s\n',strcat('*********',datestr(now,'HH:MM:SS mm/dd/yyyy')));

screen_size_4=get(0,'ScreenSize');
screen_height=screen_size_4(4);

input('Press ENTER and then choose a typical image to set grain recognition parameters.','s');
[fileName,filepath] = uigetfile(fullfile(fileDir,'*.jpg'));
if abs(brightness-1)>0.01
    Img=imread(fullfile(filepath,fileName))*brightness;
else
    Img=imread(fullfile(filepath,fileName));
end
size_I=size(Img);
fig_handle=figure();imshow(Img);set(fig_handle,'pointer','cross');ax_handle=gca;hold on
title('Step1: Click-and-drag mouse to crop image, press ENTER to continue');
set(gcf, 'position', get(0,'ScreenSize'));set(gcf,'MenuBar','none');set(gcf,'ToolBar','none');
k = waitforbuttonpress;
plotPointHandle={};
plotPointHandle2={};
ImgCropMat=[1 size_I(1) 1 size_I(2)];
while k~=1
    if k==0
        point1 = get(gca,'CurrentPoint');    % button down detected
        finalRect = rbbox;                   % return figure units
        point2 = get(gca,'CurrentPoint');    % button up detected
        point1 = point1(1,1:2);              % extract x and y
        point2 = point2(1,1:2);
        p1 = min(point1,point2);             % calculate locations
        offset = abs(point1-point2);         % and dimensions
        ImgCropMat=[p1(2) p1(2)+offset(2) p1(1) p1(1)+offset(1)];
        x = round([p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)]);
        y = round([p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)]);
        if ~isempty(plotPointHandle)
            delete(plotPointHandle{end});
            delete(plotPointHandle2{end});
            plotPointHandle(end)=[];
            plotPointHandle2(end)=[];
        end
        plotPointHandle{end+1}=plot(x, y,'r--','LineWidth',1);
        plotPointHandle2{end+1}=copyobj(plotPointHandle{end}, ax_handle);
    end
    k = waitforbuttonpress;
end
Img=Img(ImgCropMat(1):ImgCropMat(2),ImgCropMat(3):ImgCropMat(4),:);
Img_hsv=rgb2hsv(Img);

close;

fig_handle=figure();imshow(Img);hold on
%%%%%%%%%%%%start: if refMark in each image, update pixel_scale %%%%%%%%%%%%
if strcmp(Marker_status,'y') || strcmp(Marker_status,'Y') % pixel_scale==-1
    pixel_scale=10; % initial value
    title('Step1.2: Drag the sliders to recognize the reference marker, press ENTER to proceed');
    axh=gca;
    set(fig_handle,'pointer','arrow');
    set(gcf, 'position', get(0,'ScreenSize'),'MenuBar','none','ToolBar','none');
    
    min_area=1;
    max_area=size(Img,1)*size(Img,2);
    min_len=1;
    max_len=max(size(Img));
    min_wid=1;
    max_wid=max(size(Img)); % all in pixels
    
    h0_text=uicontrol('Style','text',...
        'Position',[5 0.89*screen_height 75 25],...
        'String','BWT','fontsize',10,'backgroundcolor',[0.8 0 0]);
    h0 = uicontrol('style','slider','units','pixel','position',[80 0.89*screen_height 200 25],'Min',0,'Max',1,'Value',0.5);
    h0_2_text=uicontrol('Style','edit',...
        'Position',[285 0.89*screen_height 45 25],...
        'String','0.5','fontsize',10);
    h1_text=uicontrol('Style','text',...
        'Position',[5 0.86*screen_height 75 25],...
        'String','Hue>','fontsize',10,'backgroundcolor','r');
    h1 = uicontrol('style','slider','units','pixel','position',[80 0.86*screen_height 200 25],'Min',0,'Max',1,'Value',0);
    h1_2_text=uicontrol('Style','edit',...
        'Position',[285 0.86*screen_height 45 25],...
        'String','0','fontsize',10);
    h2_text=uicontrol('Style','text',...
        'Position',[5 0.83*screen_height 75 25],...
        'String','Hue<','fontsize',10,'backgroundcolor','r');
    h2 = uicontrol('style','slider','units','pixel','position',[80 0.83*screen_height 200 25],'Min',0,'Max',1,'Value',1);
    h2_2_text=uicontrol('Style','edit',...
        'Position',[285 0.83*screen_height 45 25],...
        'String','1','fontsize',10);
    h3_text=uicontrol('Style','text',...
        'Position',[5 0.8*screen_height 75 25],...
        'String','Sat>','fontsize',10,'backgroundcolor','g');
    h3 = uicontrol('style','slider','units','pixel','position',[80 0.8*screen_height 200 25],'Min',0,'Max',1,'Value',0);
    h3_2_text=uicontrol('Style','edit',...
        'Position',[285 0.8*screen_height 45 25],...
        'String','0','fontsize',10);
    h4_text=uicontrol('Style','text',...
        'Position',[5 0.77*screen_height 75 25],...
        'String','Area>','fontsize',10,'backgroundcolor','y');
    h4 = uicontrol('style','slider','units','pixel','position',[80 0.77*screen_height 200 25],'Min',0,'Max',max_area/2,'Value',min_area); % mm2
    h4_2_text=uicontrol('Style','edit',...
        'Position',[285 0.77*screen_height 45 25],...
        'String',num2str(min_area),'fontsize',10);
    h5_text=uicontrol('Style','text',...
        'Position',[5 0.74*screen_height 75 25],...
        'String','Area<','fontsize',10,'backgroundcolor','y');
    h5 = uicontrol('style','slider','units','pixel','position',[80 0.74*screen_height 200 25],'Min',5,'Max',max_area,'Value',max_area); % mm2
    h5_2_text=uicontrol('Style','edit',...
        'Position',[285 0.74*screen_height 45 25],...
        'String',num2str(max_area),'fontsize',10);
    h6_text=uicontrol('Style','text',...
        'Position',[5 0.71*screen_height 75 25],...
        'String','Len>','fontsize',10,'backgroundcolor','m');
    h6 = uicontrol('style','slider','units','pixel','position',[80 0.71*screen_height 200 25],'Min',0,'Max',max_len/2,'Value',min_len); % mm
    h6_2_text=uicontrol('Style','edit',...
        'Position',[285 0.71*screen_height 45 25],...
        'String',num2str(min_len),'fontsize',10);
    h7_text=uicontrol('Style','text',...
        'Position',[5 0.68*screen_height 75 25],...
        'String','Len<','fontsize',10,'backgroundcolor','m');
    h7 = uicontrol('style','slider','units','pixel','position',[80 0.68*screen_height 200 25],'Min',3,'Max',max_len,'Value',max_len); % mm
    h7_2_text=uicontrol('Style','edit',...
        'Position',[285 0.68*screen_height 45 25],...
        'String',num2str(max_len),'fontsize',10);
    h8_text=uicontrol('Style','text',...
        'Position',[5 0.65*screen_height 75 25],...
        'String','Wid>','fontsize',10,'backgroundcolor','c');
    h8 = uicontrol('style','slider','units','pixel','position',[80 0.65*screen_height 200 25],'Min',0,'Max',max_wid/2,'Value',min_wid); % mm
    h8_2_text=uicontrol('Style','edit',...
        'Position',[285 0.65*screen_height 45 25],...
        'String',num2str(min_wid),'fontsize',10);
    h9_text=uicontrol('Style','text',...
        'Position',[5 0.62*screen_height 75 25],...
        'String','Wid<','fontsize',10,'backgroundcolor','c');
    h9 = uicontrol('style','slider','units','pixel','position',[80 0.62*screen_height 200 25],'Min',2,'Max',max_wid,'Value',max_wid); % mm
    h9_2_text=uicontrol('Style','edit',...
        'Position',[285 0.62*screen_height 45 25],...
        'String',num2str(max_wid),'fontsize',10);
    h10_text=uicontrol('Style','text',...
        'Position',[5 0.59*screen_height 75 25],...
        'String','Smooth','fontsize',10,'backgroundcolor','w');
    h10 = uicontrol('style','slider','units','pixel','position',[80 0.59*screen_height 200 25],'Min',0,'Max',20,'Value',0);
    h10_2_text=uicontrol('Style','edit',...
        'Position',[285 0.59*screen_height 45 25],...
        'String','0','fontsize',10);
    h11_text=uicontrol('Style','text',...
        'Position',[5 0.56*screen_height 75 25],...
        'String','Erode','fontsize',10,'backgroundcolor',[0.5 0.5 0.5]);
    h11 = uicontrol('style','slider','units','pixel','position',[80 0.56*screen_height 200 25],'Min',0,'Max',20,'Value',0);
    h11_2_text=uicontrol('Style','edit',...
        'Position',[285 0.56*screen_height 45 25],...
        'String','0','fontsize',10);
    h12_text=uicontrol('Style','text',...
        'Position',[5 0.53*screen_height 75 25],...
        'String','Val>','fontsize',10,'backgroundcolor','r');
    h12 = uicontrol('style','slider','units','pixel','position',[80 0.53*screen_height 200 25],'Min',0,'Max',1,'Value',0);
    h12_2_text=uicontrol('Style','edit',...
        'Position',[285 0.53*screen_height 45 25],...
        'String','0','fontsize',10);
    h13_text=uicontrol('Style','text',...
        'Position',[5 0.5*screen_height 75 25],...
        'String','Val<','fontsize',10,'backgroundcolor','r');
    h13 = uicontrol('style','slider','units','pixel','position',[80 0.5*screen_height 200 25],'Min',0,'Max',1,'Value',1);
    h13_2_text=uicontrol('Style','edit',...
        'Position',[285 0.5*screen_height 45 25],...
        'String','1','fontsize',10);
    
    area_low=min_area/100*pixel_scale^2; % necessary. ctg 2019/10/29
    area_high=max_area/100*pixel_scale^2; % mm2/100 -->cm2  -> pixels^2
    len_low=min_len/10*pixel_scale;
    len_high=max_len/10*pixel_scale; % mm/10 -->cm -> pixels
    wid_low=min_wid/10*pixel_scale;
    wid_high=max_wid/10*pixel_scale; % mm/10 -->cm -> pixels
    
    set(h0_2_text,'Callback',{@myEditBox_Callback,h0,0,1,Img,Img_hsv,0});
    set(h1_2_text,'Callback',{@myEditBox_Callback,h1,0,1,Img,Img_hsv,1});
    set(h2_2_text,'Callback',{@myEditBox_Callback,h2,0,1,Img,Img_hsv,2});
    set(h3_2_text,'Callback',{@myEditBox_Callback,h3,0,1,Img,Img_hsv,3});
    set(h4_2_text,'Callback',{@myEditBox_Callback,h4,0,max_area/2,Img,Img_hsv,4});
    set(h5_2_text,'Callback',{@myEditBox_Callback,h5,5,max_area,Img,Img_hsv,5});
    set(h6_2_text,'Callback',{@myEditBox_Callback,h6,0,max_len/2,Img,Img_hsv,6});
    set(h7_2_text,'Callback',{@myEditBox_Callback,h7,3,max_len,Img,Img_hsv,7});
    set(h8_2_text,'Callback',{@myEditBox_Callback,h8,0,max_wid/2,Img,Img_hsv,8});
    set(h9_2_text,'Callback',{@myEditBox_Callback,h9,2,max_wid,Img,Img_hsv,9});
    set(h10_2_text,'Callback',{@myEditBox_Callback,h10,0,20,Img,Img_hsv,10});
    set(h11_2_text,'Callback',{@myEditBox_Callback,h11,0,20,Img,Img_hsv,11});
    set(h12_2_text,'Callback',{@myEditBox_Callback,h12,0,20,Img,Img_hsv,12});
    set(h13_2_text,'Callback',{@myEditBox_Callback,h13,0,20,Img,Img_hsv,13});
    try
        %             addlistener(h0, 'ActionEvent','PostSet',{@BinaryAdjust,Img,Img_hsv,0,h0_2_text});
        addlistener(h0, 'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,0,h0_2_text));
        addlistener(h1,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,1,h1_2_text));
        addlistener(h2,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,2,h2_2_text));
        addlistener(h3,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,3,h3_2_text));
        addlistener(h4,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,4,h4_2_text));
        addlistener(h5,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,5,h5_2_text));
        addlistener(h6,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,6,h6_2_text));
        addlistener(h7,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,7,h7_2_text));
        addlistener(h8,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,8,h8_2_text));
        addlistener(h9,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,9,h9_2_text));
        addlistener(h10,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,10,h10_2_text));
        addlistener(h11,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,11,h11_2_text));
        addlistener(h12,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,12,h12_2_text));
        addlistener(h13,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,13,h13_2_text));
    catch
        addlistener(h0,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,0,h0_2_text));
        %             addlistener(h0,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,0,h0_2_text));
        addlistener(h1,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,1,h1_2_text));
        addlistener(h2,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,2,h2_2_text));
        addlistener(h3,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,3,h3_2_text));
        addlistener(h4,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,4,h4_2_text));
        addlistener(h5,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,5,h5_2_text));
        addlistener(h6,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,6,h6_2_text));
        addlistener(h7,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,7,h7_2_text));
        addlistener(h8,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,8,h8_2_text));
        addlistener(h9,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,9,h9_2_text));
        addlistener(h10,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,10,h10_2_text));
        addlistener(h11,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,11,h11_2_text));
        addlistener(h12,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,12,h12_2_text));
        addlistener(h13,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,13,h13_2_text));
    end
    object_zoomin;
    
    while 1
        char_now=input('','s');
        %char_now=get(gcf, 'CurrentCharacter');
        if strcmp(char_now,'')
            fig_handle=figure();imshow(Img);hold on  %% turn on the image
            break
        else
            return;
        end
    end
    area_low_ref=area_low;
    area_high_ref=area_high;
    len_low_ref=len_low;
    len_high_ref=len_high;
    wid_low_ref=wid_low;
    wid_high_ref=wid_high;
    hue_low_ref=hue_low;
    hue_high_ref=hue_high;
    sat_low_ref=sat_low;
    bw_threshold_ref=bw_threshold;
    medfilt_P_ref=medfilt_P;
    erode_depth_ref=erode_depth;
    val_low_ref=val_low;
    val_high_ref=val_high;
    if abs(hue_low_ref)<0.01 && abs(hue_high_ref-1)<0.01 && abs(sat_low_ref)<0.01 && abs(val_low_ref)<0.01 && abs(val_high_ref-1)<0.01
        I_bw=imbinarize(rgb2gray(Img),bw_threshold_ref); %im2bw
        if sum(sum(I_bw))>prod(size(I_bw))/2
            I_bw=imcomplement(I_bw);
        end
    else
        I_hsv=rgb2hsv(Img);
        I_bw=I_hsv(:,:,1)>hue_low_ref & I_hsv(:,:,1)<hue_high_ref & I_hsv(:,:,2)>sat_low_ref & I_hsv(:,:,3)>val_low_ref & I_hsv(:,:,3)<val_high_ref;
    end
    if erode_depth_ref>0
        se = strel('disk',erode_depth_ref);
        I_bw = imerode(I_bw,se);
    end
    if medfilt_P_ref>0
        I_bw=medfilt2(I_bw, [medfilt_P_ref medfilt_P_ref]);
    end
    I_bw = imfill(I_bw,'holes');
    I_bw=bwareafilt(I_bw, [area_low_ref area_high_ref]);
    [I_B,L]=bwboundaries(I_bw);
    num=length(I_B);
    stats = regionprops(L,'Area','Centroid','MajorAxisLength','MinorAxisLength','ConvexArea');
    ref_area=0;
    for i=1:num
        if stats(i).MajorAxisLength>len_low_ref && stats(i).MajorAxisLength<len_high_ref && stats(i).MinorAxisLength>wid_low_ref && stats(i).MinorAxisLength<wid_high_ref
            if stats(i).Area/stats(i).ConvexArea>0.75
                ref_area=max(ref_area,stats(i).Area); % the maximal one is reference marker
            end
        else
            continue;
        end
    end
    pixel_scale=sqrt(ref_area/scale_area); % pixels/cm
    disp(['Calculated scale (pixels/cm): ',num2str(pixel_scale)])
end
%%%%%%%%%%%%end: if refMark in each image, update pixel_scale %%%%%%%%%%%%

title('Step2: Drag the sliders to recognize grains, press ENTER to proceed');
axh=gca;
set(fig_handle,'pointer','arrow');
set(gcf, 'position', get(0,'ScreenSize'),'MenuBar','none','ToolBar','none');

min_area=5;
max_area=100; % mm2, in case there are some connected grains
min_len=3;
max_len=20; % mm, in case there are some connected grains
min_wid=1;
max_wid=20; % mm, in case there are some connected grains
hue_low=0;
hue_high=1;
sat_low=0;
medfilt_P=0;
edge_plot=[];
enlarge_status=0;
erode_depth=0;
bw_threshold=0.5;
val_low=0;
val_high=1;

h0_text=uicontrol('Style','text',...
    'Position',[5 0.89*screen_height 75 25],...
    'String','BWT','fontsize',10,'backgroundcolor',[0.8 0 0]);
h0 = uicontrol('style','slider','units','pixel','position',[80 0.89*screen_height 200 25],'Min',0,'Max',1,'Value',0.5);
h0_2_text=uicontrol('Style','edit',...
    'Position',[285 0.89*screen_height 45 25],...
    'String','0.5','fontsize',10);
h1_text=uicontrol('Style','text',...
    'Position',[5 0.86*screen_height 75 25],...
    'String','Hue>','fontsize',10,'backgroundcolor','r');
h1 = uicontrol('style','slider','units','pixel','position',[80 0.86*screen_height 200 25],'Min',0,'Max',1,'Value',0);
h1_2_text=uicontrol('Style','edit',...
    'Position',[285 0.86*screen_height 45 25],...
    'String','0','fontsize',10);
h2_text=uicontrol('Style','text',...
    'Position',[5 0.83*screen_height 75 25],...
    'String','Hue<','fontsize',10,'backgroundcolor','r');
h2 = uicontrol('style','slider','units','pixel','position',[80 0.83*screen_height 200 25],'Min',0,'Max',1,'Value',1);
h2_2_text=uicontrol('Style','edit',...
    'Position',[285 0.83*screen_height 45 25],...
    'String','1','fontsize',10);
h3_text=uicontrol('Style','text',...
    'Position',[5 0.8*screen_height 75 25],...
    'String','Sat>','fontsize',10,'backgroundcolor','g');
h3 = uicontrol('style','slider','units','pixel','position',[80 0.8*screen_height 200 25],'Min',0,'Max',1,'Value',0);
h3_2_text=uicontrol('Style','edit',...
    'Position',[285 0.8*screen_height 45 25],...
    'String','0','fontsize',10);
h4_text=uicontrol('Style','text',...
    'Position',[5 0.77*screen_height 75 25],...
    'String','Area>','fontsize',10,'backgroundcolor','y');
h4 = uicontrol('style','slider','units','pixel','position',[80 0.77*screen_height 200 25],'Min',0,'Max',20,'Value',min_area); % mm2
h4_2_text=uicontrol('Style','edit',...
    'Position',[285 0.77*screen_height 45 25],...
    'String',num2str(min_area),'fontsize',10);
h5_text=uicontrol('Style','text',...
    'Position',[5 0.74*screen_height 75 25],...
    'String','Area<','fontsize',10,'backgroundcolor','y');
h5 = uicontrol('style','slider','units','pixel','position',[80 0.74*screen_height 200 25],'Min',5,'Max',200,'Value',max_area); % mm2
h5_2_text=uicontrol('Style','edit',...
    'Position',[285 0.74*screen_height 45 25],...
    'String',num2str(max_area),'fontsize',10);
h6_text=uicontrol('Style','text',...
    'Position',[5 0.71*screen_height 75 25],...
    'String','Len>','fontsize',10,'backgroundcolor','m');
h6 = uicontrol('style','slider','units','pixel','position',[80 0.71*screen_height 200 25],'Min',0,'Max',5,'Value',min_len); % mm
h6_2_text=uicontrol('Style','edit',...
    'Position',[285 0.71*screen_height 45 25],...
    'String',num2str(min_len),'fontsize',10);
h7_text=uicontrol('Style','text',...
    'Position',[5 0.68*screen_height 75 25],...
    'String','Len<','fontsize',10,'backgroundcolor','m');
h7 = uicontrol('style','slider','units','pixel','position',[80 0.68*screen_height 200 25],'Min',3,'Max',50,'Value',max_len); % mm
h7_2_text=uicontrol('Style','edit',...
    'Position',[285 0.68*screen_height 45 25],...
    'String',num2str(max_len),'fontsize',10);
h8_text=uicontrol('Style','text',...
    'Position',[5 0.65*screen_height 75 25],...
    'String','Wid>','fontsize',10,'backgroundcolor','c');
h8 = uicontrol('style','slider','units','pixel','position',[80 0.65*screen_height 200 25],'Min',0,'Max',3,'Value',min_wid); % mm
h8_2_text=uicontrol('Style','edit',...
    'Position',[285 0.65*screen_height 45 25],...
    'String',num2str(min_wid),'fontsize',10);
h9_text=uicontrol('Style','text',...
    'Position',[5 0.62*screen_height 75 25],...
    'String','Wid<','fontsize',10,'backgroundcolor','c');
h9 = uicontrol('style','slider','units','pixel','position',[80 0.62*screen_height 200 25],'Min',2,'Max',30,'Value',max_wid); % mm
h9_2_text=uicontrol('Style','edit',...
    'Position',[285 0.62*screen_height 45 25],...
    'String',num2str(max_wid),'fontsize',10);
h10_text=uicontrol('Style','text',...
    'Position',[5 0.59*screen_height 75 25],...
    'String','Smooth','fontsize',10,'backgroundcolor','w');
h10 = uicontrol('style','slider','units','pixel','position',[80 0.59*screen_height 200 25],'Min',0,'Max',20,'Value',0);
h10_2_text=uicontrol('Style','edit',...
    'Position',[285 0.59*screen_height 45 25],...
    'String','0','fontsize',10);
h11_text=uicontrol('Style','text',...
    'Position',[5 0.56*screen_height 75 25],...
    'String','Erode','fontsize',10,'backgroundcolor',[0.5 0.5 0.5]);
h11 = uicontrol('style','slider','units','pixel','position',[80 0.56*screen_height 200 25],'Min',0,'Max',20,'Value',0);
h11_2_text=uicontrol('Style','edit',...
    'Position',[285 0.56*screen_height 45 25],...
    'String','0','fontsize',10);
h12_text=uicontrol('Style','text',...
    'Position',[5 0.53*screen_height 75 25],...
    'String','Val>','fontsize',10,'backgroundcolor','r');
h12 = uicontrol('style','slider','units','pixel','position',[80 0.53*screen_height 200 25],'Min',0,'Max',1,'Value',0);
h12_2_text=uicontrol('Style','edit',...
    'Position',[285 0.53*screen_height 45 25],...
    'String','0','fontsize',10);
h13_text=uicontrol('Style','text',...
    'Position',[5 0.5*screen_height 75 25],...
    'String','Val<','fontsize',10,'backgroundcolor','r');
h13 = uicontrol('style','slider','units','pixel','position',[80 0.5*screen_height 200 25],'Min',0,'Max',1,'Value',1);
h13_2_text=uicontrol('Style','edit',...
    'Position',[285 0.5*screen_height 45 25],...
    'String','1','fontsize',10);

area_low=min_area/100*pixel_scale^2;
area_high=max_area/100*pixel_scale^2; % mm2/100 -->cm2  -> pixels^2
len_low=min_len/10*pixel_scale;
len_high=max_len/10*pixel_scale; % mm/10 -->cm -> pixels
wid_low=min_wid/10*pixel_scale;
wid_high=max_wid/10*pixel_scale; % mm/10 -->cm -> pixels

set(h0_2_text,'Callback',{@myEditBox_Callback,h0,0,1,Img,Img_hsv,0});
set(h1_2_text,'Callback',{@myEditBox_Callback,h1,0,1,Img,Img_hsv,1});
set(h2_2_text,'Callback',{@myEditBox_Callback,h2,0,1,Img,Img_hsv,2});
set(h3_2_text,'Callback',{@myEditBox_Callback,h3,0,1,Img,Img_hsv,3});
set(h4_2_text,'Callback',{@myEditBox_Callback,h4,0,max_area/2,Img,Img_hsv,4});
set(h5_2_text,'Callback',{@myEditBox_Callback,h5,5,max_area,Img,Img_hsv,5});
set(h6_2_text,'Callback',{@myEditBox_Callback,h6,0,max_len/2,Img,Img_hsv,6});
set(h7_2_text,'Callback',{@myEditBox_Callback,h7,3,max_len,Img,Img_hsv,7});
set(h8_2_text,'Callback',{@myEditBox_Callback,h8,0,max_wid/2,Img,Img_hsv,8});
set(h9_2_text,'Callback',{@myEditBox_Callback,h9,2,max_wid,Img,Img_hsv,9});
set(h10_2_text,'Callback',{@myEditBox_Callback,h10,0,20,Img,Img_hsv,10});
set(h11_2_text,'Callback',{@myEditBox_Callback,h11,0,20,Img,Img_hsv,11});
set(h12_2_text,'Callback',{@myEditBox_Callback,h12,0,20,Img,Img_hsv,12});
set(h13_2_text,'Callback',{@myEditBox_Callback,h13,0,20,Img,Img_hsv,13});
try
    addlistener(h0,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,0,h0_2_text));
    addlistener(h1,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,1,h1_2_text));
    addlistener(h2,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,2,h2_2_text));
    addlistener(h3,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,3,h3_2_text));
    addlistener(h4,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,4,h4_2_text));
    addlistener(h5,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,5,h5_2_text));
    addlistener(h6,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,6,h6_2_text));
    addlistener(h7,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,7,h7_2_text));
    addlistener(h8,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,8,h8_2_text));
    addlistener(h9,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,9,h9_2_text));
    addlistener(h10,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,10,h10_2_text));
    addlistener(h11,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,11,h11_2_text));
    addlistener(h12,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,12,h12_2_text));
    addlistener(h13,'ActionEvent',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,13,h13_2_text));
catch
    addlistener(h0,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,0,h0_2_text));
    addlistener(h1,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,1,h1_2_text));
    addlistener(h2,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,2,h2_2_text));
    addlistener(h3,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,3,h3_2_text));
    addlistener(h4,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,4,h4_2_text));
    addlistener(h5,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,5,h5_2_text));
    addlistener(h6,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,6,h6_2_text));
    addlistener(h7,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,7,h7_2_text));
    addlistener(h8,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,8,h8_2_text));
    addlistener(h9,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,9,h9_2_text));
    addlistener(h10,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,10,h10_2_text));
    addlistener(h11,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,11,h11_2_text));
    addlistener(h12,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,12,h12_2_text));
    addlistener(h13,'ContinuousValueChange',@(hObject, ContinuousValueChange) BinaryAdjust(hObject,Img,Img_hsv,13,h13_2_text));
end

object_zoomin;

while 1
    char_now=input('','s');
    %char_now=get(gcf, 'CurrentCharacter');
    if strcmp(char_now,'')
        break
    else
        return;
    end
end

figList=dir(fullfile(fileDir,'*.jpg'));
%% global parameters
RGB=1;
for pic = 1:1:length(figList)
    clearvars -except fileDir figList pic RGB save_fig ImgCropMat fid fid2 hue_low hue_high sat_low area_low area_high medfilt_P pixel_scale len_low len_high wid_low wid_high bw_threshold erode_depth val_low val_high area_low_ref area_high_ref len_low_ref len_high_ref wid_low_ref wid_high_ref hue_low_ref hue_high_ref sat_low_ref bw_threshold_ref erode_depth_ref medfilt_P_ref val_low_ref val_high_ref scale_area Marker_status brightness
    fn=figList(pic).name;
    figNA=fn(1:end-4);
    fprintf('%s, the %d/%d image is in processing...\n',fn,pic,length(figList))
    fprintf(fid2,fn);
    fprintf(fid,fn);
    if abs(brightness-1)>0.01
        I_Raw=imread(fullfile(fileDir,fn))*brightness;
    else
        I_Raw=imread(fullfile(fileDir,fn));
    end
    if ~isempty(ImgCropMat)
        I_Raw=I_Raw(ImgCropMat(1):ImgCropMat(2),ImgCropMat(3):ImgCropMat(4),:);
    end
    x_max=ImgCropMat(2)-ImgCropMat(1);
    y_max=ImgCropMat(4)-ImgCropMat(3);
    %%%%%%%%%%%%start: if refMark in each image, update pixel_scale %%%%%%%%%%%%
    if strcmp(Marker_status,'y') || strcmp(Marker_status,'Y')
        if abs(hue_low_ref)<0.01 && abs(hue_high_ref-1)<0.01 && abs(sat_low_ref)<0.01 && abs(val_low_ref)<0.01 && abs(val_high_ref-1)<0.01
            I_bw=imbinarize(rgb2gray(I_Raw),bw_threshold_ref); %im2bw
            if sum(sum(I_bw))>prod(size(I_bw))/2
                I_bw=imcomplement(I_bw);
            end
        else
            I_hsv=rgb2hsv(I_Raw);
            I_bw=I_hsv(:,:,1)>hue_low_ref & I_hsv(:,:,1)<hue_high_ref & I_hsv(:,:,2)>sat_low_ref & I_hsv(:,:,3)>val_low_ref & I_hsv(:,:,3)<val_high_ref;
        end
        
        if erode_depth_ref>0
            se = strel('disk',erode_depth_ref);
            I_bw = imerode(I_bw,se);
        end
        if medfilt_P_ref>0
            I_bw=medfilt2(I_bw, [medfilt_P_ref medfilt_P_ref]);
        end
        I_bw = imfill(I_bw,'holes');
        I_bw=bwareafilt(I_bw, [area_low_ref area_high_ref]);
        [I_B,L]=bwboundaries(I_bw);
        num=length(I_B);
        stats = regionprops(L,'Area','Centroid','MajorAxisLength','MinorAxisLength','ConvexArea');
        ref_area=0;
        ref_ind=0;
        for i=1:num
            if stats(i).MajorAxisLength>len_low_ref && stats(i).MajorAxisLength<len_high_ref && stats(i).MinorAxisLength>wid_low_ref && stats(i).MinorAxisLength<wid_high_ref
                if stats(i).Area>ref_area && stats(i).Area/stats(i).ConvexArea>0.75
                    ref_area=stats(i).Area; % the maximal one is reference marker
                    ref_ind=i;
                end
            else
                continue;
            end
        end
        pixel_scale=sqrt(ref_area/scale_area); % pixels/cm
    end
    %%%%%%%%%%%%end: if refMark in each image, update pixel_scale %%%%%%%%%%%%
    if (save_fig==1)
        fh=figure;set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')');
        imshow(I_Raw);hold on;set(gcf, 'Position', get(0, 'ScreenSize'));
        if strcmp(Marker_status,'y') || strcmp(Marker_status,'Y')
            bound=I_B{ref_ind};
            temp = stats(ref_ind).Centroid;
            plot(bound(:,2),bound(:,1),'y--','LineWidth',2);
            text(temp(1), temp(2), 'RefMarker', 'color', 'y','fontsize',10,'FontWeight','b');
        end
    end
    
    if hue_low<0.01 && hue_high>0.99 && sat_low<0.01 && val_low<0.01 && val_high>0.99
        I_bw=imbinarize(rgb2gray(I_Raw),bw_threshold); %im2bw
        if sum(sum(I_bw))>prod(size(I_bw))/2
            I_bw=imcomplement(I_bw);
        end
    else
        I_hsv=rgb2hsv(I_Raw);
        I_bw=I_hsv(:,:,1)>hue_low & I_hsv(:,:,1)<hue_high & I_hsv(:,:,2)>sat_low & I_hsv(:,:,3)>val_low & I_hsv(:,:,3)<val_high;
    end
    if erode_depth>0
        se = strel('disk',erode_depth);
        I_bw = imerode(I_bw,se);
    end
    if medfilt_P>0
        I_bw=medfilt2(I_bw, [medfilt_P medfilt_P]);
    end
    I_bw = imfill(I_bw,'holes');
    I_bw=bwareafilt(I_bw, [area_low area_high]);
    [I_B,L]=bwboundaries(I_bw);
    num=length(I_B);
    stats = regionprops(L,'Area','Centroid','MajorAxisLength','MinorAxisLength');
    effect_ind=[];
    raw_areaG_list=[];
    raw_lenG_list=[];
    raw_widG_list=[];
    for i=1:num
        len_temp=stats(i).MajorAxisLength;
        wid_temp=stats(i).MinorAxisLength;
        center_pos=stats(i).Centroid;
        perimeter_temp=size(I_B{i},1);
        if len_temp>len_low && len_temp<len_high && wid_temp>wid_low && wid_temp<wid_high && center_pos(1)>=y_max*0.01 && center_pos(1)<=y_max*0.99 && center_pos(2)>=x_max*0.01 ...
                && center_pos(2)<=x_max*0.99 && stats(i).Area>0.45*len_temp*wid_temp && perimeter_temp<2.5*(len_temp+wid_temp)
            effect_ind(end+1)=i;
            raw_areaG_list(end+1)=stats(i).Area;
            raw_lenG_list(end+1)=len_temp;
            raw_widG_list(end+1)=wid_temp;
        end
    end
    raw_aveLen=median(raw_lenG_list);
    raw_aveWid=median(raw_widG_list);
    raw_aveArea=median(raw_areaG_list);
    
    grain_count=0;
    areaG_list=[];
    lenG_list=[];
    widG_list=[];
    for i=effect_ind
        len_temp=stats(i).MajorAxisLength;
        wid_temp=stats(i).MinorAxisLength;
        area_temp=stats(i).Area;
        grainNum=max(1,round(area_temp/raw_aveArea)); % multiple grains connect
        grain_count=grain_count+grainNum;
        areaG_list(end+1)=stats(i).Area;
        [lenG_list(end+1),widG_list(end+1)]=lenWidCal2(I_B{i},0);
        if (save_fig==1)
            bound=I_B{i};
            temp = stats(i).Centroid;
            if grainNum==1
                plot(bound(:,2),bound(:,1),'r-','LineWidth',0.1);
                text(temp(1), temp(2), num2str(grain_count), 'color', 'r','FontSize',10);
            else
                plot(bound(:,2),bound(:,1),'b-','LineWidth',1);
                text(temp(1), temp(2), [num2str(grain_count),'(',num2str(grainNum),')'], 'color', 'b','FontSize',12,'FontWeight','b');
            end
        end
        fprintf(fid, '\t%g',[grain_count,lenG_list(end)/pixel_scale*10,widG_list(end)/pixel_scale*10,stats(i).Area/pixel_scale^2*100,stats(i).Area/(lenG_list(end)*widG_list(end))]);
        fprintf(fid, '\n');
    end
    lenAdj=median(lenG_list)/pixel_scale*10; % mm
    widAdj=median(widG_list)/pixel_scale*10; % mm
    areaAdj=median(areaG_list)/pixel_scale^2*100; % mm2
    areaRatioAdj=median(areaG_list./(lenG_list.*widG_list)); % %
    if (save_fig==1)
        title(sprintf('Grain number: %d, length: %.2f mm, width: %.2f mm, area: %.2f mm2, area ratio: %.3f',grain_count,lenAdj,widAdj,areaAdj,areaRatioAdj), 'FontWeight', 'Bold');
        print(fh,fullfile(fileDir,[figNA,'-track']),'-dpng','-r100');
        close;
    end
    fprintf(fid2, '\t%g',[grain_count,lenAdj,widAdj,areaAdj,areaRatioAdj]);
    fprintf(fid2, '\n');
end

fclose(fid);
fclose(fid2);
warning on;
end

function BinaryAdjust(hObject,I_RGB_in,I_HSV_in,bar_num,h_text)
global hue_low;
global hue_high;
global sat_low;
global edge_plot;
global medfilt_P;
global area_low;
global area_high;
global pixel_scale;
global len_low;
global len_high;
global wid_low;
global wid_high;
global erode_depth;
global bw_threshold;
global val_high;
global val_low;

value_update = get(hObject,'Value');
if bar_num==0
    bw_threshold = value_update;
elseif bar_num==1
    hue_low = value_update;
elseif bar_num==2
    hue_high = value_update;
elseif bar_num==3
    sat_low = value_update;
elseif bar_num==4
    area_low = value_update;
    area_low = area_low/100*pixel_scale^2; %pixel^2
elseif bar_num==5
    area_high = value_update;
    area_high = area_high/100*pixel_scale^2; %pixel^2
elseif bar_num==6
    len_low = value_update;
    len_low = len_low/10*pixel_scale; %pixel
elseif bar_num==7
    len_high = value_update;
    len_high = len_high/10*pixel_scale; %pixel
elseif bar_num==8
    wid_low = value_update;
    wid_low = wid_low/10*pixel_scale; %pixel
elseif bar_num==9
    wid_high = value_update;
    wid_high = wid_high/10*pixel_scale; %pixel
elseif bar_num==10
    medfilt_P = round(value_update);
elseif bar_num==11
    erode_depth = round(value_update);
elseif bar_num==12
    val_low = value_update;
elseif bar_num==13
    val_high = value_update;
end

size_I=size(I_RGB_in);
if hue_low<0.01 && hue_high>0.99 && sat_low<0.01 && val_low<0.01 && val_high>0.99
    I_bw=imbinarize(rgb2gray(I_RGB_in),bw_threshold); % im2bw
    if sum(sum(I_bw))>size_I(1)*size_I(2)/2
        I_bw=imcomplement(I_bw); % convert 0->1 and 1->0
    end
else
    I_bw=I_HSV_in(:,:,1)>hue_low & I_HSV_in(:,:,1)<hue_high & I_HSV_in(:,:,2)>sat_low & I_HSV_in(:,:,3)>val_low & I_HSV_in(:,:,3)<val_high;
end

if erode_depth>0
    se = strel('disk',erode_depth);
    I_bw = imerode(I_bw,se);
end
if medfilt_P>0
    I_bw=medfilt2(I_bw, [medfilt_P medfilt_P]);
end
I_bw = imfill(I_bw,'holes');
I_bw=bwareafilt(I_bw, [area_low area_high]);
[B,L]=bwboundaries(I_bw);
stats = regionprops(L,'MajorAxisLength','MinorAxisLength');
edge_plot_num=length(edge_plot);
if edge_plot_num
    if (ishandle(edge_plot(1)))
        delete(edge_plot)
    end
end
edge_plot=[];

for i=1:length(B)
    perimeter_temp=size(B{i},1);
    if stats(i).MajorAxisLength>len_low && stats(i).MajorAxisLength<len_high && stats(i).MinorAxisLength>wid_low && stats(i).MinorAxisLength<wid_high
        bound=B{i};edge_plot=[edge_plot,plot(bound(:,2),bound(:,1),'b-','LineWidth',0.1)];
    end
end
%     [hue_low hue_high sat_low medfilt_P erode_depth bw_threshold]
%     [area_low area_high]/pixel_scale^2*100
%     [len_low len_high wid_low wid_high]/pixel_scale*10
set(h_text,'String',num2str(round(value_update*100)/100));
end

function myEditBox_Callback(hObject,eventdata,h_slider,val_min,val_max,Img,Img_hsv,ui_number)
val = str2double(get(hObject,'String'));
val = max(min(val,val_max),val_min);
set(h_slider,'Value',val);
guidata(hObject,h_slider);
BinaryAdjust(h_slider, Img,Img_hsv,ui_number,hObject)
end

function [len,wid]=lenWidCal2(bound,PLOT)
% lenWidCal2: Compute the minimal bounding rectangle of points in the plane
% based on longest points pair
% usage:
% [len,wid,l_start,l_end,w_start,w_end]=lenWidCal2(bound,convexHull,sizeL)
%
% arguments: (input)
%  bound - vectors of points, describing edge points of the given region,
%  notice that bound(:,2) is the x part, and bound(:,1) is the y part
%
% PLOT - for plot the detailed figure, 1 on, 0 off.
%
%
% arguments: (output)
%  len,wid - length and width of the region
%
%  l_start,l_end,w_start,w_end - position of longest points pair and the
%  other two points that determines the rectangle
%
%
% Note:
%
% Thanks to Xiaojuan for her helpful advice on this problem.
%
%
% Important bug report:
%
% this version uses the 2-seperation method, however, it would go wrong in
% some cases (such as for  parallelogram), another slower version uses
% ergodic is lenWidCal4.m
% seems no bug, -_-... goes wrong just because i didnot notice asin in [-pi/2,pi/2]

x=bound(:,2);
y=bound(:,1);
dt = DelaunayTri(x,y);
k = convexHull(dt);
convex_hull=[x(k),y(k)];
startP=3;
totalP=length(convex_hull(:,1))-1;
len=0;
for i=1:totalP
    endP=i-1+totalP;
    while 1
        if endP-startP==1
            firstP=mod(startP-1,totalP)+1;
            secondP=mod(endP-1,totalP)+1;
            neighborP=mod(i+1-1,totalP)+1;
            area1=abs(det([convex_hull(secondP,:)-convex_hull(i,:);convex_hull(neighborP,:)-convex_hull(i,:)]));
            area2=abs(det([convex_hull(firstP,:)-convex_hull(i,:);convex_hull(neighborP,:)-convex_hull(i,:)]));
            if area1<area2
                firstD=pdist([convex_hull(firstP,:);convex_hull(i,:)],'euclidean');
                secondD=pdist([convex_hull(firstP,:);convex_hull(neighborP,:)],'euclidean');
                lenNew=max(len,max(firstD,secondD));
                if lenNew==firstD
                    pairL=[i,firstP];
                else if lenNew==secondD
                        pairL=[neighborP,firstP];
                    end
                end
                len=lenNew;
            else
                firstD=pdist([convex_hull(secondP,:);convex_hull(i,:)],'euclidean');
                secondD=pdist([convex_hull(secondP,:);convex_hull(neighborP,:)],'euclidean');
                lenNew=max(len,max(firstD,secondD));
                if lenNew==firstD
                    pairL=[i,secondP];
                else if lenNew==secondD
                        pairL=[neighborP,secondP];
                    end
                end
                len=lenNew;
                startP=endP;
            end
            break;
        else
            breakP=startP+round((endP-startP)/2);
            firstP=mod(breakP-1-1,totalP)+1;
            secondP=mod(breakP-1,totalP)+1;
            neighborP=mod(i+1-1,totalP)+1;
            area1=abs(det([convex_hull(secondP,:)-convex_hull(i,:);convex_hull(neighborP,:)-convex_hull(i,:)]));
            area2=abs(det([convex_hull(firstP,:)-convex_hull(i,:);convex_hull(neighborP,:)-convex_hull(i,:)]));
            if area1<area2
                endP=breakP;
            else
                startP=breakP;
            end
        end
    end
end
l_start=min(pairL);
l_end=max(pairL);

x=x';
y=y';
mainAxis=convex_hull(l_end,:)-convex_hull(l_start,:);
sinAlpha=det([1,0;mainAxis])/norm(mainAxis);
cosAlpha=[1,0]*mainAxis'/norm(mainAxis);
if sinAlpha>0 && cosAlpha>0
    alpha=asin(sinAlpha);
else if sinAlpha>0 && cosAlpha<0
        alpha=-asin(sinAlpha);
    else if sinAlpha<0 && cosAlpha>0
            alpha=asin(sinAlpha);
        else
            alpha=-asin(sinAlpha);
        end
    end
end
xNew=cos(alpha)*x+sin(alpha)*y;
yNew=-sin(alpha)*x+cos(alpha)*y;
xNew=xNew';
yNew=yNew';

[xNewMin1,index1]=min(xNew);
[xNewMax1,index2]=max(xNew);
[yNewMin2,index3]=min(yNew);
[yNewMax2,index4]=max(yNew);

xMin1=x(index1);
yMin1=y(index1);
xMin2=x(index3);
yMin2=y(index3);
xMax1=x(index2);
yMax1=y(index2);
xMax2=x(index4);
yMax2=y(index4);

w_start=[xMin2,yMin2];
w_end=[xMax2,yMax2];
wid=yNewMax2-yNewMin2;

if PLOT==1
    %     close all;
    figure;plot(x,y,'.','MarkerSize',10)
    axis equal tight;
    exalaration=100;
    axis([min(xMax2,min(xMin2,min(xMin1,xMax1)))-exalaration max(xMax2,max(xMin2,max(xMin1,xMax1)))+exalaration ...
        min(yMax2,min(yMin2,min(yMin1,yMax1)))-exalaration max(yMax2,max(yMin2,max(yMin1,yMax1)))+exalaration]);
    
    line([xMin1,xMax1],[yMin1,yMax1],'LineWidth',2,'Color',[0,1,1]);
    line([xMin2,xMax2],[yMin2,yMax2],'LineWidth',2,'Color',[0,0,1]);
    % line([xMin2,xMax22],[yMin2,yMax22],'LineWidth',2,'Color',[0,1,1]);
    edge1=[xMin2-xMin1,yMin2-yMin1];
    dist1=det([mainAxis;edge1])/(norm(mainAxis));
    rect1=[xMin1,yMin1]+dist1*[-mainAxis(2),mainAxis(1)]/norm(mainAxis);
    rect2=[xMax1,yMax1]+dist1*[-mainAxis(2),mainAxis(1)]/norm(mainAxis);
    edge2=[xMax2-xMin1,yMax2-yMin1];
    dist2=det([mainAxis;edge2])/(norm(mainAxis));
    rect3=[xMin1,yMin1]+dist2*[-mainAxis(2),mainAxis(1)]/norm(mainAxis);
    rect4=[xMax1,yMax1]+dist2*[-mainAxis(2),mainAxis(1)]/norm(mainAxis);
    line([rect1(1),rect2(1),rect4(1),rect3(1),rect1(1)],[rect1(2),rect2(2),rect4(2),rect3(2),rect1(2)],'LineWidth',3,'Color',[1,1,0]);
    line(convex_hull(:,1),convex_hull(:,2),'LineWidth',3,'Color',[0,1,0])
    text(xMin1,yMin1,'xmin1','FontSize',18)
    text(xMax1,yMax1,'xmax1','FontSize',18)
    text(xMin2,yMin2,'xmin2','FontSize',18)
    text(xMax2,yMax2,'xmax2','FontSize',18)
    set(gcf, 'Position', get(0, 'ScreenSize'));
    %     error('break')
end

end
