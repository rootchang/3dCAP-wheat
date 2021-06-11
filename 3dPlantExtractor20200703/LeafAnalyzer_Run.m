function LeafAnalyzer_Run(fileDir,save_fig,imageRotationAngle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LeafAnalyzer_Run.m
% This program is used to extract leaf number and 2D shape from a batch of images.
% Usage: 
%   LeafAnalyzer_Run(fileDir,save_fig,imageRotationAngle)
%      fileDir: full path to images storage Dir, e.g., 'D:/test'
%      save_fig: whether to save the result images (saved in the fileDir folder) for review, =1 for saving; =0 for not saving
%      imageRotationAngle: Rotation angle (0, 90,180 or 270 degrees) of
%        images for downward-orientation of leaves
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

global gray2bwThreth;
gray2bwThreth=0.5;


switch nargin
    case 3
        [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant('leaf',fileDir,save_fig);
        imageRotationAngle=num2str(imageRotationAngle);
    case 2
        [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant('leaf',fileDir,save_fig);
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are downward oriented (default=0): ','s');
    case 1
        [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant('leaf',fileDir);
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are downward oriented (default=0): ','s');
    case 0
        [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant('leaf');
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are downward oriented (default=0): ','s');
    otherwise
        error('Too many input parameters given!');
end

wrong_input=1;
while wrong_input
    if isempty(imageRotationAngle)
        imageRotationAngle='0';
    end
    try
        imageRotationAngle=str2num(imageRotationAngle);
        if abs(imageRotationAngle-0)>0.01 && abs(imageRotationAngle-90)>0.01 && abs(imageRotationAngle-180)>0.01 && abs(imageRotationAngle-270)>0.01
            imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are downward oriented (default=0): ','s');
        else
            wrong_input=0;
        end
    catch
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are downward oriented (default=0): ','s');
    end
end

if strcmp(writeMode,'w')
    fprintf(fid2, 'fileName\n\tLeafArea\tLeafLength\tWidth1\tWidth2\tWidth3\tWidth4\tWidth5\tWidth6\tWidth7\tWidthMax\n');
end
fprintf(fid, '### Image brightness adjustment coef.: %f\n',brightness); %fileDir,save_fig,,Marker_status,,writeMode
fprintf(fid, '### Image rotation angle: %d\n',imageRotationAngle);
if scale_area>0
    fprintf(fid, '### Area of scale marker in image: %d\n',scale_area);
end


screen_size_4=get(0,'ScreenSize');
screen_height=screen_size_4(4);

input('Please click ENTER and choose a typical image for setting leaves recognition parameters.','s');
[fileName,filepath] = uigetfile(fullfile(fileDir,'*.jpg'));
if abs(brightness-1)>0.01
    Img=imread(fullfile(filepath,fileName))*brightness;
else
    Img=imread(fullfile(filepath,fileName));
end
if abs(imageRotationAngle-90)<0.01
    Img=rot90(Img,1);
elseif abs(imageRotationAngle-180)<0.01
    Img=rot90(Img,2);
elseif abs(imageRotationAngle-270)<0.01
    Img=rot90(Img,3);
end

size_I=size(Img);
fig_handle=figure();imshow(Img);set(fig_handle,'pointer','cross');ax_handle=gca;hold on
title('Step1: Click-and-drag mouse to crop image, press ENTER to proceed');
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
    
    size_I=size(Img);
    area_lower_limit=size_I(1)*size_I(2)/2;
    area_lower_default=1;
    area_upper_limit=size_I(1)*size_I(2);
    area_upper_default=size_I(1)*size_I(2);
    len_lower_limit=max(size_I)/2;
    len_lower_default=1;
    len_upper_limit=max(size_I);
    len_upper_default=max(size_I);
    wid_lower_limit=max(size_I)/2;
    wid_lower_default=1;
    wid_upper_limit=max(size_I);
    wid_upper_default=max(size_I);
    
    
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
    h4 = uicontrol('style','slider','units','pixel','position',[80 0.77*screen_height 200 25],'Min',0,'Max',area_lower_limit,'Value',area_lower_default); % mm2
    h4_2_text=uicontrol('Style','edit',...
        'Position',[285 0.77*screen_height 45 25],...
        'String',num2str(area_lower_default),'fontsize',10);
    h5_text=uicontrol('Style','text',...
        'Position',[5 0.74*screen_height 75 25],...
        'String','Area<','fontsize',10,'backgroundcolor','y');
    h5 = uicontrol('style','slider','units','pixel','position',[80 0.74*screen_height 200 25],'Min',0,'Max',area_upper_limit,'Value',area_upper_default); % mm2
    h5_2_text=uicontrol('Style','edit',...
        'Position',[285 0.74*screen_height 45 25],...
        'String',num2str(area_upper_default),'fontsize',10);
    h6_text=uicontrol('Style','text',...
        'Position',[5 0.71*screen_height 75 25],...
        'String','Len>','fontsize',10,'backgroundcolor','m');
    h6 = uicontrol('style','slider','units','pixel','position',[80 0.71*screen_height 200 25],'Min',0,'Max',len_lower_limit,'Value',len_lower_default); % mm
    h6_2_text=uicontrol('Style','edit',...
        'Position',[285 0.71*screen_height 45 25],...
        'String',num2str(len_lower_default),'fontsize',10);
    h7_text=uicontrol('Style','text',...
        'Position',[5 0.68*screen_height 75 25],...
        'String','Len<','fontsize',10,'backgroundcolor','m');
    h7 = uicontrol('style','slider','units','pixel','position',[80 0.68*screen_height 200 25],'Min',0,'Max',len_upper_limit,'Value',len_upper_default); % mm
    h7_2_text=uicontrol('Style','edit',...
        'Position',[285 0.68*screen_height 45 25],...
        'String',num2str(len_upper_default),'fontsize',10);
    h8_text=uicontrol('Style','text',...
        'Position',[5 0.65*screen_height 75 25],...
        'String','Wid>','fontsize',10,'backgroundcolor','c');
    h8 = uicontrol('style','slider','units','pixel','position',[80 0.65*screen_height 200 25],'Min',0,'Max',wid_lower_limit,'Value',wid_lower_default); % mm
    h8_2_text=uicontrol('Style','edit',...
        'Position',[285 0.65*screen_height 45 25],...
        'String',num2str(wid_lower_default),'fontsize',10);
    h9_text=uicontrol('Style','text',...
        'Position',[5 0.62*screen_height 75 25],...
        'String','Wid<','fontsize',10,'backgroundcolor','c');
    h9 = uicontrol('style','slider','units','pixel','position',[80 0.62*screen_height 200 25],'Min',0,'Max',wid_upper_limit,'Value',wid_upper_default); % mm
    h9_2_text=uicontrol('Style','edit',...
        'Position',[285 0.62*screen_height 45 25],...
        'String',num2str(wid_upper_default),'fontsize',10);
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
    
    area_low=area_lower_default*pixel_scale^2;
    area_high=area_upper_default*pixel_scale^2; % cm2  -> pixels^2
    len_low=len_lower_default*pixel_scale;
    len_high=len_upper_default*pixel_scale; % cm -> pixels
    wid_low=wid_lower_default*pixel_scale;
    wid_high=wid_upper_default*pixel_scale; % cm -> pixels
    
    set(h0_2_text,'Callback',{@myEditBox_Callback,h0,0,1,Img,Img_hsv,0});
    set(h1_2_text,'Callback',{@myEditBox_Callback,h1,0,1,Img,Img_hsv,1});
    set(h2_2_text,'Callback',{@myEditBox_Callback,h2,0,1,Img,Img_hsv,2});
    set(h3_2_text,'Callback',{@myEditBox_Callback,h3,0,1,Img,Img_hsv,3});
    set(h4_2_text,'Callback',{@myEditBox_Callback,h4,0,area_lower_limit,Img,Img_hsv,4});
    set(h5_2_text,'Callback',{@myEditBox_Callback,h5,5,area_upper_limit,Img,Img_hsv,5});
    set(h6_2_text,'Callback',{@myEditBox_Callback,h6,0,len_lower_limit,Img,Img_hsv,6});
    set(h7_2_text,'Callback',{@myEditBox_Callback,h7,3,len_upper_limit,Img,Img_hsv,7});
    set(h8_2_text,'Callback',{@myEditBox_Callback,h8,0,wid_lower_limit,Img,Img_hsv,8});
    set(h9_2_text,'Callback',{@myEditBox_Callback,h9,2,wid_upper_limit,Img,Img_hsv,9});
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
%%%%%%%%%%%%end: if refMarker in each image, update pixel_scale %%%%%%%%%%%%

title('Step2: Drag the sliders to recognize the leaves, press ENTER to proceed');
axh=gca;
set(fig_handle,'pointer','arrow');
set(gcf, 'position', get(0,'ScreenSize'),'MenuBar','none','ToolBar','none');


    area_lower_limit=20; % cm2
    area_lower_default=3;
    area_upper_limit=2000;
    area_upper_default=200;
    len_lower_limit=10; % cm
    len_lower_default=3;
    len_upper_limit=100;
    len_upper_default=100;
    wid_lower_limit=10; % cm
    wid_lower_default=0.2;
    wid_upper_limit=100;
    wid_upper_default=5;

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
h4 = uicontrol('style','slider','units','pixel','position',[80 0.77*screen_height 200 25],'Min',0,'Max',area_lower_limit,'Value',area_lower_default); % cm2
h4_2_text=uicontrol('Style','edit',...
    'Position',[285 0.77*screen_height 45 25],...
    'String',num2str(area_lower_default),'fontsize',10);
h5_text=uicontrol('Style','text',...
    'Position',[5 0.74*screen_height 75 25],...
    'String','Area<','fontsize',10,'backgroundcolor','y');
h5 = uicontrol('style','slider','units','pixel','position',[80 0.74*screen_height 200 25],'Min',0,'Max',area_upper_limit,'Value',area_upper_default); % cm2
h5_2_text=uicontrol('Style','edit',...
    'Position',[285 0.74*screen_height 45 25],...
    'String',num2str(area_upper_default),'fontsize',10);
h6_text=uicontrol('Style','text',...
    'Position',[5 0.71*screen_height 75 25],...
    'String','Len>','fontsize',10,'backgroundcolor','m');
h6 = uicontrol('style','slider','units','pixel','position',[80 0.71*screen_height 200 25],'Min',0,'Max',len_lower_limit,'Value',len_lower_default); % cm
h6_2_text=uicontrol('Style','edit',...
    'Position',[285 0.71*screen_height 45 25],...
    'String',num2str(len_lower_default),'fontsize',10);
h7_text=uicontrol('Style','text',...
    'Position',[5 0.68*screen_height 75 25],...
    'String','Len<','fontsize',10,'backgroundcolor','m');
h7 = uicontrol('style','slider','units','pixel','position',[80 0.68*screen_height 200 25],'Min',0,'Max',len_upper_limit,'Value',len_upper_default); % cm
h7_2_text=uicontrol('Style','edit',...
    'Position',[285 0.68*screen_height 45 25],...
    'String',num2str(len_upper_default),'fontsize',10);
h8_text=uicontrol('Style','text',...
    'Position',[5 0.65*screen_height 75 25],...
    'String','Wid>','fontsize',10,'backgroundcolor','c');
h8 = uicontrol('style','slider','units','pixel','position',[80 0.65*screen_height 200 25],'Min',0,'Max',wid_lower_limit,'Value',wid_lower_default); % cm
h8_2_text=uicontrol('Style','edit',...
    'Position',[285 0.65*screen_height 45 25],...
    'String',num2str(wid_lower_default),'fontsize',10);
h9_text=uicontrol('Style','text',...
    'Position',[5 0.62*screen_height 75 25],...
    'String','Wid<','fontsize',10,'backgroundcolor','c');
h9 = uicontrol('style','slider','units','pixel','position',[80 0.62*screen_height 200 25],'Min',0,'Max',wid_upper_limit,'Value',wid_upper_default); % cm
h9_2_text=uicontrol('Style','edit',...
    'Position',[285 0.62*screen_height 45 25],...
    'String',num2str(wid_upper_default),'fontsize',10);
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

area_low=area_lower_default*pixel_scale^2;
area_high=area_upper_default*pixel_scale^2; % cm2  -> pixels^2
len_low=len_lower_default*pixel_scale;
len_high=len_upper_default*pixel_scale; % cm -> pixels
wid_low=wid_lower_default*pixel_scale;
wid_high=wid_upper_default*pixel_scale; % cm -> pixels

set(h0_2_text,'Callback',{@myEditBox_Callback,h0,0,1,Img,Img_hsv,0});
set(h1_2_text,'Callback',{@myEditBox_Callback,h1,0,1,Img,Img_hsv,1});
set(h2_2_text,'Callback',{@myEditBox_Callback,h2,0,1,Img,Img_hsv,2});
set(h3_2_text,'Callback',{@myEditBox_Callback,h3,0,1,Img,Img_hsv,3});
set(h4_2_text,'Callback',{@myEditBox_Callback,h4,0,area_lower_limit,Img,Img_hsv,4});
set(h5_2_text,'Callback',{@myEditBox_Callback,h5,5,area_upper_limit,Img,Img_hsv,5});
set(h6_2_text,'Callback',{@myEditBox_Callback,h6,0,len_lower_limit,Img,Img_hsv,6});
set(h7_2_text,'Callback',{@myEditBox_Callback,h7,3,len_upper_limit,Img,Img_hsv,7});
set(h8_2_text,'Callback',{@myEditBox_Callback,h8,0,wid_lower_limit,Img,Img_hsv,8});
set(h9_2_text,'Callback',{@myEditBox_Callback,h9,2,wid_upper_limit,Img,Img_hsv,9});
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
if ~isempty(ImgCropMat)
    fprintf(fid, '### Image cropping coordinates(y1,y2;x1,x2): %.2f,%.2f;%.2f,%.2f\n',ImgCropMat);
end
figList=dir(fullfile(fileDir,'*.jpg'));
%% global parameters
RGB=1;
for pic = 1:1:length(figList)
    clearvars -except fileDir figList pic RGB save_fig ImgCropMat fid fid2 hue_low hue_high sat_low area_low area_high medfilt_P pixel_scale len_low len_high wid_low wid_high bw_threshold erode_depth val_low val_high area_low_ref area_high_ref len_low_ref len_high_ref wid_low_ref wid_high_ref hue_low_ref hue_high_ref sat_low_ref bw_threshold_ref erode_depth_ref medfilt_P_ref val_low_ref val_high_ref scale_area Marker_status brightness imageRotationAngle
    fn=figList(pic).name;
    figNA=fn(1:end-4);
    fprintf('%s, the %d/%d-th image in processing...\n',fn,pic,length(figList))
    fprintf(fid2,'%s\n',fn);
    fprintf(fid,'%s\n',fn);
    if abs(brightness-1)>0.01
        I_Raw=imread(fullfile(fileDir,fn))*brightness;
    else
        I_Raw=imread(fullfile(fileDir,fn));
    end
    if abs(imageRotationAngle-90)<0.01
        I_Raw=rot90(I_Raw,1);
    elseif abs(imageRotationAngle-180)<0.01
        I_Raw=rot90(I_Raw,2);
    elseif abs(imageRotationAngle-270)<0.01
        I_Raw=rot90(I_Raw,3);
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
    width_cal_region=0.15*pixel_scale;
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
    
    leaf_order=0;
    for i = 1:num
        LeafArea=stats(i).Area;
        LeafLen_raw=stats(i).MajorAxisLength;
        LeafWid_raw=stats(i).MinorAxisLength;
        if LeafLen_raw>=len_low && LeafWid_raw>=wid_low && LeafArea>=area_low && LeafLen_raw<=len_high && LeafWid_raw<=wid_high && LeafArea<=area_high
            leaf_order=leaf_order+1;
        else
            continue
        end
        bound=I_B{i};
        x=bound(:,2);
        y=bound(:,1);
        [rectx,recty,useless,useless2] = minboundrect(x,y,'a'); % anticlosewise direction
        LeafLength_temp=norm([rectx(3),recty(3)]-[rectx(2),recty(2)])/pixel_scale;
        LeafWidth_temp=norm([rectx(1),recty(1)]-[rectx(2),recty(2)])/pixel_scale;
        long_axis=0;
        if LeafLength_temp>=LeafWidth_temp
            LeafLength=LeafLength_temp;
            LeafWidth=LeafWidth_temp;
            long_axis=23;
        else
            LeafLength=LeafWidth_temp;
            LeafWidth=LeafLength_temp;
            long_axis=12;
        end
        plot(x,y,'-b','LineWidth',2)
        temp = stats(i).Centroid;
        text(temp(1), temp(2), num2str(leaf_order), 'color', 'r','FontSize',20);
        fprintf(fid2, '\t%g\t%g\t',LeafArea/pixel_scale^2,LeafLength);
        sort_y=sort(recty);
        if recty(1)==sort_y(1) || recty(1)==sort_y(2)
            if long_axis==12
                plot(rectx(1),recty(1),'.r','MarkerSize',30);
                wid_list=zeros(1,7);
                y_vector=([rectx(2),recty(2)]-[rectx(1),recty(1)])/norm([rectx(2),recty(2)]-[rectx(1),recty(1)]);
                x_vector=([rectx(4),recty(4)]-[rectx(1),recty(1)])/norm([rectx(4),recty(4)]-[rectx(1),recty(1)]);
            else
                plot(rectx(2),recty(2),'.r','MarkerSize',30);
                wid_list=zeros(1,7);
                y_vector=([rectx(3),recty(3)]-[rectx(2),recty(2)])/norm([rectx(3),recty(3)]-[rectx(2),recty(2)]);
                x_vector=([rectx(1),recty(1)]-[rectx(2),recty(2)])/norm([rectx(1),recty(1)]-[rectx(2),recty(2)]);
            end
        else
            if long_axis==12
                plot(rectx(3),recty(3),'.r','MarkerSize',30);
                wid_list=zeros(1,7);
                y_vector=([rectx(4),recty(4)]-[rectx(3),recty(3)])/norm([rectx(4),recty(4)]-[rectx(3),recty(3)]);
                x_vector=([rectx(2),recty(2)]-[rectx(3),recty(3)])/norm([rectx(2),recty(2)]-[rectx(3),recty(3)]);
            else
                plot(rectx(4),recty(4),'.r','MarkerSize',30);
                wid_list=zeros(1,7);
                y_vector=([rectx(1),recty(1)]-[rectx(4),recty(4)])/norm([rectx(1),recty(1)]-[rectx(4),recty(4)]);
                x_vector=([rectx(3),recty(3)]-[rectx(4),recty(4)])/norm([rectx(3),recty(3)]-[rectx(4),recty(4)]);
            end
        end
        %row_num=mod(leafIdxList{i}-1,Img_row)+1;
        %col_num=floor(leafIdxList{i}-1,Img_row)+1;
        xy_new=[x_vector;y_vector]*[x';y'];
        x_new=xy_new(1,:);
        y_new=xy_new(2,:);
        [y_new_sort,ind]=sort(y_new);
        xy_new=[x_new(ind);y_new(ind)];
        ymin=y_new_sort(1);
        ymax=y_new_sort(end);
        y_pos=linspace(ymin+width_cal_region,ymax-width_cal_region,7);
        %         [useless,ind_temp]=find(xy_new(2,:)<ymin+width_cal_region);
        raw_pos_record=[];
        for pos=1:7
            [useless,ind_temp]=find(xy_new(2,:)<y_pos(pos)+width_cal_region & xy_new(2,:)>y_pos(pos)-width_cal_region);
            xy_temp=sortrows(xy_new(:,ind_temp)')';
            max_ind=floor(length(xy_temp(1,:))*0.9);
            min_ind=floor(length(xy_temp(1,:))*0.1);
            wid_list(pos)=xy_temp(1,max_ind)-xy_temp(1,min_ind);
            raw_pos_record=[raw_pos_record,xy_temp(:,[min_ind,max_ind])];
        end
        fprintf(fid2, '%g\t',wid_list/pixel_scale);
        fprintf(fid2, '%g',max(wid_list)/pixel_scale);
        fprintf(fid2, '\n');
    end
    
    if save_fig
        print(fh,fullfile(fileDir,[figNA,'-track']),'-dpng','-r100');
        close;
    end
    
end

close all;
fclose all;
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
    area_low = area_low*pixel_scale^2; %pixel^2
elseif bar_num==5
    area_high = value_update;
    area_high = area_high*pixel_scale^2; %pixel^2
elseif bar_num==6
    len_low = value_update;
    len_low = len_low*pixel_scale; %pixel
elseif bar_num==7
    len_high = value_update;
    len_high = len_high*pixel_scale; %pixel
elseif bar_num==8
    wid_low = value_update;
    wid_low = wid_low*pixel_scale; %pixel
elseif bar_num==9
    wid_high = value_update;
    wid_high = wid_high*pixel_scale; %pixel
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
[area_low area_high]
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
        bound=B{i};edge_plot=[edge_plot,plot(bound(:,2),bound(:,1),'b-','LineWidth',1)];
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

function [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
% minboundrect: Compute the minimal bounding rectangle of points in the plane
% usage: [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
%
% arguments: (input)
%  x,y - vectors of points, describing points in the plane as
%        (x,y) pairs. x and y must be the same lengths.
%
%  metric - (OPTIONAL) - single letter character flag which
%        denotes the use of minimal area or perimeter as the
%        metric to be minimized. metric may be either 'a' or 'p',
%        capitalization is ignored. Any other contraction of 'area'
%        or 'perimeter' is also accepted.
%
%        DEFAULT: 'a'    ('area')
%
% arguments: (output)
%  rectx,recty - 5x1 vectors of points that define the minimal
%        bounding rectangle.
%
%  area - (scalar) area of the minimal rect itself.
%
%  perimeter - (scalar) perimeter of the minimal rect as found
%
%
% Note: For those individuals who would prefer the rect with minimum
% perimeter or area, careful testing convinces me that the minimum area
% rect was generally also the minimum perimeter rect on most problems
% (with one class of exceptions). This same testing appeared to verify my
% assumption that the minimum area rect must always contain at least
% one edge of the convex hull. The exception I refer to above is for
% problems when the convex hull is composed of only a few points,
% most likely exactly 3. Here one may see differences between the
% two metrics. My thanks to Roger Stafford for pointing out this
% class of counter-examples.
%
% Thanks are also due to Roger for pointing out a proof that the
% bounding rect must always contain an edge of the convex hull, in
% both the minimal perimeter and area cases.
%
%
% See also: minboundcircle, minboundtri, minboundsphere
%
%
% default for metric
if (nargin<3) || isempty(metric)
    metric = 'a';
elseif ~ischar(metric)
    error 'metric must be a character flag if it is supplied.'
else
    % check for 'a' or 'p'
    metric = lower(metric(:)');
    ind = strmatch(metric,{'area','perimeter'});
    if isempty(ind)
        error 'metric does not match either ''area'' or ''perimeter'''
    end
    
    % just keep the first letter.
    metric = metric(1);
end

% preprocess data
x=x(:);
y=y(:);

% not many error checks to worry about
n = length(x);
if n~=length(y)
    error 'x and y must be the same sizes'
end

% start out with the convex hull of the points to
% reduce the problem dramatically. Note that any
% points in the interior of the convex hull are
% never needed, so we drop them.
if n>3
    %edges = convhull(x,y,{'Qt'});  % 'Pp' will silence the warnings
    edges = convhull(x,y);
    % exclude those points inside the hull as not relevant
    % also sorts the points into their convex hull as a
    % closed polygon
    
    x = x(edges);
    y = y(edges);
    
    % probably fewer points now, unless the points are fully convex
    nedges = length(x) - 1;
elseif n>1
    % n must be 2 or 3
    nedges = n;
    x(end+1) = x(1);
    y(end+1) = y(1);
else
    % n must be 0 or 1
    nedges = n;
end

% now we must find the bounding rectangle of those
% that remain.

% special case small numbers of points. If we trip any
% of these cases, then we are done, so return.
switch nedges
    case 0
        % empty begets empty
        rectx = [];
        recty = [];
        area = [];
        perimeter = [];
        return
    case 1
        % with one point, the rect is simple.
        rectx = repmat(x,1,5);
        recty = repmat(y,1,5);
        area = 0;
        perimeter = 0;
        return
    case 2
        % only two points. also simple.
        rectx = x([1 2 2 1 1]);
        recty = y([1 2 2 1 1]);
        area = 0;
        perimeter = 2*sqrt(diff(x).^2 + diff(y).^2);
        return
end
% 3 or more points.

% will need a 2x2 rotation matrix through an angle theta
Rmat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)];

% get the angle of each edge of the hull polygon.
ind = 1:(length(x)-1);
edgeangles = atan2(y(ind+1) - y(ind),x(ind+1) - x(ind));
% move the angle into the first quadrant.
edgeangles = unique(mod(edgeangles,pi/2));

% now just check each edge of the hull
nang = length(edgeangles);
area = inf;
perimeter = inf;
met = inf;
xy = [x,y];
for i = 1:nang
    % rotate the data through -theta
    rot = Rmat(-edgeangles(i));
    xyr = xy*rot;
    xymin = min(xyr,[],1);
    xymax = max(xyr,[],1);
    
    % The area is simple, as is the perimeter
    A_i = prod(xymax - xymin);
    P_i = 2*sum(xymax-xymin);
    
    if metric=='a'
        M_i = A_i;
    else
        M_i = P_i;
    end
    
    % new metric value for the current interval. Is it better?
    if M_i<met
        % keep this one
        met = M_i;
        area = A_i;
        perimeter = P_i;
        
        rect = [xymin;[xymax(1),xymin(2)];xymax;[xymin(1),xymax(2)];xymin];
        rect = rect*rot';
        rectx = rect(:,1);
        recty = rect(:,2);
    end
end

plot(x,y,'.')
line(rectx,recty,'LineWidth',3,'Color',[1,1,0])
% disp(perimeter)
length1=pdist([rectx(1),recty(1);rectx(2),recty(2)]);
length2=pdist([rectx(1),recty(1);rectx(3),recty(3)]);

% get the final rect

% all done

end % mainline end