function TillerAnalyzer_Run(fileDir,save_fig,imageRotationAngle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TillerAnalyzer_Run.m
% This program is aimed at analyzing all leaves (and panicle) length, angle, and position on a rice tiller.
% Usage:
%   TillerAnalyzer_Run(fileDir,save_fig,imageRotationAngle)
%      fileDir: full path to images storage Dir, e.g., 'D:/test'
%      save_fig: whether to save the result images (saved in the fileDir folder) for review, =1 for saving; =0 for not saving
%      imageRotationAngle: Rotation angle (0, 90,180 or 270 degrees) of
%        images for upward-orientation of tillers
% Output:
%   file 'result_raw.txt' in directory fileDir, which records raw coordinates of selected points;
%   file 'result_clean.txt' in directory fileDir, which records value of extracted parameters for each image;
% Debugging record:
% 2018/6/28: calculation of vertical and horizontal length of leaf was corrected: project base2tip vector on stem and anti-stem direction
% 2018/6/29: calculation of leaf base angle was corrected: Angle between
%    nearest stem direction and initial leaf direction (vector of the first two picked points on the leaf)
% 2020/2/7:  code and function improvement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars -except fileDir save_fig imageRotationAngle
close all;
fclose all;
global ImgFile;
global filename;
global fid;
global fid2;
global enlarge_status;
enlarge_status=0;
global xyList2;
global totList;
global organNum;
global plotPointHandle;
global plotPointHandle2;
global plotLineHandle;
global currPt;
global organMarker;
global organMarker2;
global fig_handle;
global ax_handle;
global track_fileDir;
global save_img;
global Marker_status;
global pixel_scale;

xyList2=zeros(0,2);
totList={};
organNum=0;
plotPointHandle={};
plotPointHandle2={};
plotLineHandle={};
currPt=zeros(0,2);

organMarker={'ST','L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16','L17'};
organMarker2='P';

switch nargin
    case 3
        [fileDir,save_img,brightness,Marker_status,~,writeMode] = preprocessing_3dplant('leaf',fileDir,save_fig);
        imageRotationAngle=num2str(imageRotationAngle);
    case 2
        [fileDir,save_img,brightness,Marker_status,~,writeMode] = preprocessing_3dplant('tiller',fileDir,save_fig);
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are vertically oriented and the newest (uppermost) leaf is on the left side (default=0): ','s');
    case 1
        [fileDir,save_img,brightness,Marker_status,~,writeMode] = preprocessing_3dplant('tiller',fileDir);
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are vertically oriented and the newest (uppermost) leaf is on the left side (default=0): ','s');
    case 0
        [fileDir,save_img,brightness,Marker_status,~,writeMode] = preprocessing_3dplant('tiller');
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that leaves are vertically oriented and the newest (uppermost) leaf is on the left side (default=0): ','s');
    otherwise
        error('Too many input parameters are given!');
end

correctInput=0;
while ~correctInput
    if isempty(imageRotationAngle)
        imageRotationAngle='0';
    end
    try
        imageRotationAngle=str2num(imageRotationAngle);
        if abs(imageRotationAngle-0)>0.01 && abs(imageRotationAngle-90)>0.01 && abs(imageRotationAngle-180)>0.01 && abs(imageRotationAngle-270)>0.01
            imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that tillers are upward oriented (default=0): ','s');
        else
            correctInput=1;
        end
    catch
        imageRotationAngle=input('Rotate images for 90,180 or 270 degrees so that tillers are upward oriented (default=0): ','s');
    end
end

if strcmp(writeMode,'w')
    if pixel_scale>0
        fprintf(fid, 'File name\nImage cropping coordinates (y1,y2;x1,x2)\nStem trace coordinates\nLeaves (and panicle) trace coordinates\n');
    else
        fprintf(fid, 'File name\nImage cropping coordinates (y1,y2;x1,x2)\nScale(cm,pixels/cm: scale trace coordinates)\nStem trace coordinates\nLeaves (and panicle) trace coordinates\n');
    end
    fprintf(fid2, 'File name\n\tLeaf base height(cm)\tLeaf length(cm)\tLeaf angle(rad)\tLeaf tip vertical length(cm)\tLeaf tip horizontal length(cm)\n');
end

track_fileDir=fileDir;

figList=dir(fullfile(fileDir,'*.jpg'));
disp(['There are ', num2str(length(figList)) ,' image(s) in total.'])

correctInput=0;
while correctInput==0
    startFrom=input('Start from which image? (default is 1): ');
    if isempty(startFrom)
        startFrom=1;
    end
    if startFrom>length(figList)
        disp(['There are only ', num2str(length(figList)) ,' image(s) in total! Please type in a number <= ',num2str(length(figList)),'!'])
    else
        correctInput=1;
    end
end

total_image_for_process=length(figList)-startFrom+1;
early_stop=0;
for pic = startFrom:1:length(figList)
%     clearvars -except ImgFile filename fid fid2 enlarge_status xyList2 totList organNum plotPointHandle plotPointHandle2 plotLineHandle currPt organMarker organMarker2 fig_handle ax_handle track_fileDir save_img Marker_status pixel_scale total_image_for_process early_stop pic startFrom figList fileDir brightness imageRotationAngle 

    filename=figList(pic).name;
    ImgFile = fullfile(fileDir,filename);
    fprintf(fid, '%s\n',strcat('<<',filename));
    fprintf(fid2, '%s\n',strcat('<<',filename));
    
    xyList2=zeros(0,2);
    totList={};
    organNum=0;
    plotPointHandle={};
    plotPointHandle2={};
    plotLineHandle={};
    currPt=zeros(0,2);
    disp([filename,', the ',num2str(pic-startFrom+1),'/',num2str(total_image_for_process),' image is in processing...']);
    
    if abs(brightness-1)>0.01
        Img=imread(ImgFile)*brightness;
    else
        Img=imread(ImgFile);
    end
    if abs(imageRotationAngle-90)<0.01
        Img=rot90(Img,1);
    elseif abs(imageRotationAngle-180)<0.01
        Img=rot90(Img,2);
    elseif abs(imageRotationAngle-270)<0.01
        Img=rot90(Img,3);
    end

    fig_handle=figure;imshow(Img);axis on; ax_handle=gca; title('Step1/3: Zoom in get a better view if necessary, press ENTER to proceed');
    %     set(ax_handle,'position',[0 0 1 1]); %,'Box','On','LineWidth',2
    set(fig_handle, 'position', get(0,'ScreenSize'));
    %set(fig_handle,'MenuBar','none');
    %set(fig_handle,'ToolBar','none');
    %set(fig_handle,'Pointer','fullcrosshair');
    hold on
    zoom on;
    
    k = waitforbuttonpress;
    while k~=1
        k = waitforbuttonpress;
    end
    if ~isdeployed
        input(''); % to feed the waitforbuttonpress
    end
    zoom off;
    ImgCropMat=round(axis);
    ImgCropMat(2)=ImgCropMat(2)-1;
    ImgCropMat(4)=ImgCropMat(4)-1;
    fprintf(fid, '%d,%d;%d,%d\tImage_Crop_Coordinates\n',[ImgCropMat(3),ImgCropMat(4),ImgCropMat(1),ImgCropMat(2)]);
    figure(1);title('Step2/3: Click to assign two points with known distance, press ENTER to proceed');
    set(fig_handle,'pointer','cross')
    tiller_analyzer();
    while 1
        char_now=input('','s');
        %char_now=get(gcf, 'CurrentCharacter');
        if strcmp(char_now,'')
            break
        else
            early_stop=1;
            break
        end
    end
    if ~isempty(totList)
        stem_pts=totList{1};
        stem_pts_num=size(stem_pts,1);
        stem_segments_pts=[];
        N=200;
        for j=1:stem_pts_num-1
            for k=1:N
                stem_segments_pts(end+1,:)=stem_pts(j,1:2)+(k-1)*(stem_pts(j+1,1:2)-stem_pts(j,1:2))/N;
            end
        end
        stem_segments_num=N*(stem_pts_num-1);
        %         stem_base_pt=stem_pts(1,1:2);
        %stem_length=sum(sqrt(sum((stem_pts(2:end,:)-stem_pts(1:end-1,:)).^2,2)));
        %stem_direction=(stem_pts(end,:)-stem_pts(1,:))/norm(stem_pts(end,:)-stem_pts(1,:));
        %fprintf(fid2,'%g\t%g\n', stem_base_pt(1),stem_base_pt(2));
        %         last_leafBase_pts=stem_base_pt;
        %         leaf_base_height=0;
        for i=2:1:length(totList)
            leaf_pts=totList{i};
            leaf_length=sum(sqrt(sum((leaf_pts(2:end,1:2)-leaf_pts(1:end-1,1:2)).^2,2)));
            leaf_base_pt=leaf_pts(1,1:2);
            vec_leaf=leaf_pts(2,1:2)-leaf_base_pt;
            dist_list=sqrt(sum((repmat(leaf_base_pt,stem_segments_num,1)-stem_segments_pts).^2,2));
            [min_dist,min_ind]=min(dist_list);
            vec_stem=stem_segments_pts(min(min_ind+1,stem_segments_num),:)-stem_segments_pts(max(min_ind-1,1),:);
            len_vec_stem=norm(vec_stem);
            leaf_angle=acos(dot(vec_stem,vec_leaf)/(len_vec_stem*norm(vec_leaf)));
            leaf_base_height=0;
            if min_ind>1
                leaf_base_height=sum(sqrt(sum((stem_segments_pts(2:min_ind,1:2)-stem_segments_pts(1:min_ind-1,1:2)).^2,2)));
            end
            leaf_base_tip=leaf_pts(end,1:2)-leaf_base_pt;
            z_direction=vec_stem/len_vec_stem;
            y_direction=[z_direction(2),-z_direction(1)];
            leaf_tipH=dot(leaf_base_tip,z_direction); % vertical height, could be negative
            leaf_tipR=abs(dot(leaf_base_tip,y_direction)); % horizontal length, positive
            if length(leaf_pts(1,:))<3
                fprintf(fid2, '\t%0.2f\t%0.2f\t%0.3f\t%0.2f\t%0.2f\tleaf\n', leaf_base_height/pixel_scale, leaf_length/pixel_scale, leaf_angle, leaf_tipH/pixel_scale, leaf_tipR/pixel_scale);
            else
                fprintf(fid2, '\t%0.2f\t%0.2f\t%0.3f\t%0.2f\t%0.2f\tpanicle\n', leaf_base_height/pixel_scale, leaf_length/pixel_scale, leaf_angle, leaf_tipH/pixel_scale, leaf_tipR/pixel_scale);
            end
            %             last_leafBase_pts=leaf_base_pt;
        end
    end
    if early_stop==1
        break
    end
    if Marker_status=='y'
        pixel_scale=-1;
    end
end

close all;
fclose all;
warning on;
end