function [fileDir,save_fig,brightness,Marker_status,scale_area,writeMode] = preprocessing_3dplant(organType,fileDir,save_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing_3dplant 
% This program is used to preprocess images (file directory path, image scale, brightness setting and etc) for further structural  from a batch of images. 
% Input:
% fileDir: full path to images storage Dir, e.g., 'D:/test'
% save_fig: whether to save the processing marked images (saved in the folder fileDir) for backup check, =1 for saving; =0 for not saving
%
% Output:
%   fileDir:
%   save_fig:
%   brightness:
%   pixel_scale:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global fig_handle;
global ax_handle;
global pixel_scale;
global fid;
global fid2;

    warning off;
    %%%%%%%%%%%%%%%%%%% assign variable fileDir,save_fig %%%%%%%%%%%%%%%%%%%
    switch nargin
        case 3
            1;
        case 2
            save_fig=input('Store the result images (1) or not (0)? (default=1): ','s');
        case 1
            fileDir=input('Please set the directory storing the images for processing: ','s');
            save_fig=input('Store the result images (1) or not (0)? (default=1): ','s');
        otherwise
            error('Too many or too few input parameters are given!');
    end
    if isempty(save_fig)
        save_fig=1;
    elseif ischar(save_fig)
        while ~strcmp(save_fig,'0') && ~strcmp(save_fig,'1')
            save_fig=input('Store the result images (1) or not (0)? (default=1): ','s');
        end
        save_fig=str2num(save_fig);
    end
    while ~isdir(fileDir)
        fileDir=input('Error: invalid path! Please set the correct directory storing the images for processing: ','s');
    end
    %%%%%%%%%%%%%%%%%%% set output write mode %%%%%%%%%%%%%%%%%%%
    fnOut2=fullfile(fileDir,'result_clean.txt');
    fnOut=fullfile(fileDir,'result_raw.txt');
    correctInput=0;
    while correctInput==0
        writeMode=input('Overwrite (w) or append (a) to an existing file (default is "a"): ','s');
        if isempty(writeMode)
            writeMode='a';
        end
        if ~exist(fnOut,'file') && ~exist(fnOut2,'file')
            writeMode='w';
        end
        if strcmp(writeMode,'w') || strcmp(writeMode,'a')
            fid2 = fopen(fnOut2,writeMode);
            fid = fopen(fnOut,writeMode);
            correctInput=1;    
        else
            disp('Please enter w or a, we do not accept other characters!')
        end
        fprintf(fid, '%s\n',strcat('*********',datestr(now,'HH:MM:SS mm/dd/yyyy')));
        fprintf(fid2, '%s\n',strcat('*********',datestr(now,'HH:MM:SS mm/dd/yyyy')));
    end
    %%%%%%%%%%%%%%%%%%% set image brightness adjustment coef. %%%%%%%%%%%%%%%%%%%
    brightness_set=0;
    while ~strcmp(brightness_set,'n') && ~strcmp(brightness_set,'y')
        brightness_set=input('Do you want to re-set the brightness of images (y/n)?: ','s');
    end
    if strcmp(brightness_set,'n')
        brightness=1;
    else
        brightness=input('Type in the brightness coefficient directly, or press ENTER and then choose a typical image to set it: ','s');
        if isempty(brightness)
            [fileName,filepath] = uigetfile(fullfile(fileDir,'*.jpg'));
            Img=imread(fullfile(filepath,fileName));
            brightness=brightness_adjustment(Img);
        else
            brightness=str2double(brightness);
        end
        disp(['Brightness has been adjusted! Brightness coefficient for all the following images will be: ',num2str(brightness)]);
    end
    %%%%%%%%%%%%%%%%%%% set image scale %%%%%%%%%%%%%%%%%%%
    Marker_status='';
    scale_area=-1;
    pixel_scale=0;
    images_same_scale=input('Does every image have the same scale (y) or not (n)?: ','s');
    if strcmp(images_same_scale,'y') || strcmp(images_same_scale,'Y')
        1;
    else
        disp('REMINDER: As the images have different scales, there should be a reference marker/scale with known area/length in each image.');
        Marker_status='y';
        scale_area=-1;
    end
    if strcmp(Marker_status,'y') || strcmp(Marker_status,'Y')
        if strcmp(organType,'grain') || strcmp(organType,'leaf')
            scale_area=input('Please input the reference marker area (cm2): ');
        else
            scale_area=-1;
        end
        pixel_scale=-1;
    else
        pixel_scale=input('Type in the scale (pixels/cm) or press ENTER and then choose an image to set it: ','s');
        if pixel_scale
            while 1
                try
                    pixel_scale=str2double(pixel_scale);
                    break;
                catch
                    pixel_scale=input('Please type in the correct (numerical) scale value (pixels/cm): ');
                end
            end
        else
            [fileName,filepath] = uigetfile(fullfile(fileDir,'*.jpg'));
            if abs(brightness-1)>0.01
                Img=imread(fullfile(filepath,fileName))*brightness;
            else
                Img=imread(fullfile(filepath,fileName));
            end
            fig_handle=figure; imshow(Img);ax_handle=gca;title('Drag the mouse to zoom in if necessary, press ENTER to proceed.');
            set(gcf, 'position', get(0,'ScreenSize'));
            hold on
            zoom on
            k = waitforbuttonpress;
            while k~=1
                k = waitforbuttonpress;
            end
            if ~isdeployed
                input(''); % to feed the waitforbuttonpress
            end
            zoom off;
            figure(1);title('Left-click on two terminal points with known distance, press ENTER to proceed');
            set(fig_handle,'pointer','cross');
            scale_extractor;
            while 1
                input('','s');
                break
            end
        end
    end
    