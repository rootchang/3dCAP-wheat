function brightness_new=brightness_adjustment(Img1)
    global brightness;

%%%%%%%%%%%%%%%%%%%%%%%%% BrightAdjust %%%%%%%%%%%%%%%%%%%%%%%%%
    screen_size_4=get(0,'ScreenSize');
    screen_height=screen_size_4(4);
    brightness=1;
    fig_handle=figure;imshow(Img1);title('Drag the slider to adjust the brightness, press ENTER to continue');
    hold on;
    set(fig_handle,'pointer','arrow');set(fig_handle, 'position', get(0,'ScreenSize'));axh=gca;
    
    h1_text=uicontrol('Style','text',...
        'Position',[5 0.85*screen_height 45 20],...
        'String',num2str(brightness),'fontsize',16,'backgroundcolor','g');
    h1 = uicontrol('style','slider','units','pixel','position',[50 0.85*screen_height 300 20],'Min',0,'Max',3,'Value',brightness);
    try
        addlistener(h1,'ActionEvent',@(hObject, ContinuousValueChange) BrightAdjust(hObject, ContinuousValueChange,Img1,h1_text,axh));
    catch
        addlistener(h1,'ContinuousValueChange',@(hObject, ContinuousValueChange) BrightAdjust(hObject, ContinuousValueChange,Img1,h1_text,axh));
    end
    k = 0;
    while k~=1
        k = waitforbuttonpress;
    end
    close;
    brightness_new=brightness;
    clear brightness;
end

function BrightAdjust(hObject,ContinuousValueChange,I,h_text,axh)
    global brightness;

    brightness = get(hObject,'Value');
    I=I*brightness;
    cla(axh); % to avoid keeping increase of memory use
    imshow(I);
    %drawnow;
    set(h_text,'String',num2str(round(brightness*100)/100));
end