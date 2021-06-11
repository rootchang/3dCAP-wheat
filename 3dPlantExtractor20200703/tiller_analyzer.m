function tiller_analyzer(action)
global track_fileDir;
global filename;
global fid;
global fig_handle;
global ax_handle;
global xyList2;
global totList;
global organNum;
global plotPointHandle;
global plotPointHandle2;
global plotLineHandle;
global currPt;
global enlarge_status;
global organMarker;
global organMarker2;
global save_img;
global Marker_status;
global pixel_scale;

if nargin == 0, action = 'start'; end

switch(action)
    case 'start'
        set(fig_handle, 'WindowButtonMotionFcn', 'tiller_analyzer move1');
    case 'move1'
        H = get(fig_handle,'UserData');
        if ~isempty(H)
            f1 = H(1); a1 = H(2); a2 = H(3);
            a2_param = get(a2,'UserData');
            f_pos = get(f1,'Position');
            a1_pos = get(a1,'Position');
            
            [f_cp, a1_cp] = pointer2d(f1,a1);
            
            set(a2,'Position',[(f_cp./f_pos(3:4)) 0 0]+a2_param(2)*a1_pos(3)*[-1 -1 2 2]);
            a2_pos = get(a2,'Position');
            axis square
            
            set(a2,'XLim',a1_cp(1)+(1/a2_param(1))*(a2_pos(3)/a1_pos(3))*diff(get(a1,'XLim'))*[-0.5 0.5]);
            set(a2,'YLim',a1_cp(2)+(1/a2_param(1))*(a2_pos(4)/a1_pos(4))*diff(get(a1,'YLim'))*[-0.5 0.5]);
            %set(get(a2,'Children'),'Box','on','Color','k','LineWidth', 1.2,'Layer','top');
        end;
        set(fig_handle, 'KeyPressFcn', 'tiller_analyzer char');
        set(fig_handle, 'WindowButtonDownFcn', 'tiller_analyzer down');
        
    case 'move2'
        currPt = get(ax_handle, 'CurrentPoint');
        if ~isempty(xyList2)
            if ~isempty(plotLineHandle)
                set(plotLineHandle{end}, 'XData', [currPt(1,1),xyList2(end,1)], 'YData', [currPt(1,2),xyList2(end,2)])
            else
                plotLineHandle{end+1}=plot([currPt(1,1),xyList2(end,1)],[currPt(1,2),xyList2(end,2)],'w--','lineWidth',1);
            end
        end
    case 'down'
        currPt = get(ax_handle, 'CurrentPoint');
        x = currPt(1,1);
        y = currPt(1,2);
        xyList2 = [xyList2; [x,y]];
        plotPointHandle{end+1}=scatter(x, y,20,'ro','filled');
        plotPointHandle2{end+1}=copyobj(plotPointHandle{end}, ax_handle); % copy chosen point from zoom-in panel to original figure
        if enlarge_status==1
            set(fig_handle, 'WindowButtonMotionFcn', 'tiller_analyzer move1');
        else
            set(fig_handle, 'WindowButtonMotionFcn', 'tiller_analyzer move2');
        end
    case 'char'
        char_now=double(get(fig_handle, 'CurrentCharacter'));
        switch(char_now)
            case 'z'
                if ~isempty(plotPointHandle)
                    delete(plotPointHandle{end});
                    delete(plotPointHandle2{end});
                    delete(plotLineHandle{end});
                    plotPointHandle(end)=[];
                    plotPointHandle2(end)=[];
                    plotLineHandle(end)=[];
                    xyList2 = xyList2(1:end-1,:);
                end
            case 9 % Tab  %% for end of a leaf
                if size(xyList2,1)>1 % at least two points
                    organNum=organNum+1;
                    
                    if organNum==1
                        for i=1:length(xyList2(:,1))-1
                            plot([xyList2(i,1),xyList2(i+1,1)],[xyList2(i,2),xyList2(i+1,2)],'k-.','lineWidth',0.1);
                            delete(plotPointHandle{end});
                            delete(plotPointHandle2{end});
                            plotPointHandle(end)=[];
                            plotPointHandle2(end)=[];
                        end
                        aaa=text(xyList2(1,1), xyList2(1,2),organMarker{organNum},'FontSize',10,'FontWeight','bold','color','b');
                    else
                        aaa=text(xyList2(end,1), xyList2(end,2),organMarker{organNum},'FontSize',10,'FontWeight','bold','color','b');
                    end
                    for i=1:1:length(xyList2(:,1))
                        fprintf(fid,'%g,%g\t', xyList2(i,1),xyList2(i,2));
                    end
                    if organNum==1
                        fprintf(fid,'Stem\n');
                    else
                        fprintf(fid,'Leaf%d\n',organNum-1);
                    end
                    totList{organNum}=xyList2;
                    xyList2 = zeros(0,2);
                    copyobj(aaa, ax_handle);
                    delete(aaa);
                    delete(plotLineHandle{end});
                    plotPointHandle={};
                    plotPointHandle2={};
                    plotLineHandle={};
                else
                    msgbox('Please assign at least one more point before press TAB!');
                end
            case {113,49} % 'q' or '1', end of a panicle or last leaf
                if size(xyList2,1)>1 % at least two points
                    organNum=organNum+1;
                    for i=1:1:length(xyList2(:,1))
                        fprintf(fid,'%g,%g\t', xyList2(i,1),xyList2(i,2));
                    end
                    if char_now==49 % 1  %% for end of a panicle
                        aaa=text(xyList2(end,1), xyList2(end,2),organMarker2,'FontSize',10,'FontWeight','bold','color','b');
                        fprintf(fid,'Panicle\n');
                        xyList2(1,3)=0; % add an additional column
                    else
                        if organNum==1
                            for i=1:length(xyList2(:,1))-1
                                plot([xyList2(i,1),xyList2(i+1,1)],[xyList2(i,2),xyList2(i+1,2)],'k-.','lineWidth',0.1);
                                delete(plotPointHandle{end});
                                delete(plotPointHandle2{end});
                                plotPointHandle(end)=[];
                                plotPointHandle2(end)=[];
                            end
                            aaa=text(xyList2(1,1), xyList2(1,2),organMarker{organNum},'FontSize',10,'FontWeight','bold','color','b');
                            fprintf(fid,'Stem');
                        else
                            aaa=text(xyList2(end,1), xyList2(end,2),organMarker{organNum},'FontSize',10,'FontWeight','bold','color','b');
                            fprintf(fid,['Leaf',num2str(organNum-1)]);
                        end
                        fprintf(fid,'\n');
                    end
                    totList{organNum}=xyList2;
                    xyList2=zeros(0,2);
                    copyobj(aaa, ax_handle);
                    delete(aaa);
                    delete(plotLineHandle{end});
                    if save_img
                        img = getframe(gcf);
                        imwrite(img.cdata, fullfile(track_fileDir,[filename(1:end-4),'-track.png']));
                    end
                    if char_now=='1'
                        pause(1);
                    end
                    close;
                    disp('Continue to process the next figure, press "Enter"; stop here, press "s" and "Enter":')
                    return
                elseif size(xyList2,1)==1
                    msgbox('Please assign at least one more point before press q or 1!');
                else
                    if save_img
                        img = getframe(gcf);
                        imwrite(img.cdata, fullfile(track_fileDir,[filename(1:end-4),'-track.png']));
                    end
                    close;
                    disp('Continue to process the next figure, press "Enter"; stop here, press "s" and "Enter":')
                    return
                end
            case 'e' % 101
                if enlarge_status==0
                    enlarge_status=1;
                    f1 = fig_handle;
                    a1 = get(f1,'CurrentAxes');
                    a2 = copyobj(a1,f1);

                    set(f1,'UserData',[f1,a1,a2],'CurrentAxes',a2); %'Pointer','crosshair',
                    %iptPointerManager(a2, 'disable');
                    set(a2,'UserData',[5,0.15], ...  %magnification, frame size
                        'XTick',[],'YTick',[]); %'Box','on','LineWidth', 2,'Layer','top'
                    set(fig_handle, 'WindowButtonMotionFcn', 'tiller_analyzer move1');
                else
                    enlarge_status=0;
                    H = get(fig_handle,'UserData');
                    f1 = H(1); a1 = H(2); a2 = H(3);
                    set(a1,'Color',get(a2,'Color'));
                    set(f1,'UserData',[],'CurrentAxes',a1);
                    if ~strcmp(get(f1,'SelectionType'),'alt')
                        delete(a2);
                    end;
                    set(fig_handle, 'WindowButtonMotionFcn', 'tiller_analyzer move2');
                end
            case 13 % Enter
                if (strcmp(Marker_status,'y') || strcmp(Marker_status,'Y')) && pixel_scale<0
                    if size(xyList2,1)>1 % at least two points
                        set(fig_handle,'visible','off');
                        scale_length=input('The true length of the scale is (cm): ');
                        set(fig_handle,'visible','on');
                        title('Step3/3: Left-click to assign key points of organs on tiller, end up with pressing TAB (stem and leaves) or 1 (spike), press q after finish.');
                        for i=1:length(xyList2(:,1))-1
                            plot([xyList2(i,1),xyList2(i+1,1)],[xyList2(i,2),xyList2(i+1,2)],'r-','lineWidth',1);
                            delete(plotPointHandle{end});
                            delete(plotPointHandle2{end});
                            plotPointHandle(end)=[];
                            plotPointHandle2(end)=[];
                        end
                        aaa=text(xyList2(1,1), xyList2(1,2),['Scale (cm): ',num2str(scale_length)],'FontSize',10,'FontWeight','bold','color','r');
                        xyList_copy=xyList2;
                        xyList2 = zeros(0,2);
                        copyobj(aaa, ax_handle);
                        delete(aaa);
                        delete(plotLineHandle{end});
                        scale_in_pixel=sum(sqrt(sum((xyList_copy(2:end,1:2)-xyList_copy(1:end-1,1:2)).^2,2)));
                        pixel_scale=scale_in_pixel/scale_length;
                        fprintf(fid,'%g,%g: ', [scale_length,pixel_scale]);
                        for i=1:1:length(xyList_copy(:,1))
                            fprintf(fid,'%g,%g\t', xyList_copy(i,1),xyList_copy(i,2));
                        end
                        fprintf(fid,'Scale\n');
                        plotPointHandle={};
                        plotPointHandle2={};
                        plotLineHandle={};
                    else
                        msgbox('Please assign at least one more point before press ENTER!');
                    end
                end
        end
end
end

function [fig_pointer_pos, axes_pointer_val] = pointer2d(fig_hndl,axes_hndl)
if (nargin == 0), fig_hndl = gcf; axes_hndl = gca; end;
if (nargin == 1), axes_hndl = get(fig_hndl,'CurrentAxes'); end;

set(fig_hndl,'Units','pixels');

pointer_pos = get(0,'PointerLocation');	%pixels {0,0} lower left
fig_pos = get(fig_hndl,'Position');	%pixels {l,b,w,h}

fig_pointer_pos = pointer_pos - fig_pos([1,2]);
set(fig_hndl,'CurrentPoint',fig_pointer_pos);

if (isempty(axes_hndl))
	axes_pointer_val = [];
elseif (nargout == 2)
	axes_pointer_line = get(axes_hndl,'CurrentPoint');
	axes_pointer_val = sum(axes_pointer_line)/2;
end;
end