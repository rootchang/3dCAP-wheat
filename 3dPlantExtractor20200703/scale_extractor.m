function scale_extractor(action)
global fig_handle;
global ax_handle;
global pixel_scale;
global enlarge_status;
global xyList2;
global plotPointHandle;
global plotPointHandle2;
global plotLineHandle;
global currPt;
global fid;

if nargin == 0, action = 'start'; end
    
switch(action)
    case 'start'
        set(fig_handle, 'WindowButtonMotionFcn', 'scale_extractor move1');
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
        set(fig_handle, 'KeyPressFcn', 'scale_extractor char');
        set(fig_handle, 'WindowButtonDownFcn', 'scale_extractor down');
        
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
            set(fig_handle, 'WindowButtonMotionFcn', 'scale_extractor move1');
        else
            set(fig_handle, 'WindowButtonMotionFcn', 'scale_extractor move2');
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
                    set(fig_handle, 'WindowButtonMotionFcn', 'scale_extractor move1');
                else
                    enlarge_status=0;
                    H = get(fig_handle,'UserData');
                    f1 = H(1); a1 = H(2); a2 = H(3);
                    set(a1,'Color',get(a2,'Color'));
                    set(f1,'UserData',[],'CurrentAxes',a1);
                    if ~strcmp(get(f1,'SelectionType'),'alt')
                        delete(a2);
                    end;
                    set(fig_handle, 'WindowButtonMotionFcn', 'scale_extractor move2');
                end
            case 13 % Enter
                    if size(xyList2,1)>1 % at least two points
                        for i=1:length(xyList2(:,1))-1
                            plot([xyList2(i,1),xyList2(i+1,1)],[xyList2(i,2),xyList2(i+1,2)],'r-','lineWidth',1);
                            delete(plotPointHandle{end});
                            delete(plotPointHandle2{end});
                            plotPointHandle(end)=[];
                            plotPointHandle2(end)=[];
                        end
                        set(fig_handle,'visible','off');
                        scale_length=input('The true length of the scale is (cm): ');
                        text(xyList2(1,1), xyList2(1,2),['Scale (cm): ',num2str(scale_length)],'FontSize',10,'FontWeight','bold','color','r');
                        scale_in_pixel=sum(sqrt(sum((xyList2(2:end,1:2)-xyList2(1:end-1,1:2)).^2,2)));
                        pixel_scale=scale_in_pixel/scale_length;
                        fprintf(fid,'Scale(cm,pixels/cm: scale trace coordinates)\n');
                        fprintf(fid,'%g,%g: ', [scale_length,pixel_scale]);
                        for i=1:1:length(xyList2(:,1))
                            fprintf(fid,'%g,%g\t', xyList2(i,1),xyList2(i,2));
                        end
                        fprintf(fid,'Scale\n');
                        
                        close;
                        disp(['Scale has been set! Scale for the following images is: ',num2str(pixel_scale),' pixels/cm. Press ENTER to proceed...']);
                        return
                    else
                        msgbox('Please click on at least one more point before press ENTER!');
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