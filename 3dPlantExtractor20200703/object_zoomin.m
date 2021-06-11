function object_zoomin(action)
    global fig_handle;
    global enlarge_status;
    global hue_low;
    global hue_high;
    global sat_low;
    global medfilt_P;
    global area_low;
    global area_high;
    global len_low;
    global len_high;
    global wid_low;
    global wid_high;
    global pixel_scale;
    global erode_depth;
    global bw_threshold;
    
    if nargin == 0, action = 'start'; end

    switch(action)
        case 'start',
            set(fig_handle, 'WindowButtonMotionFcn', 'object_zoomin move');
        case 'move',
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
            set(fig_handle, 'KeyPressFcn', 'object_zoomin char');
    case 'char',
            char_now=double(get(fig_handle, 'CurrentCharacter'));
            switch(char_now)
                case 'e', % 101
                    if enlarge_status==0
                        enlarge_status=1;
                        f1 = fig_handle;
                        a1 = get(f1,'CurrentAxes');
                        a2 = copyobj(a1,f1);

                        set(f1,'UserData',[f1,a1,a2],'CurrentAxes',a2); %'Pointer','crosshair',
                        %iptPointerManager(a2, 'disable');
                        set(a2,'UserData',[5,0.15], ...  %magnification, frame size
                            'XTick',[],'YTick',[]); %'Box','on','LineWidth', 2,'Layer','top'
    %                     set(fig_handle, 'WindowButtonMotionFcn', 'object_zoomin move');
                    else
                        enlarge_status=0;
                        H = get(fig_handle,'UserData');
                        f1 = H(1); a1 = H(2); a2 = H(3);
                        set(a1,'Color',get(a2,'Color'));
                        set(f1,'UserData',[],'CurrentAxes',a1);
                        if ~strcmp(get(f1,'SelectionType'),'alt')
                            delete(a2);
                        end;
    %                     set(fig_handle, 'WindowButtonMotionFcn', 'object_zoomin move');
                    end
                case 13, % ENTER
                    close;
                    fprintf('The parameters for recognition have been set! They are hue_lower=%.2f, hue_upper=%.2f, sat_lower=%.2f, area_lower=%.1fcm2, area_upper=%.1fcm2,len_lower=%.1fcm,len_upper=%.1fcm,wid_lower=%.1fcm,wid_upper=%.1fcm,imFilter=%d,imErode=%d,bw_threshold=%.2f\n',...
                        hue_low,hue_high,sat_low,area_low/pixel_scale^2,area_high/pixel_scale^2,len_low/pixel_scale,len_high/pixel_scale,wid_low/pixel_scale,wid_high/pixel_scale,medfilt_P,erode_depth,bw_threshold);
                    disp('To apply these parameters to all the following images, press "Enter"; stop here, press "s" and "Enter":');
                    return;
            end
    end
end

function [fig_pointer_pos, axes_pointer_val] = pointer2d(fig_hndl,axes_hndl)
    %
    %pointer2d(fig_hndl,axes_hndl)
    %
    %	Returns the coordinates of the pointer (in pixels)
    %	in the desired figure (fig_hndl) and the coordinates
    %       in the desired axis (axes coordinates)
    %
    % Example:
    %  figure(1),
    %  hold on,
    %  for i = 1:1000,
    %     [figp,axp]=pointer2d;
    %     plot(axp(1),axp(2),'.','EraseMode','none');
    %     drawnow;
    %  end;
    %  hold off

    % Rick Hindman - 4/18/01

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
