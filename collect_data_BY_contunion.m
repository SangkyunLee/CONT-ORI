function collect_data_BY_contunion(seslist, path_list,sesinfo, fndata, fnsave, param4cellselect, sepdata)
% function collect_data_BY_contunion(seslist, path_list, fnsave, param4cellselect)

% pre_stimtime_ms = param4cellselect.pre_stimtime_ms;
selcontrast = param4cellselect.selcontrast;
param4cellselect.animal_state  = sesinfo.animal_state;
% visual_resp_thr = param4cellselect.visual_resp_thr;

evtsort = param4cellselect.evtsort;

for inx_session = seslist

    info.session_path = path_list{inx_session};
    info.animal_state = sesinfo.animal_state;
    if iscell(fndata) && length(fndata)==length(seslist)
        datafullfn=fullfile(path_list{inx_session},'/matlab/data',fndata{inx_session});
    else
        datafullfn=fullfile(path_list{inx_session},'/matlab/data',fndata);
    end

    % generate A OBJECT    
    CO = cont_ori_session(datafullfn,info);            
    
    % -------- set eyepar------------
    if isfield(param4cellselect,'eyethr') && ...
            ~isempty(param4cellselect.eyethr)
        
        CO = TP.get_eyeparfn(CO);
        eyepar = param4cellselect.eyethr;
        CO = CO.set_eyepar(eyepar);
        param4cellselect.beye = true;
        eye = CO.eye;
    else
        param4cellselect.beye = false;
    end
    
    
    %------ calculate_motion_rotaroadNimageframe   
    if isfield(param4cellselect,'motionthr')   && ...
            ~isempty(param4cellselect.motionthr)     
        
        motionthr = param4cellselect.motionthr;
        twin = motionthr.twin;
        channel = motionthr.channel;
        CO = CO.set_motion_rotaroad_imageframe(twin,channel);
        
%     else
%         motionthr = struct([]);
    end
    %---------------------------------------------
    
    
    
    if nargin>6,
        set_subdata = sepdata(inx_session,:);
    else
        set_subdata{1} = CO.scanId;
    end
    
    
    
    
    sIdmap = zeros(1, max(CO.scanId));    
    sIdmap(CO.scanId)=1:CO.no_scan;
    for isub = 1 :length(set_subdata)
        scanlist0 = set_subdata{isub}; 
        if isempty(scanlist0)
            continue;
        end
        scanlist = sIdmap(scanlist0);
        iscan = scanlist(1);
        stimparam = CO.scans(iscan).Params.stimparam;
        sf_NI = CO.scans(iscan).Params.samplingfreq_NI;
        stimtime_ms = stimparam.stim_samplesinNI/sf_NI*1000;

        



        %--------- collect data
        param4cellselect.stimtime_ms = stimtime_ms;
        out = collect_data(CO, scanlist,param4cellselect);
        
        
        %------- remove any scans with no events selected
        valscan = ~cellfun(@isempty,out.trialinxs);
        if sum(valscan) ~=length(scanlist)
            msg = sprintf('Session%d: from scanID=%s, scan %s has no trials selected',inx_session,num2str(CO.scanId),num2str(CO.scanId(~valscan)));
            h = msgbox(msg,'TRIAL SELECTION');
            uiwait(h);
            scanlist = scanlist(valscan);
            out = collect_data(CO, scanlist,param4cellselect);            
        end
        flds = fieldnames(out);
        for i = 1 : length(flds)
            cmdstr = sprintf('%s = out.%s;',flds{i},flds{i});
            eval(cmdstr);
        end


        %------- collect stimulus conditions
        trialcond= CO.collect_trial_conditions(trialinxs, scanlist, evtsort);


        %---------------
        cellinx_sel = inxcell_peakresptime_conts{1};
        for icont =2: length(selcontrast)
            cellinx_sel = union(cellinx_sel, inxcell_peakresptime_conts{icont});
        end

        %---------- DATA concatenation


        % gather spk data
        Xdata=MXspk(:,1);
        [Y0, events_cont_inx, events_ori_inx] = cont_ori_session.concatenate_trialMresp(Xdata,trialcond,1:length(trialcond));
        Xdata=MXspk(:,2);
        [Y] = cont_ori_session.concatenate_trialMresp(Xdata,trialcond,1:length(trialcond));
        Xsel = Y(:,:,cellinx_sel);
        rXsel = Y(:,:,cellinx_sel) - Y0(:,:,cellinx_sel);
        Xsel = reshape(Xsel,[size(Xsel,1)*size(Xsel,2), size(Xsel,3)]);
        rXsel = reshape(rXsel,[size(rXsel,1)*size(rXsel,2), size(rXsel,3)]);

        % gather dff data
        Xdata=MXc(:,1);
        [Y0] = cont_ori_session.concatenate_trialMresp(Xdata,trialcond,1:length(trialcond));
        Xdata=MXc(:,2);
        [Y] = cont_ori_session.concatenate_trialMresp(Xdata,trialcond,1:length(trialcond));
        XCsel = Y(:,:,cellinx_sel);
        rXCsel = Y(:,:,cellinx_sel) - Y0(:,:,cellinx_sel);
        XCsel = reshape(XCsel,[size(XCsel,1)*size(XCsel,2), size(XCsel,3)]);
        rXCsel = reshape(rXCsel,[size(rXCsel,1)*size(rXCsel,2), size(rXCsel,3)]);


        ORIlist = CO.get_stim_evt('ori');
        CONlist = CO.get_stim_evt('con');

        events_cont=Inf*ones(size(events_cont_inx));
        events_ORI=Inf*ones(size(events_ori_inx));
        validtrial = ~isinf(events_cont_inx(:));
        events_cont(validtrial) = CONlist(events_cont_inx(validtrial));
        events_ORI(validtrial) = ORIlist(events_ori_inx(validtrial));

        events = events_cont_inx+ length(CONlist)*(events_ori_inx-1);

        Ncell = size(Y,3);
        ORGDAT_LOC = info.session_path;
        scriptname = mfilename;
        mkdir(fileparts(fnsave));
        
        if length(set_subdata)==1
            fnsave_ses = [fnsave 'ses' num2str(inx_session) '.mat'];
            figname = [fnsave 'ses' num2str(inx_session) '.png'];            
        else
            fnsave_ses = [fnsave 'ses' num2str(inx_session) '-P' num2str(isub) '.mat'];
            figname = [fnsave 'ses' num2str(inx_session) '-P' num2str(isub) '.png'];
        end
        print(hFig,figname,'-dpng');
        if exist('eye','var')
            save(fnsave_ses,'ORGDAT_LOC','scriptname', 'datafullfn',...
            'Xsel','rXsel','XCsel','rXCsel',...
            'events','events_ORI', 'events_cont',...
            'Ncell','cellinx_sel', 'inxoutcell',...
            'inxcell_peakresptime_conts',...
            'scanlist','eye','trialinxs',...
            'param4cellselect','-v7.3');
        else
        save(fnsave_ses,'ORGDAT_LOC','scriptname', 'datafullfn',...
            'Xsel','rXsel','XCsel','rXCsel',...
            'events','events_ORI', 'events_cont',...
            'Ncell','cellinx_sel', 'inxoutcell',...
            'inxcell_peakresptime_conts',...
            'scanlist',...
            'param4cellselect','-v7.3');
        end
    end
end
end


function out = collect_data(CO, scanlist,params)

    flds = fieldnames(params);
    for i = 1 : length(flds)
        cmdstr = sprintf('%s = params.%s;',flds{i},flds{i});
        eval(cmdstr);
    end
 %---------- select cell and peakresponse_frames            
    [inxcell_peakresptime_conts,inxoutcell,...
        peakframes_spk, peakframes_dFFc, frametime, hFig] =...
        CO.select_cell_By_contrast_union(scanlist, pre_stimtime_ms, selcontrast,visual_resp_thr, motionthr, beye);

    %----------- collect data            
    ttmp1 = frametime(peakframes_dFFc);            
    twindff(1) = ttmp1(1)- diff(frametime(1:2))/2;
    twindff(2) = ttmp1(end)+ diff(frametime(1:2))/2;

    switch lower(params.animal_state)
        case {'aw','awake'}
            twins=[[-pre_stimtime_ms 0]' [0 stimtime_ms]'];
            [MXspk, ~, trialinxs] = CO.collect_trialMresp( scanlist, twins, motionthr, beye);

            twins=[[-pre_stimtime_ms 0]' twindff'*1000];
            [~, MXc] = CO.collect_trialMresp( scanlist, twins , motionthr, beye);

        case {'an','anesthesia'}
            ttmp1 = frametime(peakframes_spk);            
            twinspk(1) = ttmp1(1)- diff(frametime(1:2))/2;
            twinspk(2) = ttmp1(end)+ diff(frametime(1:2))/2;        
            twins=[[-pre_stimtime_ms 0]' twinspk'*1000];
            [MXspk, ~, trialinxs] = CO.collect_trialMresp(scanlist, twins,motionthr,beye );

            twins=[[-pre_stimtime_ms 0]' twindff'*1000];
            [~, MXc] = CO.collect_trialMresp( scanlist, twins,motionthr,beye );
        otherwise
            error('animal state correctly specified');
    end
    out.MXspk = MXspk;
    out.MXc = MXc;
    out.trialinxs = trialinxs;
    
    out.inxcell_peakresptime_conts = inxcell_peakresptime_conts;
    out.inxoutcell = inxoutcell;
    out.peakframes_spk = peakframes_spk;
    out.peakframes_dFFc = peakframes_dFFc;
    out.frametime = frametime;
    out.hFig = hFig;
    
    
end




        %---------- select cell and peakresponse_frames            
%         [inxcell_peakresptime_conts,inxoutcell,...
%             peakframes_spk, peakframes_dFFc, frametime, hFig] =...
%             CO.select_cell_By_contrast_union(scanlist, pre_stimtime_ms, selcontrast,visual_resp_thr, motionthr, beye);
%         
%         %----------- collect data            
%         ttmp1 = frametime(peakframes_dFFc);            
%         twindff(1) = ttmp1(1)- diff(frametime(1:2))/2;
%         twindff(2) = ttmp1(end)+ diff(frametime(1:2))/2;
% 
%         switch lower(info.animal_state)
%             case {'aw','awake'}
%                 twins=[[-pre_stimtime_ms 0]' [0 stimtime_ms]'];
%                 [MXspk, ~, trialinxs] = CO.collect_trialMresp( scanlist, twins, motionthr, beye);
% 
%                 twins=[[-pre_stimtime_ms 0]' twindff'*1000];
%                 [~, MXc] = CO.collect_trialMresp( scanlist, twins , motionthr, beye);
% 
%             case {'an','anesthesia'}
%                 ttmp1 = frametime(peakframes_spk);            
%                 twinspk(1) = ttmp1(1)- diff(frametime(1:2))/2;
%                 twinspk(2) = ttmp1(end)+ diff(frametime(1:2))/2;        
%                 twins=[[-pre_stimtime_ms 0]' twinspk'*1000];
%                 [MXspk, ~, trialinxs] = CO.collect_trialMresp(scanlist, twins,motionthr,beye );
% 
%                 twins=[[-pre_stimtime_ms 0]' twindff'*1000];
%                 [~, MXc] = CO.collect_trialMresp( scanlist, twins,motionthr,beye );
%             otherwise
%                 error('animal state correctly specified');
%         end


