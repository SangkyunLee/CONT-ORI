classdef cont_ori_session < TP.TPSession
    % contrast-orientation experiment
    

    methods
        function self = cont_ori_session(fn_data,info)
            self = self@TP.TPSession(fn_data,info);
        end
        
        
        function evtlist = get_stim_evt(self, evtstr,inxscan )            
            % evtlist = get_stim_evt (self, scanlist, evtsort )            
            % this function assume all stim conditions identical across
            % scans
            % evtstr
            if nargin<3
                inxscan=1;
            end
                
            stimparam = self.scans(inxscan).Params.stimparam;
            
                        
            switch lower(evtstr)
                case {'orientation', 'orientations','ori'}
                    evtlist = stimparam.orientation;
                case {'contrasts','contrast','con'}
                    evtlist = stimparam.contrast;
                case {'spatialfreq','sf'}
                    evtlist = stimparam.spatialFreq;
                case {'tf','tempofreq'}
                    evtlist = stimparam.tempoFreq;
                otherwise
                    error('not implemented');
            end
            
            
        end
            
        
        
        function [inxcell_peakresptime_conts,inxoutcell, peakframes_spk, peakframes_dFFc,tf,hFig] =...
                select_cell_By_contrast_union(self,scanlist, pre_stimtime_ms, selcontrast,visual_resp_thr, motionthr, beye)
            
            if nargin<6
                motionthr = struct([]);
            end
            if nargin<7,
                beye = false;
            end
                
            
            % function select_cell_By_contrast_union(self, pre_stimtime_ms, selcontrast, motionthr)
            % select cells based on mean response in any contrasts
            evtsort={'Contrast','Orientation'};
            parsort={'msperframe','samplingfreq_NI'};
            assert(check_stim_consistency(self, scanlist,evtsort),'not consistent across scans');
            assert(check_param_consistency(self, scanlist,parsort),'not consistent across scans');
            
            firstscan = scanlist(1);
            stimparam = self.scans(firstscan).Params.stimparam;
            param = self.scans(firstscan).Params;
                
            stimtime = stimparam.stim_samplesinNI/param.samplingfreq_NI;
            blanktime = stimparam.blank_samplesinNI/param.samplingfreq_NI;            
            msperframe = param.msperframe;
            
            ncell = self.scan_ncell(1);
            %learningvalid =ones(1,ncell);
            
            
            
            
            for inx_scan2 = 1: length(scanlist)
                inx_scan = scanlist(inx_scan2);
                
                stimparam = self.scans(inx_scan).Params.stimparam;
                param = self.scans(inx_scan).Params;
                nframe_prestim = round(pre_stimtime_ms/param.msperframe);
                
                
                %-------- collect trial response timecourse
                
                % apply motion threshold and eyemove to exclude trials with heavy
                % motion and eye movements
                [X, Xc, seltrial] = self.collect_trialresp_singlescan(inx_scan, pre_stimtime_ms, motionthr, beye);



                trial_cont = stimparam.cond.StimEventLUT(seltrial,1); 
                contrast_list = stimparam.contrast;

                for icontrast = 1:length(selcontrast)
                    Vcont = selcontrast(icontrast);
                    inx =find(contrast_list==Vcont);
                    inxtrial_contI = find(trial_cont==inx);



                    Xspk_m = squeeze(mean(X{1}(inxtrial_contI,:,:),1));
                    Xc_m = squeeze(mean(Xc{1}(inxtrial_contI,:,:),1));
                    if inx_scan2 == 1,
                        mrespN = zeros(size(Xc_m,1),ncell, length(selcontrast), length(scanlist));
                        mrespS = zeros(size(Xspk_m,1),ncell, length(selcontrast), length(scanlist));
                    end
                    mrespN(:,:,icontrast,inx_scan2) = Xc_m(:,1:ncell);
                    mrespS(:,:,icontrast,inx_scan2) = Xspk_m(:,1:ncell);
                end
                
                
%                 if isfield(self.scans(inx_scan),'info') && isfield(self.scans(inx_scan).info,'history')
%                     for icell = 1:ncell
%                         learningvalid(icell) = learningvalid(icell) &self.scans(inx_scan).info.history(icell);
%                     end
%                 else
%                     learningvalid1 = sign(sum(self.scans(inx_scan).Nhat,1));
%                     for icell = 1:ncell
%                         learningvalid(icell) = learningvalid(icell) & learningvalid1(icell);
%                     end
%                 end                
            end   
            
            learningvalid = self.get_cell_spkestval(scanlist);
            mrespS = mean(mrespS,4); %spike
            mrespN = mean(mrespN,4); %dF/F
            
            %------------------------
            inxincell = find(learningvalid);
            for icont = 1:length(selcontrast)
                smXc1 = std(mrespS(:,:,icont),0,1);  
                ab= gamfit(smXc1(learningvalid==1));
                % ecdf(smXc1);
                % xx = linspace(0,double(max(smXc1)));
                % line(xx,gamcdf(xx,ab(1),ab(2)),'Color','r')

                thr_std = gaminv(0.95,ab(1),ab(2));    
                % peakframes_spk = find(tf>=0 & tf<=stimtime);
                inxincell= intersect(inxincell,find(mean(mrespS(:,:,icont),1)>0 & smXc1<thr_std));
            end
            inxoutcell = setdiff(1:ncell, inxincell);
            
            
            
            
            
            %----- select frames for stimulation
            trialtime = stimtime+blanktime;
            nframe_trial = round(trialtime/msperframe*1000);
    
            tf = (-nframe_prestim:nframe_trial-1)*msperframe/1000;
            nframe_stim = ceil(stimtime/msperframe*1000);


            inxcell_peakresptime_conts = cell(1,length(selcontrast));
            hFig= figure;
            clr={'b','g','r'};
            legendstr = cell(1,length(selcontrast));
            for icont = 1:length(selcontrast)
                X0= mean(mrespS(1:nframe_prestim,:,icont),1);
                Xp = mean(mrespS(nframe_prestim+2:nframe_prestim+nframe_stim,:,icont),1);
                inxcell_peakresptime_icont =find(Xp>(visual_resp_thr*X0));
                inxcell_peakresptime_icont = intersect(inxcell_peakresptime_icont, inxincell);
                if icont==1,
                     inxcell_peakresptime = inxcell_peakresptime_icont;
                else
                    inxcell_peakresptime = union(inxcell_peakresptime, inxcell_peakresptime_icont);
                end
                inxcell_peakresptime_conts{icont} = intersect(inxcell_peakresptime_icont, inxincell);
                
                if icont==1,
                    avgresp_spk = zeros(length(tf),length(selcontrast));            
                    avgresp_dFFc = zeros(length(tf),length(selcontrast));            
                end
                avgresp_spk(:,icont) = mean(mrespS(:,inxcell_peakresptime_icont,icont),2);            
                avgresp_dFFc(:,icont) = mean(mrespN(:,inxcell_peakresptime_icont,icont),2);

                subplot(211);hold on;plot(tf,avgresp_spk(:,icont),'.-','Color',clr{icont});
                ht = title(['SPIKE, visualresp thr:' num2str(visual_resp_thr)]);
                subplot(212);hold on;plot(tf,avgresp_dFFc(:,icont),'.-','Color',clr{icont});
                title('dF/F');
                legendstr{icont}=num2str(selcontrast(icont));
            end
            legend(legendstr);
            ratio =length(inxcell_peakresptime)/length(inxincell);
            fprintf('No. selected cell=%d, ratio=%02f,',...
                length(inxcell_peakresptime),ratio);
            set(ht,'String',['SPK, Vresp thr:' num2str(visual_resp_thr) 'Cellratio: ' num2str(ratio)]);
            
            scale = zeros(length(selcontrast),1);
            for icont = 1 : length(selcontrast)
                scale(icont) = length(inxcell_peakresptime_conts{icont});
            end
            
            switch upper(self.animal_state)
                case {'AW','AWAKE'}
                    peakframes_spk = find(tf>0 & tf<=stimtime+0.05);
                otherwise                    
                    [~, mi1]=max(mean(avgresp_spk*diag(scale),2));
                    if mi1<nframe_prestim
                        peakframes_spk = find(tf>0 & tf<=stimtime+0.05);
                    else
                        ntf =tf-tf(mi1);
                        peakframes_spk= find(abs(ntf)<=stimtime/2);
                    end
            end

            [~, mi2]=max(mean(avgresp_dFFc*diag(scale),2));            
            ntf =tf-tf(mi2);
            peakframes_dFFc = find(abs(ntf)<=stimtime/2);
            
        end
        
        
        function trialcond = collect_trial_conditions(self, trialinxs, scanlist, evtsort)
            %function trialcond = collect_trial_conditions(self, trialinxs, scanlist, evtsort)
            assert(max(scanlist)<= self.no_scan, ['Total number of scan is ' num2str(self.no_scan)]);
            
            emptyscan = find(cellfun(@isempty,trialinxs));   
            emptyscan1 = intersect(scanlist, emptyscan);            
            assert(isempty(emptyscan1),sprintf('scan_notrial: %s',num2str(emptyscan1(:)')));
%             msg = sprintf('%s:\nNO SELECT TRIAL in SCAN: %s',self.session_path, num2str(emptyscan1(:)'));
%             warndlg(msg);
%             
            trialcond=struct([]);
            ns = length(scanlist);
            for i0 = 1 : ns
                i = scanlist(i0);
                param = self.scans(i).Params;
                stimparam = param.stimparam;
                cond = stimparam.cond;
                
                
                nsortevt = length(cond.event);
                for ies1 = 1 : length(evtsort)
                    for ies2 = 1 : nsortevt
                        
                        if strcmp(evtsort{ies1},cond.event(ies2).type)                                           
                            trialcond(i0).(evtsort{ies1}) =...
                                cond.StimEventLUT(trialinxs{i0},ies2);
                        end
                        
                    end
                end
            end    
            
            
        end
        
        
    end
  
    
    methods (Static)
            
        function [Y, events_CONT, events_ORI] = concatenate_trialMresp(X,trialcond,scanlist)
            %function [Y, events_cont, events_ORI] = concatenate_trialMresp(X,trialcond,scanlist)
            %   X: cell variable nx1 : n for scans
            %       within each cell TxC ; T: no. trial, C: no. cell
            %   trialcond: 1xn array of struct, n for scans
            %       fields in each struct: contrast and orientation
            %   OUTPUT:
            %       Y: T x n x C
            %       events_CONT: T x n
            %       events_ORI: T x n
            
            
            
            assert(length(X)==length(trialcond),'data size should be the same as trialcond size');

            maxtrial=0;
            for iscan=scanlist
                maxtrial = max(maxtrial, size(X{iscan},1));
            end
            minncell=Inf;
            for iscan=scanlist
                minncell = min(minncell, size(X{iscan},2));
            end
            
            Y = zeros(maxtrial,length(X),minncell);            
            events_CONT = Inf*ones(maxtrial,length(X));
            events_ORI = Inf*ones(maxtrial,length(X));
            
            
            for iscan =scanlist
                Xtmp = X{iscan};
                fnames = fieldnames(trialcond(1));
                
                for inx_field = 1 : length(fnames)
                    switch lower(fnames{inx_field})
                        case {'orientation', 'orientations','ori'}
                            ori_list = trialcond(iscan).(fnames{inx_field});
                            events_ORI(1:length(ori_list),iscan) = ori_list;
                            
                        case {'contrasts','contrast'}
                            cont_list = trialcond(iscan).(fnames{inx_field});
                            events_CONT(1:length(cont_list),iscan) = cont_list;
                        otherwise
                            error('not implemented');
                    end
                end      
                Y(1:size(Xtmp,1),iscan,:) = Xtmp(:,1:minncell);
            end            
        end
        
    end
    
end