function newinfer_spk(dff_fac, bonlyrim,bspatial, ndct_time, exptypestr,navg)

if ispc
    addpath(genpath('Z:\codes2P'));
elseif isunix
    addpath(genpath('/home/slee/data/codes2P'));
end
if nargin<6
    navg=1;
end



logfn='../doc/*.xls';
logfn = dir(logfn);
logfn=sprintf('../doc/%s',logfn.name);
sheet=1;
finfo = load_logfile(logfn, sheet);


% exptypestr = 'Cont-ORI' %'ORI100'%'ORI-cont0.58'%'Nat1'%''multibar'%'ORI100'% 'ori'%'ORI' % 
list_selT=[];
for ifinfo=1:length(finfo)
    strinx = strfind(lower(finfo(ifinfo).Experiment),lower(exptypestr));
    if length(finfo(ifinfo).Experiment)>=length(exptypestr)...
            && ~isempty(strinx) ...
            && strinx==1 ...
            && (strcmp(finfo(ifinfo).use,'1') || finfo(ifinfo).use(1)==1 )
       list_selT = [list_selT ifinfo];
    end
end


% dff_fac=3;
% ndct_time=100;
% bonlyrim=1;

tau = 0.85;
rel_tol = 0.05;
default_lambda1 = 10;
dffn_fac = 0;  
F0cutoffthr=50;

switch exptypestr
    case { 'ORI100'}
        if bonlyrim 
            fnprefix=['.' filesep 'data'  filesep 'RIM_ORI'];
        else
            fnprefix=['.' filesep 'data'  filesep 'DISK_ORI'];  
        end
    case {'ORI-cont0.58'}
        if bonlyrim 
            fnprefix=['.' filesep 'data'  filesep 'RIM_ORI058'];
        else
            fnprefix=['.' filesep 'data'  filesep 'DISK_ORI058'];    
        end    
    case {'Cont-ORI','Cont','Cont-ori'}
        if bonlyrim 

            fnprefix=['.' filesep 'data'  filesep 'RIM_Cont-ORI'];

        else
            fnprefix=['.' filesep 'data'  filesep 'DISK_Cont-ORI'];    
        end
    case {'Nat'}
        if bonlyrim 
            fnprefix=['.' filesep 'data'  filesep 'RIM_Nat'];
        else
            fnprefix=['.' filesep 'data'  filesep 'DISK_Nat'];    
        end
    case {'multibar'}
        if bonlyrim 
            fnprefix=['.' filesep 'data'  filesep 'RIM_Multibar'];
        else
            fnprefix=['.' filesep 'data'  filesep 'DISK_Multibar'];    
        end
        
    
    otherwise
        error('not implemented yet');
end

if navg>1
    [a, b] = fileparts(fnprefix);
    fnprefix=fullfile(a,sprintf('navg%d%s',navg,b));
end
    
%%
for ctm_scale = 0.7%[0.6 0 0.4]%[0.5 0 0.7]%[0.6 0.7 0.5 0]

    clear ndata;
    
    
    for idata = 1 : length(list_selT)    
        file_loc = finfo(list_selT(idata));
        data = load(['.' filesep 'data' filesep file_loc.Savefn]);
        
        
        
%         stimtime = data.Params.stimparam.stim_samplesinNI/data.Params.samplingfreq_NI;
%         blanktime = data.Params.stimparam.blank_samplesinNI/data.Params.samplingfreq_NI;
%         trialtime = stimtime+blanktime;
%         nframe_trial = round(trialtime/data.Params.msperframe*1000);
        
        Nseg =length(data.F);
        nCell = data.ROI{1}(1).NinitROI;
        
        % as{1,:}= pixel lists selected for dF/F calculation
        % as{2,:}= spatial weights for the selected pixels
        as = cell(2,nCell);        
        history = zeros(nCell,1);        
        if exist('navg','var')        
            T = floor(data.Params.Nframes/navg);
            data.Params.msperframe =  data.Params.msperframe*navg;
        else
            T = data.Params.Nframes;
        end
        Nhat = zeros(T, nCell);
        dFF = zeros(T, nCell);
        dFF2 = zeros(T, nCell);
        dFFn = zeros(T, nCell);
        
        if ndct_time>0,
            ndct = round(data.Params.msperframe/1000*T/ndct_time);            
        else
            ndct=0;
        end
        
        


        opts.dt = data.Params.msperframe/1000;
        opts.rel_tol = rel_tol;      
        opts.lambda1 = default_lambda1; 
        opts.tau = tau;
        opts.chunk = 1000;    %data from all chunks together  
            
        for icell = 1: nCell
          
            
            if length(data.ROI{1})<icell
                validcell=false;
            else
                validcell = true;
                for isubdat = 1 : Nseg
                    validcell = validcell & ~isempty(data.ROI{isubdat}(icell).pixel_list);
                end
            end
            
            if validcell
                fill_list_sub = cell(1,Nseg);
                pixel_list_sub = cell(1,Nseg);                
                for isubdat = 1 : Nseg
                    fill_list_sub{isubdat} = data.ROI{isubdat}(icell).fill_list;
                    pixel_list_sub{isubdat} = data.ROI{isubdat}(icell).pixel_list;
                    if isubdat==1,
                        compixs1 = fill_list_sub{isubdat};
                        compixs2 = pixel_list_sub{isubdat};
                    else                        
                        compixs1 = intersect(compixs1, fill_list_sub{isubdat});
                        compixs2 = intersect(compixs2, pixel_list_sub{isubdat});
                    end                    
                end
                if bonlyrim
                    compixs = compixs2;
                else
                    compixs = compixs1;
                end

                F=zeros(T,length(compixs));
                Fn=zeros(T,1);
                accumTseg=0;
                newpixinx_sub = zeros(length(compixs),Nseg);
                for isubdat = 1 : Nseg
                    mask=zeros(size(data.cci));
                    mask(compixs1)=1;
                    mask(compixs2)=2;
                    if bonlyrim
                        newpixinx_sub(:,isubdat) = find(mask(fill_list_sub{isubdat})==2);
                    else
                        newpixinx_sub(:,isubdat) = find(mask(fill_list_sub{isubdat})>0);
                    end
                    if navg>1 %exist('navg','var') 
                        for iavg = 1 :navg
                            if iavg
                                Ftmp = data.F{isubdat}{icell}(iavg:navg:end,newpixinx_sub(:,isubdat)); 
                                Fntmp = data.Fn{isubdat}(iavg:navg:end,icell);
                            else
                                Ftmp = Ftmp +data.F{isubdat}{icell}(iavg:navg:end,newpixinx_sub(:,isubdat));    
                                Fntmp = Fntmp + data.Fn{isubdat}(iavg:navg:end,icell);
                            end
                        end
                        Ftmp = Ftmp/navg;
                        Fntmp = Fntmp/navg;

                    else
                        Ftmp = data.F{isubdat}{icell}(:,newpixinx_sub(:,isubdat));
                        Fntmp = data.Fn{isubdat}(:,icell);
                    end


                    F(accumTseg+1:accumTseg+size(Ftmp,1),:) =Ftmp;
                    Fn(accumTseg+1:accumTseg+size(Ftmp,1),:) = Fntmp;

                    accumTseg = accumTseg + size(Ftmp,1);
                end

                [~, Beta, DCT]=applyhighDCTfilter(F,ndct);
                F = F- DCT(:,2:end)*Beta(2:end,:);

                [~, Beta, DCT]=applyhighDCTfilter(Fn,ndct);
                Fn = Fn - DCT(:,2:end)*Beta(2:end,:);

                F= bsxfun(@minus,F, ctm_scale*Fn);

                % calculation F0
                if bspatial
                    F01 = zeros(1,size(F,2));
                    for inxp=1:size(F,2)
                       [mu,sigma] = normfit(F(:,inxp));
                       inxf = (F(:,inxp)<(mu + dff_fac*sigma));
                       F01(inxp) = mean(F(inxf,inxp));
                    end

                    inxpix_posF0 = find(F01>F0cutoffthr);
                    dFF1 = bsxfun(@rdivide,bsxfun(@minus,F(:,inxpix_posF0),F01(:,inxpix_posF0)),F01(:,inxpix_posF0));
                    [Nhat1, dFF2, a, ~, history1]=getspikes4(double(dFF1),opts);
                    %figure; plot(mean(dFFi,2))
                else
                    mF1 =mean(F,1);
                    inxpix_posF0 = find(mF1>F0cutoffthr);
                    mF = mean(F(:,inxpix_posF0),2);
                    [mu,sigma] = normfit(mF);
                    inxf = (mF<(mu + dff_fac*sigma));
                    F01 = mean(mF(inxf));
                    dFF1 = (mF-F01)/F01;
                    opts.MAX_ITER=10;
                    [Nhat1, dFF2, a, ~, history1]=getspikes4(double(dFF1),opts);                 
                end
                atmp = zeros(size(F01));
                atmp(inxpix_posF0)=a;
                a1 = atmp;


                if ~isnan(history1.posts(end))
                    history(icell) = 1;
                end



                [mu,sigma] = normfit(Fn);
                inxf = (Fn<(mu+dffn_fac*sigma));
                Fn0 = mean(Fn(inxf));

                dFFni = bsxfun(@rdivide,bsxfun(@minus,Fn,Fn0),Fn0);
                %figure; plot(dFFni);

                dFFn(:,icell) = dFFni(1:T);
                as{1,icell} = compixs;
                as{2,icell} = a1;            
                Nhat(:,icell) = Nhat1(1:T);
                dFF(:,icell) = mean(dFF1(1:T,:),2);
                dFF2(:,icell) =mean(dFF2(1:T,:),2);
            end

        end

        
        data.Nhat = Nhat;
        data.dFF = dFF;
        data.dFF2 = dFF2;
        data.dFFn =dFFn;
        data.as = as;
        data.info.history = history;
        data.info.dff_fac=dff_fac;
        data.info.dffn_fac= dffn_fac;
        data.info.F0cutoffthr =F0cutoffthr;

        Params =data.Params;
        if exist('navg','var') 
            mirror_start_time = Params.mirror_start_time;
            newL = floor(Params.Nframes/navg);
            mirror_start_time = mirror_start_time(1:navg:newL*navg);
            Params.Nframes = newL;
            Params.mirror_start_time = mirror_start_time;
            Params.mirror_frame_start = zeros(size(Params.mirror_frame_start));
            Params.mirror_frame_start(mirror_start_time)=10;
            timeinfo = gen_stimtime_GratingExperiment2PhotonbySang(Params);
            data.timeinfo =timeinfo;
        end
            
        newParams.samplingfreq_NI = Params.samplingfreq_NI;
        newParams.msperframe = Params.msperframe;
        newParams.stimparam =Params.stimparam; 
        newParams.Nframes = Params.Nframes;
        newParams.files = Params.files;
        newParams.stimparam = Params.stimparam;
        newParams = orderfields(newParams);


        data.info.ndct_time = ndct_time;
        data.info.nhat_opts =opts;
        data= rmfield(data,'F');
        data = rmfield(data,'Fn');
        data = rmfield(data,'Params');
        data.Params = newParams;
        data.exptype = finfo(list_selT(idata)).Experiment;
    

        newdata(idata)=orderfields(data); 

    end
    data = newdata;  
    if bspatial
    fn_save = sprintf('%s_ctm%.2f_HPF%d_tau%.2f_F0-%.1fsigma.mat',fnprefix,ctm_scale,ndct_time,tau,dff_fac)
    else
        fn_save = sprintf('nospatial_%s_ctm%.2f_HPF%d_tau%.2f_F0-%.1fsigma.mat',fnprefix,ctm_scale,ndct_time,tau,dff_fac)
    end
    
    save(fn_save,'data','-v7.3')
    pause(5);
    clear data newdata;
    
end