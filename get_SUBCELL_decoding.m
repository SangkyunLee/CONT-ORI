clear all
close all


if ispc
    addpath(genpath('D:\packages\codes2P'));
    addpath(genpath('Z:\codes2P'));
elseif isunix
    addpath(genpath('/home/slee/data/codes2P'));
end




%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE' };
fntmp = {'AN1-16_','AN17-22_','','','AW23-40_'};


%----------------
iexp_type=2; % for iexp_type
ctm=0.6; %Ubuntu awake
npool =6;
% CELLSEL_THRs=0:0.05:0.3;
% SEL_METHOD='weightl1'
CELLSEL_THRs=0:5:30;
SEL_METHOD='ncell';

bpart=true;
if bpart
    fnpf='P1'
else
    fnpf=''
end

cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5';
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
bnormaldata = true;

switch iexp_type
    case {1}
        contrasts=[100 40 20];
        ORI_list=[-15 0 30 90];
        ORI_compindexset = 1:6;
        nses=16;
        seslist =[(1:8) 11 12 15 16];
    case {2}
        contrasts=[100 40];
        ORI_list=[-10 0 30 90];
        ORI_compindexset = 1:6;
        nses=22;
        seslist =[17:22];
    case {3}
        contrasts=[100 40];    
        ORI_list=[0 30 35 60 90 120 150];
        ORI_compindexset = [1 2 3 4 7 12 16 19 21];   
        nses=7;
        seslist =1:7;
    case {4}
        contrasts=[100 40 20];    
        ORI_list=[-10 0 30 90];
        ORI_compindexset = 1:6;
        nses=17;
        seslist =1:17;
    case {5}
        contrasts=[100 40];    
        ORI_list=[-10 0 30 90];
        ORI_compindexset = 1:6;
        nses=40;
        seslist =[23 (25:30) 32 33 (35:37) 40];       
        
end
    
% if iexp_type==1 || iexp_type==2
%     contrasts=[100 40 20];
%     ORI_list=[-15 0 30 90];
%     ORI_compindexset = 1:6;
%     nses=12;
%     seslist =[(1:8) 11 12];
% elseif iexp_type==3
%     
%     contrasts=[100 40];    
%     ORI_list=[0 30 35 60 90 120 150];
%     ORI_compindexset = [1 2 3 4 7 12 16 19 21];   
%     nses=7;
%     seslist =1:7;
% else
%     contrasts=[100 40 20];    
%     ORI_list=[-10 0 30 90];
%     ORI_compindexset = 1:6;
%     nses=17;
%     seslist =1:17;
% end

currentFolder = pwd;
modelfname = sprintf('%s-%sSMLR_L1_ctm%0.2f.mat',fnpf,fntmp{iexp_type},ctm); %L1norm -regularization
fullmodelfname = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, modelfname);


fnsave = sprintf('%s-%sCELLSEL_%s_SMLR_L2_ctm%0.2f.mat',fnpf,fntmp{iexp_type},upper(SEL_METHOD),ctm); %L2norm -regularization
fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fnsave);
lambda1s=[1e-3 1e-2 1e-1  1 10];
lambda2s=0;








ORI_condset=combnk(ORI_list,2)';
MODEL = load(fullmodelfname);



DEC_SELCELL = cell(length(ORI_compindexset),max(seslist));
CELLSEL_L1DECORDER = cell(length(ORI_compindexset),max(seslist));
CELLSEL_INX = cell(max(seslist),1);

if exist('npool','var')==1
    delete(gcp('nocreate'))
    parpoolid = parpool(npool);
    pause(1);
end


%% ---------------------------------------------------
for ises = seslist 
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
    
    if ises>size(DEC_SELCELL,2)
        DEC_SELCELL{length(contrasts),ises}=[];
        CELLSEL_L1DECORDER{length(ORI_compindexset),length(contrasts),ises}=[];
        CELLSEL_INX {ises}=[];
    end
    
    if ~isempty(DEC_SELCELL{1,ises})
        continue;
    end
    
    fndata = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf);
    fullfndata = fullfile(data_path,fndata);
    fndats{ises}=fullfndata;
    disp(['fndata= ' fullfndata]);
    load(fullfndata);
    Ncell = size(Xsel,2);
    
    if strcmp(fntmp{iexp_type},'AW23-40_')
        events_ORI(events_ORI(:)==-15)=-10;
    end
    
    
    
    
    
    M4c.CELLINX4DECODING= MODEL.CELLINX4DECODING{ises};
    M4c.CONT4DECODING = MODEL.CONT4DECODING{ises};
    M4c.ORI_compindexset = MODEL.sel_compset;
    M4c.ORI_condset = MODEL.condset;    
    M4c.contrasts = MODEL.contrasts;
    M4c.ORI_list = MODEL.ORI_list;
    
    scale1 = max(M4c.CONT4DECODING,[],1);
    M4c.rCONT4DECODING = bsxfun(@rdivide, M4c.CONT4DECODING,scale1);
    %bTHR_cellcontr = M4c.rCONT4DECODING>CELLSEL_THR;
    
   CELLSEL_INX{ises}=cellinx_sel;
    
    
        
    data_de0 = Xsel; % absolute spike rates
    if bnormaldata
        scale = sqrt(sum(data_de0.^2,1));
        data_de0=bsxfun(@rdivide,data_de0,scale);
    end


    evts=[events_cont(:) events_ORI(:)];
    for icomp0 = 1 : length(ORI_compindexset)
        
        icomp = ORI_compindexset(icomp0);      
        


        cellinx_TR = cell(length(CELLSEL_THRs),length(contrasts));
        switch lower(SEL_METHOD)
            case {'weightl1'}
                for idx1 = 1: length(CELLSEL_THRs)
                    bTHR_cellcontr = M4c.rCONT4DECODING>CELLSEL_THRs(idx1);            
                    for idx2 = 1: length(contrasts)

                        cellinx1_TR= M4c.CELLINX4DECODING(bTHR_cellcontr(:,icomp,idx2),icomp,idx2);
                        cellinxmap = zeros(max(cellinx_sel),1);
                        cellinxmap(cellinx_sel)=1:length(cellinx_sel);
                        cellinx_TR{idx1,idx2} = cellinxmap(sort(cellinx1_TR));
                    end
                end
            case {'ncell'}        
        
                for idx1 = 1: length(CELLSEL_THRs)

                    for idx2 = 1: length(contrasts)
                        CONT = M4c.rCONT4DECODING(:,icomp,idx2);
                        bTHR_cellcontr = false(size(CONT));
                        %---ncell==0, all cells are used
                        if CELLSEL_THRs(idx1)==0,
                            bTHR_cellcontr(:)=true;
                        else
                            bTHR_cellcontr(1:min(CELLSEL_THRs(idx1),length(CONT)))=true; 
                            
                        end

                        cellinx1_TR= M4c.CELLINX4DECODING(bTHR_cellcontr,icomp,idx2);
                        cellinxmap = zeros(max(cellinx_sel),1);
                        cellinxmap(cellinx_sel)=1:length(cellinx_sel);
                        cellinx_TR{idx1,idx2} = cellinxmap(sort(cellinx1_TR));
                    end
                end
            otherwise
                error('not implemented');
        end
        
        

        for icont = 1 : length(contrasts)
            sel_conts =contrasts(icont);
            fprintf('icont:%d, icomp:%d\n',icont, icomp);

            %---------- select samples ----------------
            sevts{1} = sel_conts;
            sevts{2} = ORI_condset(:,icomp);
            inxsample = TP.select_subdata(evts,sevts);
            % similar samples between two conditions,not exceed 10% more samples               
            inx2 = TP.select_samples_evtratio(events_ORI(inxsample),sevts{2}, 1.1);
            inxsample = inxsample(inx2);
            
            %----------------------------------
      
            
            copts.sel_cl = ORI_condset(:,icomp)';
            copts.mode =1; % mean response decoder
            copts.inxsample = 1;            
            copts.Ncv = 100;    
            copts.Nch = Ncell;
            copts.out_ch = [];                        
            copts.cvmode = 'random';            
            copts.pertest = 0.1;   
            copts.classifier = 'smlr';
            copts.lambda1s = lambda1s;                        
            copts.lambda2s = lambda2s;

            data = data_de0(inxsample,:)';
            evtsel = events_ORI(inxsample); 


            opts = copts;                            
            opts.perval = 0.1;
            opts.pertrain = 0.8;
            data = reshape(data,[size(data,1) 1 size(data,2)]);             
    
            for ithr = 1 : length(CELLSEL_THRs)
                for icont2 = 1 : length(contrasts)                     
                    % excluding cell based on training data contrast
                    
                   
                    out_ch = setdiff(1:opts.Nch,cellinx_TR{ithr,icont2});
                    opts.out_ch=out_ch;
                    if length(out_ch) == opts.Nch
                        DEC_SELCELL{icomp0,ises}(ithr,icont2,icont)=Inf;
                    else
                        [out ] = par_cv_classification(data,evtsel,opts);
                        DEC_SELCELL{icomp0,ises}(ithr,icont2,icont) = mean(out.acc_test);
                    end
                end
            end
            

        end % for icont  
        CELLSEL_L1DECORDER{icomp0,ises}= cellinx_TR;
       

    end %for icomp0 
    
    copts.sel_cl=[];
    copts.Nch=[];
    scriptname = mfilename('fullpath');

    mkdir(fileparts(fullfnsav));
    save(fullfnsav,'DEC_SELCELL','CELLSEL_L1DECORDER',...
        'ORI_compindexset',...
        'ORI_condset','ORI_list',...
        'CELLSEL_INX','CELLSEL_THRs',...
        'copts','fullmodelfname',...
        'bnormaldata','scriptname');

end % ises

if exist('npool','var')==1,
    delete(parpoolid)
    pause(3);
end
