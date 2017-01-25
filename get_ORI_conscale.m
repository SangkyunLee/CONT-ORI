clear all
close all


if ispc
    addpath(genpath('D:\packages\codes2P'));
    addpath(genpath('Z:\codes2P'));
elseif isunix
    addpath(genpath('/home/slee/data/codes2P'));
end




%----------------
exp_type={'AN','AN_0TO150','AWAKE','AN_FULLORI'};

%----------------
iexp_type=4;
ctm=0.6; %Ubuntu awake
% npool =39;
%-------------
%----------------
% iexp_type=3;
% ctm=0.6; %KLUSTER
% npool =12;
%-------------
cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5';
calfunstr='fit_ab4'
% pprotype=['DATA_RIM_' cell_sel_method];
% data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
% fnsave = sprintf('NEW1_ORIsc_ctm%0.2f.mat',ctm); 


pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
fnsave = sprintf('AN1_DISK_ORIsc_ctm%0.2f_%s.mat',ctm,calfunstr); 

bnormaldata = true;


if iexp_type==1
    contrasts=[100 40 20];    
    nses=12;
    seslist =[(1:8) 11 12];
elseif iexp_type==2
    
    contrasts=[100 40];        
    nses=7;
    seslist =1:7;
elseif iexp_type==4
    contrasts=[100 40];        
    nses=1;
    seslist =1;
else
    contrasts=[100 40 20];        
    nses=17;
    seslist =1:17;
end

currentFolder = pwd;
% modelfname = sprintf('SMLR_L1_ctm%0.2f.mat',ctm); %L1norm -regularization
% fullmodelfname = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, modelfname);




fullfnsav = fullfile(fileparts(currentFolder),'GRP_data',exp_type{iexp_type},DATA_thr_str, fnsave);








% MODEL = load(fullmodelfname);
% 
% if exist('npool','var')==1
%     delete(gcp('nocreate'))
%     parpoolid = parpool(npool);
%     pause(1);
% end

ORIsc =cell(max(seslist),length(contrasts)-1);
M4c(max(seslist))=struct;
ORItun(max(seslist))=struct;


%% ---------------------------------------------------
for ises = seslist 
   
    fndata = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    fullfndata = fullfile(data_path,fndata);
    fndats{ises}=fullfndata;
    disp(['fndata= ' fullfndata]);
    data=load(fullfndata);
  
    
    
    
    
%     
%     M4c(ises).CELLINX4DECODING= MODEL.CELLINX4DECODING{ises};
%     M4c(ises).CONT4DECODING = MODEL.CONT4DECODING{ises};
%     M4c(ises).ORI_compindexset = MODEL.sel_compset;
%     M4c(ises).ORI_condset = MODEL.condset;    
%     M4c(ises).contrasts = MODEL.contrasts;
%     M4c(ises).ORI_list = MODEL.ORI_list;
%     
%     scale1 = max(M4c(ises).CONT4DECODING,[],1);
%     M4c(ises).rCONT4DECODING = bsxfun(@rdivide, M4c(ises).CONT4DECODING,scale1);
    
        
    X0= data.Xsel; % absolute spike rates
    if bnormaldata
        scale = sqrt(sum(X0.^2,1));
        X0=bsxfun(@rdivide,X0,scale);
    end

    
    

%--------- collecting data
    inxs_valtrial = data.events(:)>0 & ~isinf(data.events(:));
    unique_evt = unique(data.events(inxs_valtrial)');
    unique_evt = setdiff(unique_evt,0);
    
    evts = data.events(:);
    inxsample = cell(1,max(unique_evt));
    cons = zeros(1,max(unique_evt));
    oris = zeros(1,max(unique_evt));
    for ievt = unique_evt         
        inxsample{ievt} = TP.select_subdata(evts,{ievt});
        cons(ievt)=data.events_cont(inxsample{ievt}(1));
        oris(ievt)=data.events_ORI(inxsample{ievt}(1));
    end
    
    %------- est scale and bias between ori tunings of two contrasts
    evt_cond.inxsample = inxsample;
    evt_cond.cons = cons;
    evt_cond.oris = oris;
    evt_cond.conref = 100;
    evt_cond.concom = 40;
    ORIsc{ises,1} = cal_wori(X0,evt_cond,calfunstr);
    
    if find(contrasts==20)
        evt_cond.conref = 100;
        evt_cond.concom = 20;
        ORIsc{ises,2}= cal_wori(X0,evt_cond,calfunstr);
    end
    
    %---- est tuning curves
    mX=zeros(max(unique_evt),size(X0,2));
    sX1=zeros(max(unique_evt),size(X0,2));
    sX2=zeros(max(unique_evt),size(X0,2));
    for ievt = unique_evt
        mX(ievt,:) = mean(X0(inxsample{ievt},:),1);
        sX2(ievt,:) = std(X0(inxsample{ievt},:),0,1)/sqrt(length(inxsample{ievt}));
        sX1(ievt,:) = std(X0(inxsample{ievt},:),0,1);
    end
    evtlist = [cons' oris'];
    orders = {'descend','ascend'};
    [Out, J0] = TP.sort_evtorder(evtlist,orders);

    mresp = reshape(mX(J0,:),[length(unique(oris)) length(unique(cons)) size(mX,2)]);
    sresp2 = reshape(sX2(J0,:),[length(unique(oris)) length(unique(cons)) size(sX2,2)]);
    sresp1 = reshape(sX1(J0,:),[length(unique(oris)) length(unique(cons)) size(sX1,2)]);
    evtord = zeros(length(unique(oris)),length(unique(cons)), size(Out,2)); 
    for ievt = 1 : size(Out,2)
        evtord(:,:,ievt) = reshape(Out(:,ievt),[length(unique(oris)) length(unique(cons))]);
    end
    ORItun(ises).mresp = mresp;
    ORItun(ises).semresp = sresp2;
    ORItun(ises).stdresp = sresp1;
    ORItun(ises).evtord = evtord;
end % ises

save(fullfnsav,'ORIsc','ORItun','-v7.3') ;

if exist('npool','var')==1,
    delete(parpoolid)
    pause(3);
end
