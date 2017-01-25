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
iexp_type=1;
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

% pprotype=['DATA_RIM_' cell_sel_method];
% data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
% fnsave = sprintf('NEW1_ORIsc_ctm%0.2f.mat',ctm); 


pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
fnsave = sprintf('AN17-22_DISK_ORIsc_ctm%0.2f.mat',ctm); 

bnormaldata = true;


if iexp_type==1
% for session no. 17 throgh 22, 
    contrasts=[100 40];    
    nses=22;
    seslist =17:22;
    subsets = cell(nses,2);
    subsets(17,:) = {1:3,4:6};
    subsets(18,:) = {1:3,4:6};
    subsets(19,:) = {1:3,4:6};
    subsets(20,:) = {1:3,4:6};
    subsets(21,:) = {1:3,4:6};
    subsets(22,:) = {1:3,4:6};
    
elseif iexp_type==2    
    
elseif iexp_type==4
    contrasts=[100 40];        
    nses=2;
    seslist =2;
    %subsets = cell(nses,3);
    %subsets(1,:) = {1:3 ,4:6,7:10};
    subsets = cell(nses,2);
    subsets(2,1) = {1:3};
    
else
    contrasts=[100 40 20];        
    nses=17;
    seslist =1:17;
end

currentFolder = pwd;
% modelfname = sprintf('SMLR_L1_ctm%0.2f.mat',ctm); %L1norm -regularization
% fullmodelfname = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, modelfname);




fullfnsav = fullfile(fileparts(currentFolder),'GRP_data',exp_type{iexp_type},DATA_thr_str, fnsave);







if exist('npool','var')==1
    delete(gcp('nocreate'))
    parpoolid = parpool(npool);
    pause(1);
end

ORIsc =cell(max(seslist),3);
M4c(max(seslist),3)=struct;
ORItun(max(seslist),3)=struct;


%% ---------------------------------------------------
for ises = seslist 
   
    fndata = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    fullfndata = fullfile(data_path,fndata);
    fndats{ises}=fullfndata;
    disp(['fndata= ' fullfndata]);
    data=load(fullfndata);
  
    
    ss = subsets(ises,:);
    
    for isub = 1 : length(ss)
        subdata = data;
        inx_ns = setdiff(cell2mat(ss),ss{isub});
        subdata.events(:,inx_ns)=NaN;
        subdata.events_ORI(:,inx_ns)=NaN;
        subdata.events_cont(:,inx_ns)=NaN;
    
    
        
        X0= subdata.Xsel; % absolute spike rates
        if bnormaldata
            inxnnan = ~isnan(subdata.events(:));
            scale = sqrt(sum(X0(inxnnan,:).^2,1));
            X0=bsxfun(@rdivide,X0,scale);
        end




    %--------- collecting data
        inxs_valtrial = subdata.events(:)>0 & ~isinf(subdata.events(:)) & ~isnan(subdata.events(:));
        unique_evt = unique(subdata.events(inxs_valtrial)');
        unique_evt = setdiff(unique_evt,0);

        evts = subdata.events(:);
        inxsample = cell(1,max(unique_evt));
        cons = zeros(1,max(unique_evt));
        oris = zeros(1,max(unique_evt));
        for ievt = unique_evt         
            inxsample{ievt} = TP.select_subdata(evts,{ievt});
            cons(ievt)=subdata.events_cont(inxsample{ievt}(1));
            oris(ievt)=subdata.events_ORI(inxsample{ievt}(1));
        end

        %------- est scale and bias between ori tunings of two contrasts
        evt_cond.inxsample = inxsample;
        evt_cond.cons = cons;
        evt_cond.oris = oris;
        evt_cond.conref = 100;
        evt_cond.concom = 40;
        ORIsc{ises,isub} = cal_wori(X0,evt_cond);

%         if find(contrasts==20)
%             evt_cond.conref = 100;
%             evt_cond.concom = 20;
%             ORIsc{ises,2,isub}= cal_wori(X0,evt_cond);
%         end

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
        ORItun(ises,isub).mresp = mresp;
        ORItun(ises,isub).semresp = sresp2;
        ORItun(ises,isub).stdresp = sresp1;
        ORItun(ises,isub).evtord = evtord;
    end
end % ises

save(fullfnsav,'ORIsc','ORItun','-v7.3') ;

if exist('npool','var')==1,
    delete(parpoolid)
    pause(3);
end
