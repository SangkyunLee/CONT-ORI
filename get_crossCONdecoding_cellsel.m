clear all



if ispc
    addpath(genpath('D:\packages\codes2P'));
    addpath(genpath('Z:\codes2P'));
elseif isunix
    addpath(genpath('/home/slee/data/codes2P'));
end

exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE' };
fntmp = {'AN1-16_','AN17-22_','','','AW23-40_'};

%----------------
iexp_type=2;
ctm=0.6; %Ubuntu awake, ctm0.6
npool =6
 
%-------------------

cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5';
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
bnormaldata = true;

bpart=true;
if bpart
    fnpf='P1'
else
    fnpf=''
end


if iexp_type==1
    nses=22;
    ORI_compindexset = 1:6;
    
    contrasts=[100 40 20];
    CONT_condset={{100, 40},{100, 20},{40,100},{40, 20},{20,100},{20, 40},...
       {[100 40],100},{[100 40],40},{[100 40 20],100},{[100 40 20],40},{[100 40 20],20}};        
    ORI_list=[-15 0 30 90];
    seslist =15:16%[(1:8) 11 12 15 16];
elseif iexp_type==2
    nses=22;
    ORI_compindexset = 1:6;
    contrasts=[100 40];
    CONT_condset={{100, 40},{},{40,100},{},{},{},...
        {[100 40],100},{[100 40],40},{},{},{}};    
    ORI_list=[-10 0 30 90];
    seslist = 17 : 22;

elseif iexp_type==3
    
    contrasts=[100 40];
    CONT_condset={{100, 40},{40,100},{[100 40],100},{[100 40],40}};
    ORI_list=[0 30 35 60 90 120 150];
    ORI_compindexset = [1 2 3 4 7 12 16 19 21];   
    nses=7;
    seslist =1:7;
elseif iexp_type==4
    contrasts=[100 40 20];
    CONT_condset={{100, 40},{100, 20},{40,100},{40, 20},{20,100},{20, 40},...
        {[100 40],100},{[100 40],40},{[100 40 20],100},{[100 40 20],40},{[100 40 20],20}};
    ORI_list=[-10 0 30 90];
    ORI_compindexset = 1:6;
    nses=17;
    seslist =1:17;
elseif iexp_type==5
    contrasts=[100 40];
    CONT_condset={{100, 40},{100, 20},{40,100},{40, 20},{20,100},{20, 40},...
        {[100 40],100},{[100 40],40},{[100 40 20],100},{[100 40 20],40},{[100 40 20],20}};
    ORI_list=[-10 0 30 90];
    ORI_compindexset = 1:6;
    nses=40;
    seslist =[23 (25:30) 32 33 (35:37) 40]; 
    
end
if bpart
    fnsave = sprintf('%s-%sCRCON_SMLR_L1_ctm%0.2f.mat',fnpf,fntmp{iexp_type},ctm); 
else
    fnsave = sprintf('%sCRCON_SMLR_L1_ctm%0.2f.mat',fntmp{iexp_type},ctm); 
end
fullfnsav = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fnsave);
lambda1s=0;
lambda2s=[1e-3 1e-2 1e-1  1 10];




fndats = cell(nses,1);
CELLINX4DECODING = cell(nses,1);
CONT4DECODING = cell(nses,1);
deaccALL = cell(nses,1);
validCV = cell(nses,1);


ORI_condset=combnk(ORI_list,2)';
copts.mode =1; % mean response decoder
copts.inxsample = 1;            
copts.Ncv = 100;    
copts.out_ch = [];                        
copts.cvmode = 'random';            
copts.pertest = 0.1;
copts.classifier = 'smlr';

%% ------------- run actual calculation--------------------------------------

for ises = seslist    
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
    
    if ises>length(CELLINX4DECODING)
        fndats{ises}=[];
        CELLINX4DECODING{ises}=[];
        CONT4DECODING{ises}=[];
        deaccALL{ises}=[];
        validCV{ises}=[];
    end
    
    if ~isempty(CELLINX4DECODING{ises})
        continue;
    end
    
    if exist('npool','var')==1
        delete(gcp('nocreate'))
        parpoolid = parpool(npool);
        pause(1);
    end
        
    if bpart
        fndata = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf);
    else
        fndata = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    end
    fullfndata = fullfile(data_path,fndata);
    fndats{ises}=fullfndata;
    disp(['fndata= ' fullfndata]);
    load(fullfndata);
    Ncell = size(Xsel,2);
    
    
    
    
    CELLINX4DECODING_isub = zeros(Ncell, length(ORI_condset), length(CONT_condset));
    CONT4DECODING_isub = zeros(Ncell, length(ORI_condset), length(CONT_condset));
    deaccALL_isub = zeros(length(CONT_condset),length(ORI_condset));
    validCV_isub = zeros(length(CONT_condset),length(ORI_condset));
    

    
  
    data_de0 = Xsel; % absolute spike rates
    if bnormaldata
        scale = sqrt(sum(data_de0.^2,1));
        data_de0=bsxfun(@rdivide,data_de0,scale);
    end



    evts=[events_cont(:) events_ORI(:)];
    for icont = 1 : length(CONT_condset)
        sel_conts = CONT_condset{icont};
        if isempty(sel_conts), continue; end
        for icomp = ORI_compindexset           
            
            fprintf('icont:%d, icomp:%d',icont, icomp);
            
            

            sevts{1} = sel_conts{1};
            sevts{2} = ORI_condset(:,icomp);
            inx_train = TP.select_subdata(evts,sevts);
            % similar samples between two conditions,not exceed 10% more samples               
            inx2 = TP.select_samples_evtratio(events_ORI(inx_train),sevts{2}, 1.1);
            inx_train = inx_train(inx2);
            data_train = data_de0(inx_train,:)';
            evtsel_train = events_ORI(inx_train);
            
            
            sevts{1} = sel_conts{2};
            sevts{2} = ORI_condset(:,icomp);
            inx_test = TP.select_subdata(evts,sevts);
            data_test = data_de0(inx_test,:)';
            evtsel_test = events_ORI(inx_test);
            

            common_sampleinx = intersect(inx_train, inx_test);
            inxmap1 = zeros(max(inx_train),1);                
            inxmap1(inx_train,1)=1:length(inx_train);

            inxmap2 = zeros(max(inx_test),1);                
            inxmap2(inx_test,1)=1:length(inx_test);
            
            % prevent not to reuse same samples in training and testing
            common_sampleinx_train = inxmap1(common_sampleinx,1);
            common_sampleinx_test = inxmap2(common_sampleinx,1);





            data_train = reshape(data_train,[size(data_train,1) 1 size(data_train,2)]); 
            data_test = reshape(data_test,[size(data_test,1) 1 size(data_test,2)]); 

            opts = copts;
            opts.sel_cl = ORI_condset(:,icomp)';
            opts.perval = min(0.1*length(evtsel_test)/length(evtsel_train),0.1);
            opts.pertrain = min(0.8*length(evtsel_test)/length(evtsel_train),0.8);
            opts.common_sampleinx_train = common_sampleinx_train;
            opts.common_sampleinx_test = common_sampleinx_test;
            opts.Nch = Ncell;

            
            if exist('npool','var')==1               
                [out ] = par_cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts);   
            else
                [out ] = cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts);   
            end
            [I J inx_valid]=select_cell_byweight(out);            

            CELLINX4DECODING_isub(:,icomp,icont) = cellinx_sel(J);
            CONT4DECODING_isub(:,icomp,icont) =I; % contribution of single cells in population decoding            
            deaccALL_isub(icont, icomp) = mean(out.acc_test(inx_valid)); 
            validCV_isub(icont, icomp) = length(inx_valid);            

        end
    end
    
    CELLINX4DECODING{ises} = CELLINX4DECODING_isub;
    CONT4DECODING{ises} = CONT4DECODING_isub;
    deaccALL{ises}=deaccALL_isub;
    validCV{ises} = validCV_isub;
    opts.sel_cl=[];
    opts.Nch=[];    
    mkdir(fileparts(fullfnsav));
    save(fullfnsav,'opts','CELLINX4DECODING','CONT4DECODING','deaccALL',...
        'validCV','fndats',...
        'CONT_condset','ORI_condset','ORI_list','contrasts',...
        'ORI_compindexset','bnormaldata');
    delete(parpoolid)
pause(3);

    
end

