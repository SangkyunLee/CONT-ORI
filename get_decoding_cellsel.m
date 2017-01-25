clear all



if ispc
    addpath(genpath('D:\packages\codes2P'));
    addpath(genpath('Z:\codes2P'));
    machinename = getenv('COMPUTERNAME');
elseif isunix
    addpath(genpath('/home/slee/data/codes2P'));
end

exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};

bpart=false;
if bpart
%     fnpf='P1';
%     fnpf0='P2';
    
    fnpf='P2';
    fnpf0='P1';
else
    fnpf=''
end

ctm=0.6; 
iexp_type=5;
npool=7;
bshuffle=false;
%----------------------
cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5_eyethr1';
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);

switch iexp_type
    case {1}
        contrasts=[100 40 20];
        ORI_list=[-15 0 30 90];
        sel_compset = 1:6;
        nses=22;
        seslist = [(1:8) (11:12) 15 16];
        if bshuffle
            fnsave = sprintf('SUFFLE_SMLR_L1_ctm%0.2f.mat',ctm); %L1norm -regularization
        else
            fnsave = sprintf('AN1-16_SMLR_L1_ctm%0.2f.mat',ctm); %L1norm -regularization
        end
        
    case {2}
        contrasts=[100 40];
        ORI_list=[-10 0 30 90];
        sel_compset = 1:6;
        nses=22;
        seslist =17:22;%[(1:8) (11:12)];
        %
        if bpart
            fnsave = sprintf('%s-AN17-22_SMLR_L2_ctm%0.2f.mat',fnpf,ctm); %L1norm -regularization
        else
            fnsave = sprintf('AN17-22_SMLR_L1_ctm%0.2f.mat',ctm); %L1norm -regularization
        end
        %fnsave = sprintf('TEST_%s-AN17-22_SMLR_L1_ctm%0.2f.mat',fnpf,ctm); %L1norm -regularization
    case {3}
        contrasts=[100 40];
        ORI_list=[0 30 35 60 90 120 150];
        sel_compset = [1 2 3 4 7 12 16 19 21];   
        nses=7;
        seslist =1:7;
        fnsave = sprintf('SMLR_L1_ctm%0.2f.mat',ctm); 
        
    case {4}
        contrasts=[100 40 20];
        ORI_list=[-15 0 30 90];
        sel_compset = 1:6;
        nses=13;
        seslist =[10 11 22 23];
        fnsave = sprintf('AW21-22_SMLR_L1_ctm%0.2f.mat',ctm); 
    case {5}
        contrasts=[100 40];
        ORI_list=[-10 0 30 90];
        sel_compset = 1:6;
        nses=40;

        %seslist=[11 22 23 26 25 27 28 29 30 32 33 (35:37) 40]
        seslist=[23 25 26 27 29 30 32 33 36 40];
        if bpart
            fnsave = sprintf('%s-SMLR_L2_ctm%0.2f.mat',fnpf,ctm); %L1norm -regularization
        else
            fnsave = sprintf('SMLR_L2_ctm%0.2f.mat',ctm);  
        end
        
        
end



% fullfnsav = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fnsave);
% lambda1s=0;
% lambda2s=10.^(linspace(-10,1,20));

fullfnsav = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fnsave);
lambda1s = 10.^(linspace(-10,1,20));
lambda2s = 0;


%------------------------------------------------



copts.bshuffle = bshuffle;
copts.mode =1; % mean response decoder
copts.inxsample = 1;            
copts.Ncv = 100;    
copts.out_ch = [];                        
copts.cvmode = 'random';            
copts.pertest = 0.1;   
copts.perval = 0.1;
copts.classifier = 'smlr';
copts.lambda1s = lambda1s;                        
copts.lambda2s = lambda2s;
condset=combnk(ORI_list,2)';


fndats = cell(nses,1);
CELLINX4DECODING = cell(nses,1);
CONT4DECODING = cell(nses,1);
deaccALL = cell(nses,1);
validCV = cell(nses,1);

if exist('npool','var')==1
    delete(gcp('nocreate'))
    parpoolid = parpool(npool);
    pause(1);
end

%% ---------------------------------------------------
for ises = seslist
    disp(['ises: ' num2str(ises)]);
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
    if ises>length(CELLINX4DECODING)
        CELLINX4DECODING{ises}=[];
        CONT4DECODING{ises}=[]; 
        deaccALL{ises}=[];
        validCV{ises}=[];
    end
    if ~isempty(CELLINX4DECODING{ises})
        continue;
    end
        
    if bpart
        fndata = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf);
        
        fullfndata = fullfile(data_path,fndata);
        fndats{ises}=fullfndata;
        disp(['fndata= ' fullfndata]);
        load(fullfndata);
        
        % scripts below identify subcells commonly obtained from two
        % sessions
        fndata0 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf0);
        fullfndata0 = fullfile(data_path,fndata0);
        tmpdata=load(fullfndata0);
        imap = zeros(1,max(cellinx_sel));
        imap(cellinx_sel)=1:length(cellinx_sel);
        cellinx_sel1 = intersect(cellinx_sel,tmpdata.cellinx_sel);
        ininx = imap(cellinx_sel1);
        Xsel = Xsel(:,ininx);
        cellinx_sel = cellinx_sel1;
    else
        fndata = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
        fullfndata = fullfile(data_path,fndata);
        fndats{ises}=fullfndata;
        disp(['fndata= ' fullfndata]);
        load(fullfndata);
    end
    
    data_de0 = Xsel; % absolute spike rates  
    Ncell = size(Xsel,2);
    
    
    
    
    
    len_dtype=1;
    CELLINX4DECODING_isub = zeros(Ncell, length(condset), length(contrasts));
    CONT4DECODING_isub = zeros(Ncell, length(condset), length(contrasts));
    deaccALL_isub = zeros(length(contrasts),length(condset));
    validCV_isub = zeros(length(contrasts),length(condset));
    


    
    %------------------ to set compatibility between -10 deg and -15deg
    if  iexp_type==5,
        events_ORI(events_ORI(:)==-15)=-10;
    end
    evts=[events_cont(:) events_ORI(:)];
    for icont = 1:length(contrasts)
        test_cont = contrasts(icont);
        condlabel = cell(1,size(condset,2));
        for icomp = 1 : size(condset,2)    
            condlabel{icomp}=sprintf('Cont:%d,dir=%dvs%d',...
                test_cont,condset(1,icomp),condset(2,icomp));
        end

        sevts{1} = test_cont;
        for icomp = sel_compset
            fprintf('icomp:%d',icomp);

            sevts{2} = condset(:,icomp);
            inxsample = TP.select_subdata(evts,sevts);


            % similar samples between two conditions,not exceed 10% more samples               
            inx2 = TP.select_samples_evtratio(events_ORI(inxsample),sevts{2}, 1.1);
            inxsample = inxsample(inx2);
            data_de = data_de0(inxsample,:)';
            evtsel = events_ORI(inxsample);


            % data nomralization
            data_de = bsxfun(@rdivide, data_de, sqrt(sum(data_de.^2,2)));
            data_de = reshape(data_de,[size(data_de,1) 1 size(data_de,2)]); 

            opts = copts;
            opts.sel_cl = condset(:,icomp)';
            opts.Nch = Ncell;
            
            if exist('npool','var')==1                
                [out ] = par_cv_classification(data_de,evtsel,opts);
            else
                [out ] = cv_classification(data_de,evtsel,opts);
            end

            [I, J, inx_valid]=select_cell_byweight(out);
%             a = squeeze(out.Ws(1:end-1,1,:));
%             inx_valid = find(sum(abs(a),1)>0);
%             a = a(:,inx_valid);
%             a1 = bsxfun(@rdivide, abs(a), sum(abs(a),1));
%             [I J]=sort(sum(a1,2),'descend');

            

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
        'validCV','fndats','condset',...
        'contrasts','ORI_list',...
        'sel_compset');
end
if exist('npool','var')==1,
    delete(parpoolid)
    pause(3);
end


