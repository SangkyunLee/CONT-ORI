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
    fnpf='P1';
    fnpf0='P2';
    
%     fnpf='P2';
%     fnpf0='P1';
else
    fnpf=''
end

%--- ubuntu
ctm=0.6; 
iexp_type=5;
npool=40;
   

%----------------------
cell_sel_method = 'UNION_CONTRSP'; 
if iexp_type==5
DATA_thr_str = 'thr5_eyethr_xy1_p1';
else
DATA_thr_str = 'thr5'
end
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);



if iexp_type==1
    contrasts=[100 40 20];
    ORI_list=[-15 0 30 90];
    sel_compset = 1:6;
    nses=22;
    
    seslist = [(1:8) 11 12 15 16]; 
%     seslist=[12 15 16];%
elseif iexp_type==2    
    contrasts=[100 40];
    ORI_list=[-10 0 30 90];
    sel_compset = 1:6;
    nses=22;
   
    seslist = 17:22; 
elseif iexp_type==3
    
    contrasts=[100 40];
    ORI_list=[0 30 35 60 90 120 150];
    sel_compset = [1 2 3 4 7 12 16 19 21];   
    nses=7;
%     seslist =[1 7];
elseif iexp_type== 4
%     contrasts=[100 40 20];
%     ORI_list=[-10 0 30 90];
%     sel_compset = 1:6;
%     nses=17;
    
    contrasts=[100 40];
    ORI_list=[-10 0 30 90];
    sel_compset = 1:6;
    nses=40;
    % seslist =[23 (25:30) 32 33 (35:37) 40];    
%seslist =[ 25 28 29 30 32 33 (35:37) 40]; 
seslist=[11 22 23 26 27]

elseif iexp_type==5
    contrasts=[100 40];
    ORI_list=[-10 0 30 90];
    sel_compset = 1:6;
    nses=40;
    seslist = 23%[23 25 26 27 29 30 32 33 36 40];
    
end


ncellgrp = [5 4 3];
ncellgrp_current = [3 4 5];






acc_precision =0;
fndats = cell(nses,1);
%% ---------------------------------------------------
for ises = seslist
    condset=combnk(ORI_list,2)'; 
    deaccALL = Inf*ones(size(condset,2),length(contrasts));
    deaccN=cell(max(ncellgrp),size(condset,2),length(contrasts));    
    opts_grp_list = cell(size(condset,2),length(contrasts));    
    opts_ac = cell(size(condset,2),length(contrasts));  
    
    if bpart
        fnsave = sprintf('%s-DEC_CELLGRP_ctm%0.2f_ses%d.mat',fnpf,ctm,ises); %L1norm -regularization
    else
        fnsave = sprintf('DEC_CELLGRP_ctm%0.2f_ses%d.mat',ctm,ises); %L1norm -regularization
    end
    fullfnsav = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fnsave);
    %---------------------
    disp(['ises: ' num2str(ises)]);
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
    %---------------
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
    

        
    %----------------------
%     fndata = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
%     fullfndata = fullfile(data_path,fndata);
%     fndats{ises}=fullfndata;
%     disp(['fndata= ' fullfndata]);
%     load(fullfndata);    
    Ncell = length(cellinx_sel);
    
    
    
    
    
    if strfind(exp_type{iexp_type},'AWAKE')
        events_ORI(events_ORI(:)==-15)=-10;
    end
    
    
    fullfnsav1 = strrep(fullfnsav,'DEC_CELLGRP','CELLGRP');
    if exist(fullfnsav1,'file')
        load(fullfnsav1);
    else
        grplistALL = cell(1, max(ncellgrp));
        for incell = ncellgrp
            disp(incell)    
            grplistALL{incell} = VChooseK(uint8(1:Ncell), incell);
        end
        if ~exist(fileparts(fullfnsav1),'dir')
            mkdir(fileparts(fullfnsav1));
        end
        save(fullfnsav1,'grplistALL','cellinx_sel','fullfndata','-v7.3')
    end
    


%     % %------------ check whether results were saved --------------------
    check_completed = zeros(max(ncellgrp),size(condset,2),size(contrasts,2));
    if exist(fullfnsav,'file')
        load(fullfnsav)    
        for incell = ncellgrp
            for icont = 1: length(contrasts)
                for icomp = 1 : size(condset,2)
                    if icomp<= size(deaccN,2) && ~isempty(deaccN{incell,icomp,icont})
                        check_completed(incell,icomp,icont)=1;
                    end
                end
            end
        end
    end



    % %------------------------------------------------
%     if isunix  
%         delete(gcp('nocreate'))
%         parpoolid = parpool(npool);
%         pause(1);
%     elseif ispc
%         delete(gcp)
%         parpoolid = parpool(npool);
%         pause(1);
%     end
deaccALL = Inf*ones(size(condset,2),length(contrasts));
   
    evts=[events_cont(:) events_ORI(:)];
    data_de0 = Xsel;
    for incell = ncellgrp_current
        fprintf('\nincell:%d',incell)        
        grplist = grplistALL{incell};

        condlabel = cell(size(condset,2), length(contrasts));
        for icont = 1 : length(contrasts)
            fprintf('icont:%d',icont);
            test_cont = contrasts(icont);            
            for icomp = 1 : size(condset,2)    
                condlabel{icomp,icont}=sprintf('Cont:%d,dir=%dvs%d',...
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

                

                data_de = reshape(data_de,[size(data_de,1) 1 size(data_de,2)]); 


                copts.sel_cl = condset(:,icomp)';
                copts.mode =1; % mean response decoder
                copts.inxsample = 1;            
                copts.Ncv = 5;    
                copts.Nch = Ncell;
                copts.out_ch = [];                        
                copts.cvmode = 'kfold';            
                copts.pertest = 0.2;
                copts.perval = 0;

                if isinf(deaccALL(icomp,icont))

                    opts = copts;
                    opts.classifier = 'SL2LDA';
                    opts.lambda1s = 10.^(linspace(-10,1,20));
                    opts.lambda2s = 0;
                    opts.perval = 0.1;
                    opts.pertest = 0.1;
                    opts.Ncv =10;
                    
                    try
                        [out b] = cv_classification(data_de,evtsel,opts);
                    catch ME
                        out.acc_test = Inf;                        
                    end
                    deaccALL(icomp,icont) = mean(out.acc_test);        
                    if isempty(opts_ac{icomp,icont})
                        opts_ac{icomp,icont} = opts;
                    end
                end
                

                
                if check_completed(incell,icomp,icont)
                    continue;
                end
                
                opts_grp = copts;
                opts_grp.classifier = 'LDA_opt1';
                opts_grp.acc_precision = acc_precision;


                out1 = cv_grouping_classification(data_de,evtsel,opts_grp,grplist);            
                % to save space, I save uint16 with 3 digits, later on for real
                % calculation, it should be divided by 10
                if acc_precision>0
                    deaccN{incell,icomp,icont}=uint16(mean(out1.acc_test,2));
                else
                    deaccN{incell,icomp,icont}=uint8(mean(out1.acc_test,2));
                end
                
                if isempty(opts_grp_list{icomp,icont})
                    opts_grp_list{icomp,icont} = opts_grp;
                end
                


            end
        end
        
    end
    save(fullfnsav,'deaccN','deaccALL','Ncell','condset','ncellgrp','cellinx_sel','condlabel','opts_grp_list','opts_ac','-v7.3')
    %delete(parpoolid)
    pause(3);
end




