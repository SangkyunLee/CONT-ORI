clear all
close all


if ispc
    addpath(genpath('D:\packages\codes2P'));
    addpath(genpath('Z:\codes2P'));
elseif isunix
    addpath(genpath('/home/slee/data/codes2P'));
end




%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE'};
fntmp = {'AN1-16_','AN17-22_','',''};

%----------------sh
iexp_type=2; % for iexp_type
ctm=0.6; %Ubuntu awake
npool =40;

CELLSEL_THRs=0:10:30;
SEL_METHOD='ncell';
Niter=1000;



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
        seslist =17:22;
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
end
    


currentFolder = pwd;

fnsave = sprintf('%sN-CELL_%s_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type},upper(SEL_METHOD),ctm); %L2norm -regularization
fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fnsave);
lambda1s=[1e-3 1e-2 1e-1  1 10];
lambda2s=0;








ORI_condset=combnk(ORI_list,2)';



DEC_SELCELL = cell(length(ORI_compindexset),max(seslist));



if exist('npool','var')==1
    delete(gcp('nocreate'))
    parpoolid = parpool(npool);
    pause(1);
end


%% ---------------------------------------------------
for ises = 18%seslist 
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
    
    if ises>size(DEC_SELCELL,2)
        DEC_SELCELL{length(contrasts),ises}=[];
    end
    
%     if ~isempty(DEC_SELCELL{1,ises})
%         continue;
%     end
    
    fndata = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    fullfndata = fullfile(data_path,fndata);
    fndats{ises}=fullfndata;
    disp(['fndata= ' fullfndata]);
    load(fullfndata);
    Ncell = size(Xsel,2);
    
    
    
    

    
        
    data_de0 = Xsel; % absolute spike rates



    evts=[events_cont(:) events_ORI(:)];
    if ises==18,
        listcomp=3:6;
    else
        listcomp = 1:6;
    end
    for icomp0 = listcomp % 1:length(ORI_compindexset)
    %for icomp0 = 1 : length(ORI_compindexset)
        
        icomp = ORI_compindexset(icomp0);      
        

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
                ncell = CELLSEL_THRs(ithr);
                if ncell>opts.Nch
                    DEC_SELCELL{icomp0,ises}(ithr,1:Niter,icont)=NaN;
                    
                
                elseif ncell==0,
                    opts.out_ch=[];    
                    [out ] = par_cv_classification(data,evtsel,opts);
                    DEC_SELCELL{icomp0,ises}(ithr,1:Niter,icont)= mean(out.acc_test); 
                        
                else   
                    tmpacc = zeros(Niter,1);
                    parfor iter=1:Niter                    
                        % excluding cell based on training data contrast
                        opt1 = opts;
                        incell = randperm(opts.Nch);                    
                        incell = incell(1:ncell);                   
                        out_ch = setdiff(1:opts.Nch, incell);
                        opt1.out_ch=out_ch;

                        [out ] = cv_classification(data,evtsel,opt1);
                         tmpacc(iter)= mean(out.acc_test);
                    end
                    DEC_SELCELL{icomp0,ises}(ithr,:,icont)=tmpacc;
                end
            end
            

        end % for icont  
        
        copts.sel_cl=[];
        copts.Nch=[];
        scriptname = mfilename('fullpath');

        mkdir(fileparts(fullfnsav));
        save(fullfnsav,'DEC_SELCELL',...
        'ORI_compindexset',...
        'ORI_condset','ORI_list',...
        'CELLSEL_THRs',...
        'copts',...
        'scriptname');

    end %for icomp0 
    
    

end % ises

if exist('npool','var')==1,
    delete(parpoolid)
    pause(3);
end
