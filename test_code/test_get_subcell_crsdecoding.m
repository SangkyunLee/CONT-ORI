% clear all
% close all

function get_subcell_crsdecoding(iexp_type,fnpf1,fnpf2)



parpoolid = set_env(true);

%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE'};
fntmp = {'AN1-16','AN17-22','','','AW23-40'};
%----------------

ctm=0.6; 
% CELLSEL_THRs=0:0.05:0.3;
% SEL_METHOD='weightl1'
CELLSEL_THRs=0:5:30;
SEL_METHOD='ncell';
%-------------------
% iexp_type=2;
% fnpf1='P1';
% fnpf2='P2';
igrp=3;


cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5';
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);

[contrasts, ORI_list, ORI_compindexset, ~, seslist] =get_expinfo(iexp_type);

%SUMMARY_P1-AN17-22_DECSCR_CELLGRP_ctm0.60
Mfn0 = sprintf('SUMMARY_%s-%s_DECSCR_CELLGRP_ctm%0.2f.mat',fnpf1,fntmp{iexp_type},ctm);
currentFolder = pwd;
Mfn = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, Mfn0);


fnsave = sprintf('%s-%sCELLSEL_%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fnpf,fntmp{iexp_type},upper(SEL_METHOD),ctm); %L2norm -regularization
fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fnsave);
lambda1s=10.^(linspace(-10,1,20));
lambda2s=0;







ORI_condset=combnk(ORI_list,2)';


load(Mfn);
M0=S;

DEC_SELCELL = cell(length(ORI_compindexset),max(seslist));
CELLSEL_NSCR = cell(length(ORI_compindexset),max(seslist));
CELLSEL_INX = cell(max(seslist),1);



%% ---------------------------------------------------
for ises = seslist 
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
    
    if ises>size(DEC_SELCELL,2)
        DEC_SELCELL{length(contrasts),ises}=[];
        CELLSEL_NSCR{length(ORI_compindexset),length(contrasts),ises}=[];
        CELLSEL_INX {ises}=[];
    end
    
    if ~isempty(DEC_SELCELL{1,ises})
        continue;
    end

    fndata1 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf1);
    fndata2 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf2);
    [D1, D2]= loadData(data_path,fndata1,fndata2);
    M=M0(ises);
    
    [D1,D2] = subdata(M.cellinx_sel,D1,D2,{'Xsel'});
    
    
    Ncell = size(M.cellinx_sel,2);
    if strcmp(fntmp{iexp_type},'AW23-40_')
        events_ORI(events_ORI(:)==-15)=-10;
    end
    
      
    
    CELLSEL_INX{ises}=M.cellinx_sel;
    A = max(M.C,[],1);
    b = min(M.C,[],1);
    C = bsxfun(@minus,M.C,b)./bsxfun(@times, ones(size(M.C)),A-b);
        



    evt1 = [D1.events_cont(:) D1.events_ORI(:)];
    evt2 = [D2.events_cont(:) D2.events_ORI(:)];
    for icomp0 = 1 : length(ORI_compindexset)
        
        icomp = ORI_compindexset(icomp0);      
        


        cellinx_TR = cell(length(CELLSEL_THRs),length(contrasts));
        switch lower(SEL_METHOD)
            case {'weightl1'}
                for idx1 = 1: length(CELLSEL_THRs)
%                     bTHR_cellcontr = M4c.rCONT4DECODING>CELLSEL_THRs(idx1);            
%                     for idx2 = 1: length(contrasts)
% 
%                         cellinx1_TR= M4c.CELLINX4DECODING(bTHR_cellcontr(:,icomp,idx2),icomp,idx2);
%                         cellinxmap = zeros(max(cellinx_sel),1);
%                         cellinxmap(cellinx_sel)=1:length(cellinx_sel);
%                         cellinx_TR{idx1,idx2} = cellinxmap(sort(cellinx1_TR));
%                     end
                end
            case {'ncell'}        
        
                for idx1 = 1: length(CELLSEL_THRs)
                    for idx2 = 1: length(contrasts)                        
                        c = C(:,igrp,icomp,idx2);
                        
                        %---ncell==0, all cells are used
                        if CELLSEL_THRs(idx1)==0,
                            cidx = 1:Ncell;
                        else
                            [~,cidx] =sort(c,'descend');
                            cidx = cidx(1:min(CELLSEL_THRs(idx1),Ncell));
                        end
                        cellinx_TR{idx1,idx2} = sort(cidx);
                    end
                end
            otherwise
                error('not implemented');
        end
        
        

        for icont = 1 : length(contrasts)
            %sel_conts =contrasts(icont);
            fprintf('icont:%d, icomp:%d\n',icont, icomp);

            %---------- select train samples ----------------
            sevts{1} = contrasts(icont);
            sevts{2} = ORI_condset(:,icomp);
            inxsample = TP.select_subdata(evt1,sevts);
            % similar number of samples between two conditions,not exceed 10% more samples               
            inx2 = TP.select_samples_evtratio(D1.events_ORI(inxsample),sevts{2}, 1.1);
            inx_train = inxsample(inx2);
            
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

            data_train = D1.Xsel(inx_train,:)';
            evtsel_train = D1.events_ORI(inx_train); 


            opts = copts;                            
            opts.perval = 0.1;
            opts.pertrain = 0.8;
            data_train = reshape(data_train,[size(data_train,1) 1 size(data_train,2)]);             
    
            for ithr = 1 : length(CELLSEL_THRs)
                for icont2 = 1 : length(contrasts)                     
                    
                    % excluding cell based on training data contrast
                    
                    
                    sevts{1} = contrasts(icont2);
                    sevts{2} = ORI_condset(:,icomp);
                    inxsample = TP.select_subdata(evt2,sevts);                        
                    inx2 = TP.select_samples_evtratio(D2.events_ORI(inxsample),sevts{2}, 1.1);
                    inx_test = inxsample(inx2);

                    data_test = D2.Xsel(inx_test,:)';
                    evtsel_test = D2.events_ORI(inx_test);
                    data_test = reshape(data_test,[size(data_test,1) 1 size(data_test,2)]); 
                    opts1 = copts;
                    opts1.common_sampleinx_train = [];
                    opts1.common_sampleinx_test = [];
                    opts1.perval = min(0.1*length(evtsel_test)/length(evtsel_train),0.1);
                    opts1.pertrain = min(0.8*length(evtsel_test)/length(evtsel_train),0.8);

                    % excluding cell based on training data contrast
                    out_ch = setdiff(1:opts1.Nch,cellinx_TR{ithr,icont});
                    opts1.out_ch=out_ch;
                    if exist('npool','var')==1
                        [out1 ] = par_cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts1); 
                    else
                        [out1 ] = cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts1); 
                    end
                    if length(out_ch) == opts1.Nch
                        DEC_SELCELL{icomp0,ises}(ithr,icont2,icont)=Inf;
                    else
                        DEC_SELCELL{icomp0,ises}(ithr,icont2,icont) = mean(out1.acc_test);
                    end
                    
                end
            end
            

        end % for icont  
        CELLSEL_NSCR{icomp0,ises}= cellinx_TR;
       

    end %for icomp0 
    
    copts.sel_cl=[];
    copts.Nch=[];
    scriptname = mfilename('fullpath');

    mkdir(fileparts(fullfnsav));
    save(fullfnsav,'DEC_SELCELL','CELLSEL_NSCR',...
        'ORI_compindexset',...
        'ORI_condset','ORI_list',...
        'CELLSEL_INX','CELLSEL_THRs',...
        'copts','Mfn',...
        'bnormaldata','scriptname');

end % ises

if exist('npool','var')==1,
    delete(parpoolid)
    pause(3);
end
end
%-------------------------------

function varargout = loadData(data_path,varargin)
    ninput = length(varargin);
    varargout = cell(1,ninput);
    for i = 1: ninput
        fndata = varargin{i};
        fullfndata = fullfile(data_path,fndata);
        disp(['fndata= ' fullfndata]);
        varargout{i} = load(fullfndata);
    end
end

%-------------------------------
function [D1,D2] = subdata(subinx,D1,D2,dtyps,dr)
% function [D1,D2] = subdata(M,D1,D2,dtyps)
% retrieve subdata with subinx in given dtyps
% NOTE: it does not retrieve all data automatically, 
% it retrive data only with dyps (cell type variable)
% dr specifies dimension inx of Dx.(dyps{ityp}) 
    
    if nargin<5
        dr = 2;
    end
      
    inx1 = D1.cellinx_sel;
    inx2 = D2.cellinx_sel;

    [~,~,ni1]=intersect(subinx,inx1);
    [~,~,ni2]=intersect(subinx,inx2);

    D1.cellinx_sel =D1.cellinx_sel(ni1);
    D2.cellinx_sel =D2.cellinx_sel(ni2);
    for ityp = 1 : length(dtyps)
        dtyp = dtyps{ityp};
        if dr==1,
            D1.(dtyp) = D1.(dtyp)(ni1,:);
            D2.(dtyp) = D2.(dtyp)(ni2,:);
        elseif dr==2,
            D1.(dtyp) = D1.(dtyp)(:,ni1);
            D2.(dtyp) = D2.(dtyp)(:,ni2);
        end
    end
end

