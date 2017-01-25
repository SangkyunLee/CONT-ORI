function get_subcell_crsdecoding(iexp_type,fnpf1,fnpf2,DATA_thr_str,bshuffle,zmean)
% function get_subcell_crsdecoding(iexp_type,fnpf1,fnpf2)
% get subcell decoding cross sessions
if nargin<5
    bshuffle = false;
end

parpoolid = set_env(true);
[contrasts, ORI_list, ORI_compindexset, nses, seslist] =get_expinfo(iexp_type);

%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
fntmp = {'AN1-16','AN17-22','','','AW23-40'};
%----------------
dtype ='Xsel';
ctm=0.6; 

CELLSEL_THRs=[0 1 3 5 10 20];
SEL_METHOD={'smlrW','ncell'};

% CELLSEL_THRs=[0 1  3 5 10 20];
% SEL_METHOD={'DECSCR','ncell'};
% igrp=3;
%-------------------

bias = false;
if ~exist('zmean','var') || isempty(zmean)
    zmean = false;
end
%--------------------



cell_sel_method = 'UNION_CONTRSP'; 
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);

currentFolder = pwd;




if bshuffle
    fnsave = sprintf('SHUFFLE_%s-%s_%sCELLSEL_%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fnpf1,fnpf2,fntmp{iexp_type},upper(SEL_METHOD{1}),ctm); %L2norm -regularization
else
    fnsave = sprintf('%s-%s_%sCELLSEL_%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fnpf1,fnpf2,fntmp{iexp_type},upper(SEL_METHOD{1}),ctm); %L2norm -regularization
end
if bias
    fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fnsave);
else
    if ~bias && zmean
        fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING_NOBIAS_ZMEAN',exp_type{iexp_type},DATA_thr_str, fnsave);
    else
        fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING_NOBIAS',exp_type{iexp_type},DATA_thr_str, fnsave);
    end
end

% fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fnsave);
lambda1s=10.^(linspace(-10,1,20));
lambda2s=0;



copts.bshuffle = bshuffle;
copts.sel_cl = [];
copts.mode =1; % mean response decoder
copts.inxsample = 1;            
copts.Ncv = 100;    
copts.Nch = [];
copts.out_ch = [];                        
copts.cvmode = 'random';            
copts.pertest = 0.1;   
copts.perval = 0.1;
copts.pertrain = 0.8;
copts.classifier = 'smlr';
copts.lambda1s = lambda1s;                        
copts.lambda2s = lambda2s;
copts.bias =bias; % false; %----------->NEW_DECODING_NOBIAS
copts.zmean = zmean; % zeromean for ex. w'*(x-x0)



ORI_condset=combnk(ORI_list,2)';

if strcmp(SEL_METHOD{1},'DECSCR') && strcmp(dtype,'Xsel')
    %SUMMARY_P1-AN17-22_DECSCR_CELLGRP_ctm0.60
    Mfn0 = sprintf('SUMMARY_%s-%s_DECSCR_CELLGRP_ctm%0.2f.mat',fnpf1,fntmp{iexp_type},ctm);
    Mfn = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, Mfn0);
    load(Mfn);
    M0=S;
else
    Mfn=[];
end

DEC_SELCELL = cell(length(ORI_compindexset),nses);
CELLSEL_NSCR = cell(length(ORI_compindexset),nses);
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
    sesinfo = get_scanspersession_eye2(exp_type{iexp_type});
    sesinfo = sesinfo(ises,:);
    if ~isempty(DEC_SELCELL{1,ises}) || all(cellfun(@isempty,sesinfo))
        continue;
    end

    fndata1 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf1);
    fndata2 = sprintf('%s_ctm%0.2fses%d-%s.mat',pprotype,ctm,ises,fnpf2);
    [D1, D2]= loadData(data_path,fndata1,fndata2);        
    [D1,D2] = subdata([],D1,D2,{dtype});
    
    if strcmp(SEL_METHOD{1},'DECSCR')
        M=M0(ises);
    end
    CELLSEL_INX{ises}=D1.cellinx_sel;       
    Ncell = length(D1.cellinx_sel);
    
    if iexp_type==5
        D1.events_ORI(D1.events_ORI(:)==-15)=-10;
        D2.events_ORI(D2.events_ORI(:)==-15)=-10;
    end
    
      
    
    
        



    evt1 = [D1.events_cont(:) D1.events_ORI(:)];
    evt2 = [D2.events_cont(:) D2.events_ORI(:)];
    for icomp0 = 1 : length(ORI_compindexset)
        
        icomp = ORI_compindexset(icomp0);      
        switch SEL_METHOD{1}
            case {'DECSCR'}
                cellinx_TR = get_clinx_decscore(M, SEL_METHOD{2},CELLSEL_THRs, igrp, icomp); 
            case {'smlrW'}
                Mopt = copts;
                Mopt.lambda1s = 0;
                Mopt.lambda2s = 10.^(linspace(-10,1,20));
                nc = length(contrasts);
                cellinx_TR = cell(length(CELLSEL_THRs),nc);
                for icont = 1: nc
                    cont = contrasts(icont);
                    oris = ORI_condset(:,icomp)';
                    out = learn_smlr(D1,dtype,Mopt,cont,oris,parpoolid,zmean);
                    
                    W = mean(out.Ws(1:Ncell,1,:),3);
                    cellinx = get_clinx_smlrW(W, CELLSEL_THRs, SEL_METHOD{2});
                   cellinx_TR(:,icont)=cellinx;
                end
            otherwise
                error('%s: not implemented yet',SEL_METHOD{1});
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
      


            data_train = D1.(dtype)(inx_train,:)';
            if zmean
                data_train = bsxfun(@minus, data_train, mean(data_train,2));                
            end
            evtsel_train = D1.events_ORI(inx_train); 


                          

            data_train = reshape(data_train,[size(data_train,1) 1 size(data_train,2)]);             
    
            for ithr = 1 : length(CELLSEL_THRs)
                for icont2 = 1 : length(contrasts)                     
                    
                    % excluding cell based on training data contrast
                    
                    
                    sevts{1} = contrasts(icont2);
                    sevts{2} = ORI_condset(:,icomp);
                    inxsample = TP.select_subdata(evt2,sevts);                        
                    inx2 = TP.select_samples_evtratio(D2.events_ORI(inxsample),sevts{2}, 1.1);
                    inx_test = inxsample(inx2);

                    data_test = D2.(dtype)(inx_test,:)';
                    if zmean                        
                        data_test = bsxfun(@minus, data_test, mean(data_test,2));                    
                    end 
                    evtsel_test = D2.events_ORI(inx_test);
                    data_test = reshape(data_test,[size(data_test,1) 1 size(data_test,2)]); 
                    
                    
                    
                    opts1 = copts;
                    opts1.common_sampleinx_train = [];
                    opts1.common_sampleinx_test = [];
                    opts1.perval = 0.1;
                    opts1.pertrain = 0.8;
                    opts1.sel_cl = ORI_condset(:,icomp)';
                    opts1.Nch = Ncell;

                    
                    
                    % excluding cell based on training data contrast
                    out_ch = setdiff(1:opts1.Nch,cellinx_TR{ithr,icont});
                    opts1.out_ch=out_ch;
                    if ~isempty(parpoolid)
                        [out1 ] = par_cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts1); 
                    else
                        %[out1 ] = cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts1); 
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
    

    scriptname = mfilename('fullpath');

    mkdir(fileparts(fullfnsav));
    save(fullfnsav,'DEC_SELCELL','CELLSEL_NSCR',...
        'ORI_compindexset',...
        'ORI_condset','ORI_list',...
        'CELLSEL_INX','CELLSEL_THRs',...
        'copts','Mfn','dtype',...
        'scriptname');

end % ises

if ~isempty(parpoolid),
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
    if isempty(subinx)
        subinx = intersect(inx1,inx2);
    end

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


function cellinx = get_clinx_decscore(M, SEL_METHOD,CELLSEL_THRs, igrp, icomp)
% function cellinx = sel_cells_decscore(M, SEL_METHOD,CELLSEL_THRs, igrp, icomp)
% M:  get_decscore    
% cellinx: cell(CELLSEL_THRs, contrasts)
        
    A = max(M.C,[],1);
    b = min(M.C,[],1);
    % C 4D (cell x group x oricomp x contrast)
    C = bsxfun(@minus,M.C,b)./bsxfun(@times, ones(size(M.C)),A-b);

    Ncell = size(M.cellinx_sel,2);
    nc = size(C,4);
    cellinx = cell(length(CELLSEL_THRs),nc);
    switch lower(SEL_METHOD)
        case {'weightl1'}

        case {'ncell'}        

            for idx1 = 1: length(CELLSEL_THRs)
                for idx2 = 1: nc
                    c = C(:,igrp,icomp,idx2);

                    %---ncell==0, all cells are used
                    if CELLSEL_THRs(idx1)==0,
                        cidx = 1:Ncell;
                    else
                        [~,cidx] =sort(c,'descend');
                        cidx = cidx(1:min(CELLSEL_THRs(idx1),Ncell));
                    end
                    cellinx{idx1,idx2} = sort(cidx);
                end
            end
        otherwise
            error('not implemented');
    end
end


function out = learn_smlr (D1,dtype,copts,cont,oris,parpoolid, zmean)

    X = D1.(dtype);
    evt1 = [D1.events_cont(:) D1.events_ORI(:)];
    
    opt = copts;    
    opt.Nch = size(X,2);
    opt.sel_cl = oris(:)';
    
    sevts{1} = cont;
    sevts{2} = oris';
    inxsample = TP.select_subdata(evt1,sevts);
    % similar number of samples between two conditions,not exceed 10% more samples
    inx2 = TP.select_samples_evtratio(D1.events_ORI(inxsample),sevts{2}, 1.1);
    inx_train = inxsample(inx2);
    
    
    data_train = X(inx_train,:)';
    if zmean
        data_train = bsxfun(@minus, data_train, mean(data_train,2));                
    end
    sc = sqrt(sum(data_train.^2,2));
    data_train = bsxfun(@rdivide, data_train,sc); % normalization
    data_train = reshape(data_train,[size(data_train,1) 1 size(data_train,2)]);
    
    evtsel_train = D1.events_ORI(inx_train);
    


    if ~isempty(parpoolid)
        out = par_cv_classification(data_train,evtsel_train,opt);
    else
        out = cv_classification(data_train,evtsel_train,opt);
    end
    
end



function cellinx = get_clinx_smlrW(W, CELLSEL_THRs, SEL_METHOD)

    [~, inxc]=sort(abs(W),'descend');

    Ncell = size(W,1);

    cellinx = cell(length(CELLSEL_THRs),1);
    switch lower(SEL_METHOD)
        case {'weightl1'}

        case {'ncell'}        

            for idx1 = 1: length(CELLSEL_THRs)                

                %---ncell==0, all cells are used
                if CELLSEL_THRs(idx1)==0,
                    cidx = 1:Ncell;
                else

                    cidx =inxc(1:min(CELLSEL_THRs(idx1),Ncell));
                end
                cellinx{idx1} = sort(cidx);

            end
        otherwise
            error('not implemented');
    end
end


