
function get_subcell_gamnoise_crcon_decoding1(iexp_type, DATA_thr_str, comcont, cellselopt, bshuffle,ctm,seslist,zmean,getDfun, parpoolid)
% function get_subcell_gamnoise_crcon_decoding1(iexp_type, DATA_thr_str, comcont, cellselopt, bshuffle,ctm,seslist,zmean,getDfun, parpoolid)
% % This function is identical to get_subcell_gamnoise_crcon_decoding,
% but to faciliate larger parallel pool, code was modified.





if nargin<5
    bshuffle = false;
end

if ~ exist('parpoolid','var')
    parpoolid = set_env(true);
    interparpool = true;
else
    interparpool = false;
end
[~, ORI_list, ORI_compindexset, nses, ses0] =get_expinfo(iexp_type);
if ~exist('seslist','var') || isempty(seslist)
    seslist = ses0;
end



%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
switch getDfun{1}
    case 'get_gamdat1'
         fntmp = {'GAM1_AN1-16','GAM1_AN17-22','','','GAM1_'};
    case {'get_gamdat2'}
        fntmp = {'GAM2_AN1-16','GAM2_AN17-22','','','GAM2_'};
    case {'get_gamdat3'}
        fntmp = {'GAM3_AN1-16','GAM3_AN17-22','','','GAM3_'};
    otherwise
end


CELLSEL_THRs = cellselopt.CELLSEL_THRs;
SEL_METHOD = cellselopt.SEL_METHOD;
dtype='Xsel';
bias = false;
if ~exist('zmean','var') || isempty(zmean)
    zmean = false;
end
%----------------
if ~exist('ctm','var')
ctm=0.6; 
end

cell_sel_method = 'UNION_CONTRSP'; 
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);


currentFolder = pwd;
if bshuffle
    fnsave = sprintf('SHUFFLE_%s_SUBCELL-%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type}, dtype,ctm); %L2norm -regularization
else
    fnsave = sprintf('%s_SUBCELL-%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type},dtype,ctm); %L2norm -regularization
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

%--------- L2-norm regularization
lambda1s=10.^(linspace(-10,1,20));
lambda2s=0;


%-------L1-norm regularization
M2_lambda1s=0;
M2_lambda2s=10.^(linspace(-10,1,20));


%------- common option
opts.bshuffle = bshuffle;
opts.mode =1; % mean response decoder
opts.inxsample = 1;            
opts.Ncv = 10;    
opts.out_ch = [];                        
opts.cvmode = 'random';            
opts.classifier = 'smlr';
opts.bias =bias; % false; %----------->NEW_DECODING_NOBIAS
opts.zmean = zmean; % zeromean for ex. w'*(x-x0)
 

% specific opts for cell selection
mopts = opts;
mopts.lambda1s = M2_lambda1s;                        
mopts.lambda2s = M2_lambda2s;
mopts.pertest = 0.1;   
mopts.perval = 0.1;
mopts.pertrain = 0.8;

% common options to be used for testing decorders
topt0 = opts;
topt0.lambda1s = lambda1s;                        
topt0.lambda2s = lambda2s;
topt0.pertest = [];   
topt0.perval = [];
topt0.pertrain = [];


ORI_condset=combnk(ORI_list,2)';


DEC_SELCELL = zeros(length(CELLSEL_THRs),length(comcont),length(ORI_compindexset),nses);
DEC_M = zeros(length(comcont),length(ORI_compindexset),nses);
CELLSEL_INX = cell(length(comcont),length(ORI_compindexset),nses);
SAMPLE_INX = cell(length(comcont),length(ORI_compindexset),nses);
WEIGHT = cell(length(comcont),length(ORI_compindexset),nses);
dim_order={'CELLSEL_THRs','comcont','ORI_compindexset'};

scriptname = mfilename('fullpath');
%% ---------------------------------------------------
for ises = seslist 
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
    
    dtmp = DEC_SELCELL(1,:,:,ises);
    if all(dtmp(:))
        fprintf('ises: %d done.\n',ises);
        continue;
    end

    fndata1 = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    [D1]= loadData(data_path,fndata1);   
    
    Ncell = size(D1.cellinx_sel,2);
    if iexp_type==5
        D1.events_ORI(D1.events_ORI(:)==-15)=-10;
    end
    
      



    evt1 = [D1.events_cont(:) D1.events_ORI(:)];
    
    ncont = length(comcont);
    noricom = length(ORI_compindexset);
    np = ncont* noricom;
    [lori, lcomcont]= meshgrid(1:noricom,1:ncont);
    
    DEC_SELCELL1 = zeros(length(CELLSEL_THRs),np);
    WEIGHT1  = cell(np,1);
    SAMPLE_INX1  = cell(np,1);
    CELLSEL_INX1  = cell(np,1);
    DEC_M1 = zeros(np,1);
    parfor j = 1: np
        icomp0 = lori(j);
        icont = lcomcont(j);
        icomp = ORI_compindexset(icomp0);     

         
            
        fprintf('ises: %d, icont:%d, icomp:%d\n',ises,icont, icomp);


        %------- prepare data
        contcond = comcont{icont};
        ORIsel = ORI_condset(:,icomp);            
        [dtr_gam,dte_gam,evttr,evtte] = get_gamdata(D1,evt1, dtype,contcond,ORIsel, getDfun); 
   


        if zmean
            dtr_gam = bsxfun(@minus, dtr_gam, mean(dtr_gam,2));                
        end
        sc = sqrt(sum(dtr_gam.^2,2));
        dtr_gam0 = bsxfun(@rdivide, dtr_gam,sc); % normalization
        dtr_gam = reshape(dtr_gam,[size(dtr_gam,1) 1 size(dtr_gam,2)]);  
        dtr_gam0 = reshape(dtr_gam0,[size(dtr_gam0,1) 1 size(dtr_gam0,2)]);


        if ~isempty(dte_gam)
            if zmean
                dte_gam = bsxfun(@minus, dte_gam, mean(dte_gam,2));                    
            end 
            
            dte_gam0 = bsxfun(@rdivide, dte_gam,sc); % normalization
            dte_gam = reshape(dte_gam,[size(dte_gam,1) 1 size(dte_gam,2)]);                  
            dte_gam0 = reshape(dte_gam0,[size(dte_gam0,1) 1 size(dte_gam0,2)]);
        end


       %-------------


        if any(comcont{icont}{1}==comcont{icont}{2}) && length(comcont{icont}{1})>1
            error('not implemented for contrast-independent decoder');

        else                   
            opts1 = mopts;
            opts1.Nch =Ncell;
            opts1.sel_cl = ORIsel';
            opts1.common_sampleinx_train = [];
            opts1.common_sampleinx_test = [];     

        end

        if ~all(contcond{1}==contcond{2})
            opts1.perval = min(opts1.perval*length(evtte)/length(evttr),opts1.perval);
            opts1.pertrain = min(opts1.pertrain*length(evtte)/length(evttr),opts1.pertrain);            
        end
        opts1.out_ch=[];

        % selecting cell with L1 - smlr
        try

%             if all(contcond{1}==contcond{2})
%                 [out1 ] = NOpar_cv_classification(dtr_gam0,evttr,opts1);
%             else
%                 [out1 ] = NOpar_cv_classification2(dtr_gam0,evttr,dte_gam0,evtte,opts1); 
%             end 
%             DEC_M1(j) = mean(out1.acc_test);
% 
%             Ws = squeeze(out1.Ws(1:Ncell,1,:));
%             cellinx_TR_ses = select_cell(Ws, CELLSEL_THRs, SEL_METHOD);

            %-------------
            opts2 = topt0;
            opts2.Nch =Ncell;
            opts2.sel_cl = ORIsel';
            
            opts2.cvmode = 'random';
            opts2.pertest = opts1.pertest;
            opts2.pertrain = opts1.pertrain;
            opts2.perval = opts1.perval;
            
            %opts2.cvmode = 'precal';
            %opts2.pertest = [];
            %opts2.pertrain = [];
            %opts2.perval = [];
            %opts2.inxs_cv = out1.inxs_cv;
            opts2.common_sampleinx_train = [];
            opts2.common_sampleinx_test = [];
            %----------------
            DECSEL = zeros(length(CELLSEL_THRs),1);
            for ithr = 1 %: length(CELLSEL_THRs)   
                %gen_outch = @(y)setdiff(1:opts2.Nch,y);
                %out_ch = cellfun(gen_outch,cellinx_TR_ses(ithr,:),'UniformOutput',false);
                out_ch =[];
                
                opts2.out_ch=out_ch;

                if all(comcont{icont}{1}==comcont{icont}{2})
                    [out2 ] = NOpar_cv_classification(dtr_gam,evttr,opts2);
                else
                    [out2 ] = NOpar_cv_classification2(dtr_gam,evttr,dte_gam,evtte,opts2);
                end

                DECSEL(ithr) = mean(out2.acc_test);
                
            end
            DEC_SELCELL1(:,j) = DECSEL;

        catch
            continue;
        end            
 
        SAMPLE_INX1{j} = out2.inxs_cv;
        %SAMPLE_INX1{j} = opts2.inxs_cv;
        
        %WEIGHT1{j} = Ws;
        %CELLSEL_INX1{j} = cellinx_TR_ses;    

    end 
     
    %DEC_M(:,:,ises) = reshape(DEC_M1, [ncont noricom]);
    DEC_SELCELL(:,:,:,ises) = reshape(DEC_SELCELL1, [length(CELLSEL_THRs) ncont noricom]);
    %WEIGHT(:,:,ises) = reshape(WEIGHT1, [ncont noricom]);
    SAMPLE_INX(:,:,ises) = reshape(SAMPLE_INX1, [ncont noricom]);
    %CELLSEL_INX(:,:,ises) = reshape(CELLSEL_INX1, [ncont noricom]);
    
    mkdir(fileparts(fullfnsav));
    
    if exist(fullfnsav,'file')
       M0 = load(fullfnsav);
       pause(0.1);
        M0.DEC_SELCELL(:,:,:,ises)= DEC_SELCELL(:,:,:,ises);
        M0.CELLSEL_INX(:,:,ises) = CELLSEL_INX(:,:,ises);
        M0.DEC_M(:,:,ises) = DEC_M(:,:,ises);
        M0.WEIGHT = WEIGHT;

        DEC_SELCELL = M0.DEC_SELCELL;
        CELLSEL_INX = M0.CELLSEL_INX;
        DEC_M = M0.DEC_M;
        WEIGHT = M0.WEIGHT;

    end

    save(fullfnsav,'DEC_SELCELL',...
        'CELLSEL_INX','SAMPLE_INX','DEC_M','dim_order',...
        'ORI_compindexset','WEIGHT',...
        'ORI_condset','ORI_list',...
        'comcont',...
        'mopts','topt0',...
        'getDfun',...
        'scriptname');

end % ises


if interparpool && ~isempty(parpoolid),
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


function cellinx = select_cell(W, CELLSEL_THRs, SEL_METHOD)

    [~, inxc]=sort(abs(W),'descend');

    [Ncell, Ncv] = size(W);
    

    cellinx = cell(length(CELLSEL_THRs),Ncv);
    switch lower(SEL_METHOD)
        case {'weightl1'}

        case {'ncell'}        

            for idx1 = 1: length(CELLSEL_THRs)    
                for icv = 1:Ncv
                    %---ncell==0, all cells are used
                    if CELLSEL_THRs(idx1)==0,
                        cidx = 1:Ncell;
                    else

                        cidx =inxc(1:min(CELLSEL_THRs(idx1),Ncell),icv);
                        cidx = cidx(:)';
                    end
                    cellinx{idx1,icv} = sort(cidx);
                end

            end
        otherwise
            error('not implemented');
    end
end



