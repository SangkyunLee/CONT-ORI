function get_subcell_gamnoise_crcon_decoding(iexp_type, DATA_thr_str, comcont, cellselopt, bshuffle,ctm,seslist,zmean,getDfun, parpoolid)
% function get_crcon_decoding(iexp_type, DATA_thr_str, comcont, bshuffle)
% test cross-contrast or contrast-independent models
% comcont={{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5_eyethr1';



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
         fntmp = {'GAM1_AN1-16','GAM1_AN17-22','','','GAM1'};
    case {'get_gamdat2'}
        fntmp = {'GAM2_AN1-16','GAM2_AN17-22','','','GAM2'};
    case {'get_gamdat3'}
        fntmp = {'GAM3_AN1-16','GAM3_AN17-22','','','GAM3'};
    case {'get_gamdat4'}
        fntmp = {'GAM4_AN1-16','GAM4_AN17-22','','','GAM4'};
    case {'get_gamdat5'}
        fntmp = {'GAM5_AN1-16','GAM5_AN17-22','','','GAM5'};
    case {'get_gamdat6'}
        fntmp = {'GAM6_AN1-16','GAM6_AN17-22','','','GAM6'};
    otherwise
end


CELLSEL_THRs = cellselopt.CELLSEL_THRs;
%SEL_METHOD = cellselopt.SEL_METHOD;
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
M2_lambda2s=0.01;%10.^(linspace(-10,1,20));


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
    


    fndata1 = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    [D1]= loadData(data_path,fndata1);

    
    
    Ncell = size(D1.cellinx_sel,2);
    if iexp_type==5
        D1.events_ORI(D1.events_ORI(:)==-15)=-10;
    end
    
      



    evt1 = [D1.events_cont(:) D1.events_ORI(:)];
    for icomp0 =  1:length(ORI_compindexset) % length(ORI_compindexset):-1:1 %
        
        icomp = ORI_compindexset(icomp0);      

        for icont = 1 : length(comcont)
            
            fprintf('ises: %d, icont:%d, icomp:%d\n',ises,icont, icomp);

            dtmp = DEC_SELCELL(1,icont,icomp0,ises);
            if all(dtmp(:))
                fprintf('ises: %d done.\n',ises);
                continue;
            end

%             if ~isempty(WEIGHT{icont,icomp0,ises}) && ~isempty(SAMPLE_INX{icont,icomp0,ises})
%                 fprintf('ises: %d done.\n',ises);
%                 continue;
%             end
            
            
            %------- prepare data
            contcond = comcont{icont};
            ORIsel = ORI_condset(:,icomp);            
            [dtr_gam,dte_gam,evttr,evtte] = get_gamdata(D1,evt1, dtype,contcond,ORIsel, getDfun); 

           
           
            if zmean
                dtr_gam = bsxfun(@minus, dtr_gam, mean(dtr_gam,2));                
            end
            %sc = sqrt(sum(dtr_gam.^2,2));
            %dtr_gam0 = bsxfun(@rdivide, dtr_gam,sc); % normalization
            dtr_gam = reshape(dtr_gam,[size(dtr_gam,1) 1 size(dtr_gam,2)]);  
            %dtr_gam0 = reshape(dtr_gam0,[size(dtr_gam0,1) 1 size(dtr_gam0,2)]);
            
            
            if ~isempty(dte_gam)
                if zmean
                    dte_gam = bsxfun(@minus, dte_gam, mean(dte_gam,2));                    
                end 
                %dte_gam0 = bsxfun(@rdivide, dte_gam,sc); % normalization
                dte_gam = reshape(dte_gam,[size(dte_gam,1) 1 size(dte_gam,2)]);                  
                %dte_gam0 = reshape(dte_gam0,[size(dte_gam0,1) 1 size(dte_gam0,2)]);
            end
           
           
           %-------------
            
            
            if any(comcont{icont}{1}==comcont{icont}{2}) && length(comcont{icont}{1})>1
                error('not implemented for contrast-independent decoder');
                % opts1 = mopts;
                % opts1.Nch =Ncell;
                % opts1.sel_cl = ORIsel';
                % common_inx = intersect(inx_train, inx_test);
                % inxmap1 = zeros(max(inx_train),1);                
                % inxmap1(inx_train,1)=1:length(inx_train);
                % inxmap2 = zeros(max(inx_test),1);                
                % inxmap2(inx_test,1)=1:length(inx_test);
                % common_inx_train = inxmap1(common_inx,1);
                % common_inx_test = inxmap2(common_inx,1);
                % opts1.common_sampleinx_train = common_inx_train;
                % opts1.common_sampleinx_test = common_inx_test;
                % 
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
%             if ~isempty(parpoolid)
%                 if all(contcond{1}==contcond{2})
%                     [out1 ] = par_cv_classification(dtr_gam0,evttr,opts1);
%                 else
%                     [out1 ] = par_cv_classification2(dtr_gam0,evttr,dte_gam0,evtte,opts1); 
%                 end
%             else
% 
%             end
%             DEC_M(icont,icomp0,ises) = mean(out1.acc_test);
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
            for ithr = 1% : length(CELLSEL_THRs)          

                
                
                %gen_outch = @(y)setdiff(1:opts2.Nch,y);
                %out_ch = cellfun(gen_outch,cellinx_TR_ses(ithr,:),'UniformOutput',false);
                out_ch = [];
                
                opts2.out_ch=out_ch;
                if ~isempty(parpoolid)
                    if all(comcont{icont}{1}==comcont{icont}{2})
                        [out2 ] = par_cv_classification(dtr_gam,evttr,opts2);
                    else
                        [out2 ] = par_cv_classification2(dtr_gam,evttr,dte_gam,evtte,opts2);
                    end
                else
%                     if all(comcont{icont}{1}==comcont{icont}{2})
%                         [out2 ] = cv_classification(dtr_gam,evttr,opts2);
%                     else
%                         [out2 ] = cv_classification2(dtr_gam,evttr,dte_gam,evtte,opts2);
%                     end
                end
                DEC_SELCELL(ithr,icont,icomp0,ises) = mean(out2.acc_test);
            end
            
            catch
                continue;
            end
            
            
            
            
            SAMPLE_INX{icont,icomp0,ises} = out2.inxs_cv;
            %WEIGHT{icont,icomp0,ises} = Ws;
            %CELLSEL_INX{icont,icomp0,ises} = cellinx_TR_ses;    
            
            mkdir(fileparts(fullfnsav));
            
            if exist(fullfnsav,'file')
               M0 = load(fullfnsav);
               pause(0.1);
                M0.DEC_SELCELL(:,icont,icomp0,ises)= DEC_SELCELL(:,icont,icomp0,ises);
                M0.CELLSEL_INX{icont,icomp0,ises} = CELLSEL_INX{icont,icomp0,ises};
                M0.DEC_M(icont,icomp0,ises) = DEC_M(icont,icomp0,ises);
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
                'opts1','opts2',...
                'getDfun',...
                'scriptname');
            pause(1);


        end % for icont  
    end %for icomp0 

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



