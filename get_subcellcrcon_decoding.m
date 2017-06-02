function get_subcellcrcon_decoding(iexp_type, DATA_thr_str, comcont, cellselopt, bshuffle,ctm,seslist,zmean)
% function get_crcon_decoding(iexp_type, DATA_thr_str, comcont, bshuffle)
% test cross-contrast or contrast-independent models
% comcont={{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5_eyethr1';



if nargin<5
    bshuffle = false;
end

parpoolid = set_env(true);
[~, ORI_list, ORI_compindexset, nses, ses0] =get_expinfo(iexp_type);
if ~exist('seslist','var') || isempty(seslist)
    seslist = ses0;
end


%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
fntmp = {'AN1-16','AN17-22','','',''};
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
    fnsave = sprintf('SHUFFLE_%sSUBCELL-%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type}, dtype,ctm); %L2norm -regularization
else
    fnsave = sprintf('%sSUBCELL-%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type},dtype,ctm); %L2norm -regularization
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
opts.Ncv = 100;    
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


% common options to be used for testing decorders in selected cells;
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
    for icomp0 = 1 : length(ORI_compindexset)
        
        icomp = ORI_compindexset(icomp0);      

        for icont = 1 : length(comcont)
            
            fprintf('ises: %d, icont:%d, icomp:%d\n',ises,icont, icomp);

            if DEC_SELCELL(1,icont,icomp0,ises)~=0
                fprintf('ises: %d done.\n',ises);
                continue;
            end

%             if ~isempty(WEIGHT{icont,icomp0,ises}) && ~isempty(SAMPLE_INX{icont,icomp0,ises})
%                 fprintf('ises: %d done.\n',ises);
%                 continue;
%             end
            
            
            
            %---------- select train samples ----------------
            sevts{1} = comcont{icont}{1};
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
            sc = sqrt(sum(data_train.^2,2));
            data_train0 = bsxfun(@rdivide, data_train,sc); % normalization
            evtsel_train = D1.events_ORI(inx_train);             
            data_train = reshape(data_train,[size(data_train,1) 1 size(data_train,2)]);  
            data_train0 = reshape(data_train0,[size(data_train0,1) 1 size(data_train0,2)]);  
    
            
                         
            if all(comcont{icont}{1}==comcont{icont}{2})        
                inx_test=[];
            else
                sevts{1} = comcont{icont}{2};
                sevts{2} = ORI_condset(:,icomp);
                inxsample = TP.select_subdata(evt1,sevts);                        
                inx2 = TP.select_samples_evtratio(D1.events_ORI(inxsample),sevts{2}, 1.1);
                inx_test = inxsample(inx2);

                data_test = D1.(dtype)(inx_test,:)';         
                
                if zmean
                    data_test = bsxfun(@minus, data_test, mean(data_test,2));                    
                end 
                data_test0 = bsxfun(@rdivide, data_test,sc); % normalization
                
                evtsel_test = D1.events_ORI(inx_test);
                data_test = reshape(data_test,[size(data_test,1) 1 size(data_test,2)]);                  
                data_test0 = reshape(data_test0,[size(data_test0,1) 1 size(data_test0,2)]); 
            end
            
            
            common_inx = intersect(inx_train, inx_test);
            inxmap1 = zeros(max(inx_train),1);                
            inxmap1(inx_train,1)=1:length(inx_train);
            inxmap2 = zeros(max(inx_test),1);                
            inxmap2(inx_test,1)=1:length(inx_test);
            common_inx_train = inxmap1(common_inx,1);
            common_inx_test = inxmap2(common_inx,1);
            
            opts1 = mopts;
            opts1.Nch =Ncell;
            opts1.sel_cl = ORI_condset(:,icomp)';
    
            opts1.common_sampleinx_train = common_inx_train;
            opts1.common_sampleinx_test = common_inx_test;
            
            if ~all(comcont{icont}{1}==comcont{icont}{2})
                opts1.perval = min(opts1.perval*length(evtsel_test)/length(evtsel_train),opts1.perval);
                opts1.pertrain = min(opts1.pertrain*length(evtsel_test)/length(evtsel_train),opts1.pertrain);            
            end
            opts1.out_ch=[];
            
%             %-------- temporal code
% %             opts1.inxs_cv = SAMPLE_INX{icont,icomp0,ises};
% %             opts1.common_sampleinx_train = [];
% %             opts1.common_sampleinx_test = [];
% %             opts1.cvmode = 'precal';
%             %----------------------

            
            % selecting cell with L1 - smlr
            try
            if ~isempty(parpoolid)
                if all(comcont{icont}{1}==comcont{icont}{2})
                    [out1 ] = par_cv_classification(data_train0,evtsel_train,opts1);
                else
                    [out1 ] = par_cv_classification2(data_train0,evtsel_train,data_test0,evtsel_test,opts1); 
                end
            else
%                 if all(comcont{icont}{1}==comcont{icont}{2})
%                     [out1 ] = cv_classification(data_train0,evtsel_train,opts1);
%                 else
%                     [out1 ] = cv_classification2(data_train0,evtsel_train,data_test0,evtsel_test,opts1);
%                 end
            end
            DEC_M(icont,icomp0,ises) = mean(out1.acc_test);
            
            
            
            Ws = squeeze(out1.Ws(1:Ncell,1,:));
            cellinx_TR_ses = select_cell(Ws, CELLSEL_THRs, SEL_METHOD);
            
            %-------------
            opts2 = topt0;
            opts2.Nch =Ncell;
            opts2.sel_cl = ORI_condset(:,icomp)';
            opts2.cvmode = 'precal';
            opts2.pertest = [];
            opts2.pertrain = [];
            opts2.perval = [];
            opts2.inxs_cv = out1.inxs_cv;
            opts2.common_sampleinx_train = [];
            opts2.common_sampleinx_test = [];
            %----------------
            for ithr = 1 : length(CELLSEL_THRs)          

                
                
                gen_outch = @(y)setdiff(1:opts2.Nch,y);
                out_ch = cellfun(gen_outch,cellinx_TR_ses(ithr,:),'UniformOutput',false);
                
                
                opts2.out_ch=out_ch;
                if ~isempty(parpoolid)
                    if all(comcont{icont}{1}==comcont{icont}{2})
                        [out2 ] = par_cv_classification(data_train,evtsel_train,opts2);
                    else
                        [out2 ] = par_cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts2);
                    end
                else
%                     if all(comcont{icont}{1}==comcont{icont}{2})
%                         [out2 ] = cv_classification(data_train,evtsel_train,opts2);
%                     else
%                         [out2 ] = cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts2);
%                     end
                end
                DEC_SELCELL(ithr,icont,icomp0,ises) = mean(out2.acc_test);
            end
            
            catch
                continue;
            end
            
            
            
            WEIGHT{icont,icomp0,ises} = Ws;
            SAMPLE_INX{icont,icomp0,ises} = opts2.inxs_cv;
            CELLSEL_INX{icont,icomp0,ises} = cellinx_TR_ses;    
            
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
                'scriptname');
            pause(1);


        end % for icont  
    end %for icomp0 

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



