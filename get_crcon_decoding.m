function get_crcon_decoding(iexp_type, DATA_thr_str, comcont, bshuffle)
% function get_crcon_decoding(iexp_type, DATA_thr_str, comcont, bshuffle)
% test cross-contrast or contrast-independent models
% comcont={{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5_eyethr1';



if nargin<4
    bshuffle = false;
end

parpoolid = set_env(true);
[contrasts, ORI_list, ORI_compindexset, ~, seslist] =get_expinfo(iexp_type);

%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};

%----------------

ctm=0.6; 


cell_sel_method = 'UNION_CONTRSP'; 
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);


currentFolder = pwd;
if bshuffle
    fnsave = sprintf('SHUFFLE_ACELL_CRSCON_SMLR_L2_ctm%0.2f.mat',ctm); %L2norm -regularization
else
    fnsave = sprintf('ACELL_CRSCON_SMLR_L2_ctm%0.2f.mat',ctm); %L2norm -regularization
end
fullfnsav = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fnsave);
lambda1s=10.^(linspace(-10,1,20));
lambda2s=0;







ORI_condset=combnk(ORI_list,2)';


DEC_SELCELL = zeros(length(ORI_compindexset),length(comcont),max(seslist));

CELLSEL_INX = cell(max(seslist),1);



%% ---------------------------------------------------
for ises = seslist 
    if exist(fullfnsav,'file')
        load(fullfnsav);
    end
     if ~isempty(CELLSEL_INX{ises})
        continue;
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

            
            
            
            
            
            
            %---------- select train samples ----------------
            sevts{1} = comcont{icont}{1};
            sevts{2} = ORI_condset(:,icomp);
            inxsample = TP.select_subdata(evt1,sevts);
            % similar number of samples between two conditions,not exceed 10% more samples               
            inx2 = TP.select_samples_evtratio(D1.events_ORI(inxsample),sevts{2}, 1.1);
            inx_train = inxsample(inx2);
            
            %----------------------------------
      
            copts.bshuffle = bshuffle;
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
    
            
                         
                    
            % excluding cell based on training data contrast
            sevts{1} = comcont{icont}{2};
            sevts{2} = ORI_condset(:,icomp);
            inxsample = TP.select_subdata(evt1,sevts);                        
            inx2 = TP.select_samples_evtratio(D1.events_ORI(inxsample),sevts{2}, 1.1);
            inx_test = inxsample(inx2);

            data_test = D1.Xsel(inx_test,:)';
            evtsel_test = D1.events_ORI(inx_test);
            data_test = reshape(data_test,[size(data_test,1) 1 size(data_test,2)]); 
            opts1 = copts;
            
            
            common_inx = intersect(inx_train, inx_test);
            inxmap1 = zeros(max(inx_train),1);                
            inxmap1(inx_train,1)=1:length(inx_train);
            inxmap2 = zeros(max(inx_test),1);                
            inxmap2(inx_test,1)=1:length(inx_test);
            common_inx_train = inxmap1(common_inx,1);
            common_inx_test = inxmap2(common_inx,1);
            
            
            
            opts1.common_sampleinx_train = common_inx_train;
            opts1.common_sampleinx_test = common_inx_test;
            opts1.perval = min(0.1*length(evtsel_test)/length(evtsel_train),0.1);
            opts1.pertrain = min(0.8*length(evtsel_test)/length(evtsel_train),0.8);

            % excluding cell based on training data contrast
            opts1.out_ch=[];
            if ~isempty(parpoolid)
                [out1 ] = par_cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts1); 
            else
                %[out1 ] = cv_classification2(data_train,evtsel_train,data_test,evtsel_test,opts1); 
            end
            
            DEC_SELCELL(icomp0,icont,ises) = mean(out1.acc_test);
            
                    
                


        end % for icont  

       

    end %for icomp0 
    
    copts.sel_cl=[];

    scriptname = mfilename('fullpath');

    mkdir(fileparts(fullfnsav));
    save(fullfnsav,'DEC_SELCELL',...
        'ORI_compindexset',...
        'ORI_condset','ORI_list',...
        'comcont',...
        'copts',...
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

