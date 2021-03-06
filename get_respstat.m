function get_respstat(iexp_type,DATA_thr_str)
% iexp_type=1
% DATA_thr_str = 'thr5';



[contrasts, ORI_list, ~, nses, seslist] =get_expinfo(iexp_type);

%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
fntmp = {'AN1-16','AN17-22','','','AW23-40'};
%----------------

ctm=0.6; 


cell_sel_method = 'UNION_CONTRSP'; 
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);

% 
currentFolder = pwd;

fnsave = sprintf('BASIC_SUMMARY_%s_ctm%0.2f.mat',fntmp{iexp_type},ctm); %L2norm -regularization

fullfnsav = fullfile(fileparts(currentFolder),'GRP_data',exp_type{iexp_type},DATA_thr_str, fnsave);

% 



S.cellist=[];
S.cr1dinx=[];
S.Ncell =0;
S.V=[];
S.NR=[];
S.m = [];
S.SR=[];
S.Raw=[];
S = repmat(S,[1 nses]);

%% ---------------------------------------------------
for ises = seslist 



    fndata1 = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    [D1]= loadData(data_path,fndata1);

    
    
    nc = size(D1.cellinx_sel,2);
    if iexp_type==5
        D1.events_ORI(D1.events_ORI(:)==-15)=-10;
    end
    
      



    evt1 = [D1.events_cont(:) D1.events_ORI(:)];
    raw = cell(length(ORI_list),length(contrasts)); % collect all trials
    V = zeros(nc,length(ORI_list),length(contrasts)); % variance
    NR = zeros(nc*(nc-1)/2,length(ORI_list),length(contrasts)); %noise correlation
    m = zeros(nc,length(ORI_list),length(contrasts)); % mean
    SR = zeros(nc*(nc-1)/2,1,length(contrasts)); %noise correlation

    for icont = 1 : length(contrasts)
        for icomp = 1 : length(ORI_list)
            
            fprintf('ises: %d, icont:%d, icomp:%d\n',ises,icont, icomp);
            
            %---------- select samples ----------------
            sevts{1} = contrasts(icont);
            sevts{2} = ORI_list(icomp);
            inxsample = TP.select_subdata(evt1,sevts);
            D = D1.Xsel(inxsample,:);
            
            V(:,icomp,icont) = var(D,0,1);
            y1d = get_paircorr(D,1);
            NR(:,icomp,icont) = y1d;
            m(:,icomp,icont) = mean(D,1);
            raw{icomp,icont} =D;
        end %for icomp0 
        [y1d,inx] = get_paircorr(m(:,:,icont),2);
        SR(:,1,icont) = y1d;
    end % for icont  
    
    S(ises).cellist = D1.cellinx_sel;
    S(ises).cr1dinx = inx;
    S(ises).Ncell =nc;
    S(ises).V=V;
    S(ises).NR=NR;
    S(ises).m = m;
    S(ises).SR=SR;
    S(ises).Raw =raw;
end % ises

scriptname = mfilename('fullpath');

mkdir(fileparts(fullfnsav));
save(fullfnsav,'S',...    
    'contrasts','ORI_list',...    
    'scriptname');
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

