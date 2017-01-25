function get_ORI_conscale3(iexp_type,DATA_thr_str,ctm,calfunstr,parpoolid)
% data fitting within contrast within session


%----------------
exp_type={'AN','AN','AN_0TO150','AWAKE','AWAKE_EYE'};
datainxstr = {'AN1-16','AN17-22','','','AW23-40'};


cell_sel_method = 'UNION_CONTRSP'; 


pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
fnsave = sprintf('L_WSWC_%s_ORIsc_ctm%0.2f_%s.mat',datainxstr{iexp_type},ctm,calfunstr); 

if ~exist('parpoolid','var')
    parpoolid = set_env(true);
end
[contrasts, ~, ~, ~, seslist] =get_expinfo(iexp_type);

% bnormaldata = true;


currentFolder = pwd;
fullfnsav = fullfile(fileparts(currentFolder),'GRP_data',exp_type{iexp_type},DATA_thr_str, fnsave);


ORIsc =cell(max(seslist),length(contrasts));
ORItun(max(seslist))=struct;


%% ---------------------------------------------------
for ises = seslist     
    fndata1 = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    [data]= loadData(data_path,fndata1);
    
    X0=data.Xsel;
    if mean(X0(:))<1,
        X0 = X0*100;
    end
%     if bnormaldata
%         scale = sqrt(sum(X0.^2,1));
%         X0=bsxfun(@rdivide,X0,scale);
%     end

    [inxsample, cons, oris] = collect(data);

    %------- est scale and bias between ori tunings within the same
    %contrast
    evt_cond.inxsample = inxsample;
    evt_cond.cons = cons;
    evt_cond.oris = oris;
    evt_cond.conref = 100;
    evt_cond.concom = 100;
    ORIsc{ises,1} = cal_wori(X0,evt_cond,calfunstr);
    
    evt_cond.inxsample = inxsample;
    evt_cond.cons = cons;
    evt_cond.oris = oris;
    evt_cond.conref = 40;
    evt_cond.concom = 40;
    ORIsc{ises,2} = cal_wori(X0,evt_cond,calfunstr);

    
    
    
    if find(contrasts==20)
        evt_cond.conref = 20;
        evt_cond.concom = 20;
        ORIsc{ises,3}= cal_wori(X0,evt_cond,calfunstr);
    end

    %---- est tuning curves
    nevt = length(cons);
    mX=zeros(nevt,size(X0,2));
    sX1=zeros(nevt,size(X0,2));
    sX2=zeros(nevt,size(X0,2));
    for ievt = 1:nevt
        mX(ievt,:) = mean(X0(inxsample{ievt},:),1);
        sX2(ievt,:) = std(X0(inxsample{ievt},:),0,1)/sqrt(length(inxsample{ievt}));
        sX1(ievt,:) = std(X0(inxsample{ievt},:),0,1);
    end
    evtlist = [cons' oris'];
    orders = {'descend','ascend'};
    [Out, J0] = TP.sort_evtorder(evtlist,orders);

    mresp = reshape(mX(J0,:),[length(unique(oris)) length(unique(cons)) size(mX,2)]);
    sresp2 = reshape(sX2(J0,:),[length(unique(oris)) length(unique(cons)) size(sX2,2)]);
    sresp1 = reshape(sX1(J0,:),[length(unique(oris)) length(unique(cons)) size(sX1,2)]);
    evtord = zeros(length(unique(oris)),length(unique(cons)), size(Out,2)); 
    for ievt = 1 : size(Out,2)
        evtord(:,:,ievt) = reshape(Out(:,ievt),[length(unique(oris)) length(unique(cons))]);
    end
    ORItun(ises).mresp = mresp;
    ORItun(ises).semresp = sresp2;
    ORItun(ises).stdresp = sresp1;
    ORItun(ises).evtord = evtord;
    
end % ises
scn = mfilename('fullpath');
save(fullfnsav,'ORIsc','ORItun','scn','-v7.3') ;

% if ~isempty(parpoolid),
%     delete(parpoolid)
%     pause(3);
% end

end


    %--------- collecting data
function [inxsample, cons, oris] = collect(data)


        inxs_valtrial = data.events(:)>0 & ~isinf(data.events(:));
        unique_evt = unique(data.events(inxs_valtrial)');
        unique_evt = setdiff(unique_evt,0);
        nevt = length(unique_evt);
        evts = data.events(:);
        inxsample = cell(1,nevt);
        cons = zeros(1,nevt);
        oris = zeros(1,nevt);
        for ievt0 = 1:nevt         
            ievt = unique_evt(ievt0);
            inxsample{ievt0} = TP.select_subdata(evts,{ievt});
            cons(ievt0)=data.events_cont(inxsample{ievt0}(1));
            oris(ievt0)=data.events_ORI(inxsample{ievt0}(1));
        end
end


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


% 
% %-------------------------------
% function [D1,D2] = subdata(subinx,D1,D2,dtyps,dr)
% % function [D1,D2] = subdata(M,D1,D2,dtyps)
% % retrieve subdata with subinx in given dtyps
% % NOTE: it does not retrieve all data automatically, 
% % it retrive data only with dyps (cell type variable)
% % dr specifies dimension inx of Dx.(dyps{ityp}) 
%     
%     if nargin<5
%         dr = 2;
%     end
%       
%     inx1 = D1.cellinx_sel;
%     inx2 = D2.cellinx_sel;
% 
%     [~,~,ni1]=intersect(subinx,inx1);
%     [~,~,ni2]=intersect(subinx,inx2);
% 
%     D1.cellinx_sel =D1.cellinx_sel(ni1);
%     D2.cellinx_sel =D2.cellinx_sel(ni2);
%     for ityp = 1 : length(dtyps)
%         dtyp = dtyps{ityp};
%         if dr==1,
%             D1.(dtyp) = D1.(dtyp)(ni1,:);
%             D2.(dtyp) = D2.(dtyp)(ni2,:);
%         elseif dr==2,
%             D1.(dtyp) = D1.(dtyp)(:,ni1);
%             D2.(dtyp) = D2.(dtyp)(:,ni2);
%         end
%     end
% end