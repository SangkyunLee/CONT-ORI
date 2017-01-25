clear all
exp_type={'AN','AN_0TO150','AWAKE','AN_FULLORI'};

%----------------
iexp_type =4;
ctm = 0.6; 
if iexp_type==1,
    seslist = [(1:8) (11:16)]
elseif iexp_type==3
    seslist= 22; % awake
elseif iexp_type==4
    seslist=2
else
    seslist = 1:7;
end

%-------------------

cell_sel_method = 'UNION_CONTRSP'; 
DATA_thr_str = 'thr5';
%pprotype=['NEW1_DATA_RIM_' cell_sel_method];

pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);
hfig1= figure;
hfig2=figure;
for ises1 = 1:length(seslist)
    ises = seslist(ises1);
    fndata = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    fullfndata = fullfile(data_path,fndata);

    data=load(fullfndata);
    X = data.Xsel;
    scale=1./sqrt(sum(X.^2,1));
    X = bsxfun(@times, X, scale);

    inxvalid = ~isinf(data.events(:));

    events = data.events(inxvalid);
    evtlist = unique(events)';
    evt_order =zeros(2,length(evtlist));

    mX = zeros(max(evtlist),size(data.Xsel,2));
    stdX = zeros(max(evtlist),size(data.Xsel,2));
    for ievt = evtlist
        inxtrial = find(events==ievt);
        datasub = X(inxtrial,:);
        mX(ievt,:) = mean(datasub,1);
        stdX(ievt,:) = std(datasub,0,1);
        evt_order(1,ievt)=data.events_ORI(find(events==ievt,1));
        evt_order(2,ievt)=data.events_cont(find(events==ievt,1));
    end
    
    inx_cont100 = find(evt_order(2,:)==100);
    inx_cont40 = find(evt_order(2,:)==40);
    
    cont100 = mX(inx_cont100,:)./stdX(inx_cont100,:);
    cont40 = mX(inx_cont40,:)./stdX(inx_cont40,:);
    figure(hfig1);
    if iexp_type==3,
        subplot(4,5,ises1);
    else
        subplot(4,4,ises1);
    end
    plot_linearreg([cont100(:),cont40(:)]);
    
    nORI = length(inx_cont100);
    mX100 =mX(inx_cont100,:);
    stdX100=stdX(inx_cont100,:);
    [~,inx] = max(mX100,[],1);
    inx = inx + 0 :nORI : (size(mX100,2)-1)*nORI;
    cont100 = mX100(inx)./stdX100(inx);
    
    mX40 =mX(inx_cont40,:);
    stdX40=stdX(inx_cont40,:);
    [~,inx] = max(mX40,[],1);
    inx = inx + 0 :nORI : (size(mX40,2)-1)*nORI;
    cont40 = mX40(inx)./stdX40(inx);
    figure(hfig2);
    if iexp_type==3,
        subplot(4,5,ises1);
    else
        subplot(4,4,ises1);
    end
    plot_linearreg([cont100(:),cont40(:)]);
    
    
    axis equal
    xlabel('Cont100 SNR');
    ylabel('Cont40 SNR');
end
    