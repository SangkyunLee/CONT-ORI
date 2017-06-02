% plot eye statistics
iexp_type = 5;
ctm=0.6
exp_type={'','','','','AWAKE_EYE'};


DATA_thr_str = 'thr5_eyethr_xy1_p1';
cell_sel_method = 'UNION_CONTRSP'; 
pprotype=['DATA_DISK_' cell_sel_method];
data_path = fullfile('../GRP_data/', exp_type{iexp_type},DATA_thr_str);

[~, ~, ~, ~,ses] =get_expinfo(iexp_type);

N = zeros(length(ses),1);
k=1;
Ps=zeros(length(ses),3);
PALL =cell(length(ses),2);
XALL = PALL;
YALL = PALL;
for ises = ses
    fndata1 = sprintf('%s_ctm%0.2fses%d.mat',pprotype,ctm,ises);
    [D1]= loadData(data_path,fndata1);
    
    scanlist =D1.scanlist;
    ns = length(scanlist);
    cont = cell(ns,1); PR = cell(ns,1); XY = cell(ns,1);
    for is = 1 : ns
        trinx = D1.trialinxs{is};
        eye = D1.eye(scanlist(is));
        cont0 = D1.events_cont(:,is);
        cont{is} = cont0(~isinf(cont0));
        PR{is} = cellfun(@mean,eye.r(trinx));
        XY{is} = cell2mat(cellfun(@mean,eye.xy(trinx),'UniformOutput',false));
    end
    cont = cell2mat(cont);
    PR = cell2mat(PR);
    XY = cell2mat(XY);
    inx1 = cont ==100;
    inx2 = cont ==40;
    p1 = ranksum(PR(inx1), PR(inx2));
    p21 = ranksum(XY(inx1,1), XY(inx2,1));
    p22 = ranksum(XY(inx1,2), XY(inx2,2));
    if ises==32,
        u = 3.4/410; % I assume the eye horizontal 3.4mm and the pixel : 410
        %pupil
        Y1 = PR*u; binP=0.1:0.05:0.5;
        %x,y
        Y2 = XY*u; binX=-0.5:0.1:0.5;
        
        Y2(1:659,:) = bsxfun(@minus, Y2(1:659,:),median(Y2(1:659,:),1));
        Y2(660:end,:) = bsxfun(@minus, Y2(660:end,:),median(Y2(660:end,:),1));
        Y=[Y1 Y2];
        for ii=1:3
            figure; hold on;
            if ii==1, bin = binP; 
            else bin =binX; end
                
            h1 = hist(Y(inx1,ii),bin);
            h2 = hist(Y(inx2,ii),bin);
            plot(bin,h1,'.-r','linewidth',2,'markersize',20);                
            plot(bin,h2,'.-b','linewidth',2,'markersize',20);                                        
            set(gca,'FontSize',22,'linewidth',2); box off
        end
        r=10.8/0.1; % convert mm to deg
        Y2deg = sqrt(sum(Y2.^2,2))*r;
    end
    PALL{k,1}=PR(inx1)/max(PR);
    PALL{k,2}=PR(inx2)/max(PR);
    XY0 = bsxfun(@minus, XY,mean(XY,1));
    XALL{k,1}= XY0(inx1,1)/max(abs(XY0(:,1)));
    XALL{k,2}= XY0(inx2,1)/max(abs(XY0(:,1)));
    YALL{k,1}= XY0(inx1,2)/max(abs(XY0(:,2)));
    YALL{k,2}= XY0(inx2,2)/max(abs(XY0(:,2)));
    
    
    N(k)=length(PR);
    Ps(k,:)=[p1 p21 p22];
    k = k +1;
end
Ps=[ses' Ps];


bin=0:0.05:1;
M1 = cell2mat(PALL(:,1));
h1 = hist(M1,bin);
M2 = cell2mat(PALL(:,2));
h2 = hist(M2,bin);
figure; hold on;
plot(bin,h1,'.-r','linewidth',2,'markersize',20);                
plot(bin,h2,'.-b','linewidth',2,'markersize',20);                                        
set(gca,'FontSize',22,'linewidth',2); box off

ranksum(M1,M2)


