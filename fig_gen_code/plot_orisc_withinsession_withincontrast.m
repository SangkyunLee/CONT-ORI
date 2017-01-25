% function sc=plot_orisc(thr,bdisp)
clear all
thr=0.5,bdisp=true;
% %- fit_ab
% M{1}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\AN1-16_ORIsc_ctm0.60.mat');
% M{2}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\AN17-22_ORIsc_ctm0.60.mat');
% M{3}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\AW23-40_ORIsc_ctm0.60.mat');
% 
% %- fit_ab3 J = || r40-r100*w||^2 -t*log(W)
% M{1}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\AN1-16_ORIsc_ctm0.60_fit_ab3.mat');
% M{2}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\AN17-22_ORIsc_ctm0.60_fit_ab3.mat');
% M{3}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\AW23-40_ORIsc_ctm0.60_fit_ab3.mat');

%- fit_ab4 J = || r40-r100*w||^2 +lambda|w| -t*log(W)
M{1}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_WSWC_AN1-16_ORIsc_ctm0.60_fit_ab4.mat');
M{2}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_WSWC_AN17-22_ORIsc_ctm0.60_fit_ab4.mat');
M{3}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\L_WSWC_AW23-40_ORIsc_ctm0.60_fit_ab4.mat');




%---------- plot tuning function
% thr = 0.5;
selinx = cell(1,3);
nses =0;
for iexp = 1:3
    a =cellfun(@isempty, M{iexp}.ORIsc(:,1));
    selinx{iexp} = find(a==0)';
    nses = nses + length(selinx{iexp});
end
clear sc    
sc(nses)=struct('a',[],'b',[],'rb',[],'cell',[],'badfit',[],'Mev',[],'iexp',[],'ises',[],'ncell',[]);
k = 1;
for iexp = 1:3      
    for ises = selinx{iexp}
        Mf1 = M{iexp}.ORIsc{ises,1}; % 1: for 40 vs 100, 2: for 20vs 100
        a = squeeze(Mf1.as(1,:,:));  % scale  C40 = a*C100 +b
        b = squeeze(Mf1.as(2,:,:));  % bias
        Mr = M{iexp}.ORIsc{ises,1}.Mresp;
        i100 = find(M{iexp}.ORIsc{ises,1}.cons==100,1,'first');
        lcon = M{iexp}.ORIsc{ises,1}.lcon;
        inx40 = find(lcon == 40);
        inx100 = find(lcon == 100);
        
        
        MMr = max(squeeze(Mr(:,i100,:)),[],1)';
        %Mr = mean(reshape(Mr,[size(Mr,1)*size(Mr,2) size(Mr,3)]),1);
        
        rb = bsxfun(@rdivide,b,MMr);
       
        %-- this is for fit_ab
%         Mf1.ev(a(:)<0)=0;
%         [Mev, mi ]= max(Mf1.ev,[],1);
%         mi = sub2ind(size(Mf1.ev),mi,1:size(Mf1.ev,2));
%         a = a(mi);
%         b = b(mi);
        
        %-- this is for fit_ab3
        Mev = Mf1.ev;
        a = a';
        b = b';
        rb = rb';
        %----------


        CLIST = find(Mev>thr );
        nC= length(CLIST);
        sc(k).a = a(CLIST);
        sc(k).b = b(CLIST);
        sc(k).rb = rb(CLIST);
        sc(k).badfit = abs(sc(k).rb)>1;      
        sc(k).Mev = Mev(CLIST);
        sc(k).mresp = Mr(:,[inx40 inx100],CLIST);        
        sc(k).cell = CLIST;
        sc(k).iexp = iexp;
        sc(k).ises =ises;
        sc(k).ncell = length(a);
        k = k+1;
    end
end

a=[]; b = []; rb =[]; gd=[]; Mev=[]; ncell=[]; Mresp=NaN*ones(4,2,1);
for i = 1 : nses
    sz = size(sc(i).mresp);
%     Mresp(1:sz(1),1:sz(2),length(a)+1:length(a)+sz(3)) = sc(i).mresp;
    
    a = [a sc(i).a];
    b = [b sc(i).b];
    rb = [rb sc(i).rb];
    gd = [gd length(find(~sc(i).badfit))];    
    Mev = [Mev sc(i).Mev];
    
    ncell = [ncell sc(i).ncell];
end

    

if bdisp
    %------------- plot scaling factor
    figure; hist(a,-3:0.1:4);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    set(gca,'FontSize',20,'XTick',[0 1 2]);
    xlim([-0.1 4]); box off

        figure; hist(rb,-3:0.1:3)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    set(gca,'FontSize',20,'XTick',[-1 0 1]);
    set(gca,'FontSize',20,'YTick',[0:100:300]);
    xlim([-1 1])
    box off



    % ----- relationship 
    figure;
    plot(a,rb,'.'); 
    hold on; plot([0 3],[0 0],'k--')
    hold on; plot([1 1],[-1 1],'k--')
    xlim([-0.1 5]); box off
    set(gca,'FontSize',20,'XTick',[0 1 2]);
    set(gca,'YTick',[-0.5 0 0.5]);

end


%% ---- get number of cells used in the analysis 
ncell=zeros(28,1);
j=1;
for iexp = 1:3
    nsub = size(M{iexp}.ORIsc,1)
    for isub = 1 : nsub
        K= M{iexp}.ORIsc{isub,1};
        if ~isempty(K)
            ncell(j) = length(K.ev);
            j = j+1;
        end
    
    end
end

