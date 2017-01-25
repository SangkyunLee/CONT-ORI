clear all
close all


ctm = 0.6;
ORIcomp = '0vs90'; %'-10vs0' %'-10vs0' %      %

% anestheisa
explist=[1 2 ];
%AW
% explist=[3];

ncont_plot=2;


exp_type={'','AN','','AWAKE'};
fntmp0 = {'','P2-AN17-22_','',''};
fntmp = {'','P1_P2-AN17-22_','',''};
DATA_thr_str = 'thr5';
contrast={'100', '40', '20'};
% 
% SEL_METHOD='WEIGHTL1'
% CELLSEL_THR=0:0.05:0.3;

CELLSEL_THR=0:5:30;
SEL_METHOD='NCELL';
iexp_type = 2;
seslist=17:22;
nthr = length(CELLSEL_THR);
SUBCELL_DEC = cell(1,3);
deaccALL = cell(1,3);
Ncell = zeros(1,22);
Ncell_sel = -100*ones(nthr,length(contrast),20,3);

clc


currentFolder = pwd;



if iexp_type ==2,
    ncont = 2;
    icomp0 = 2; % to get dec from 0vs90
    switch ORIcomp
        case {'0vs30'}
            icomp = 3;
        case {'0vs90'}
            icomp = 2;
        case {'-10vs0','-15vs0'}
            icomp = 6;
        otherwise
    end

elseif iexp_type ==5,

    ncont = 3;
    icomp0 = 2; % to get dec from 0vs90
    switch ORIcomp
        case {'0vs30'}
            icomp = 3;
        case {'0vs90'}
            icomp = 2;
        case {'-10vs0','-15vs0'}
            icomp = 6;
        otherwise
    end
end



fn1 = sprintf('%sSMLR_L1_ctm%0.2f.mat',fntmp0{iexp_type},ctm); %L1norm -regularization
fullfn1 = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fn1);
data0=load(fullfn1);


fndata = sprintf('%sCELLSEL_%s_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type},SEL_METHOD,ctm); %L2norm -regularization
fullfndata = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fndata);
data1=load(fullfndata);



for ises0 = 1 : length(seslist)
    ises = seslist(ises0);
    if isempty(data0.deaccALL{ises})
        continue;
    end            
    deaccALL(:,ises) = data1.DEC_SELCELL{icomp0,ises}(:,icont2,icont);


    for icont = 1:ncont
        for icont2 = 1: ncont
            SUBCELL_DEC(:,icont2,icont,ises) = data1.DEC_SELCELL{icomp,ises}(:,icont2,icont);
        end            
    end
    CELLSEL_L1DEC = data1.CELLSEL_L1DECORDER{icomp,ises};
    for icont=1:ncont
        for ithr = 1: length(data1.CELLSEL_THRs)
            Ncell_sel(ithr,icont,ises) = length(CELLSEL_L1DEC{ithr,icont});
        end
    end
    Ncell(ises) = length(data1.CELLSEL_INX{ises});


end


%% -------------- plot

M = Inf*ones(nthr,ncont_plot,ncont_plot,100);
Nc = Inf*ones(nthr,ncont_plot,100);
Nc0 = Inf*ones(1,100);

inx0=0;
for iexp_type=explist
    valid_ses = intersect(find(deaccALL{iexp_type}(2,:)>0.65 & deaccALL{iexp_type}(1,:)>0.65), find(Ncell(iexp_type,:)>30));
    if size(deaccALL{iexp_type},1)<ncont_plot
        continue;
    end
    M(:,:,:,inx0+1:inx0+length(valid_ses)) = SUBCELL_DEC{iexp_type}(:,1:ncont_plot,1:ncont_plot,valid_ses);
    
    Nc(1:nthr,1:ncont_plot,inx0+1:inx0+length(valid_ses)) = Ncell_sel(1:nthr,1:ncont_plot,valid_ses,iexp_type);
    Nc0(inx0+1:inx0+length(valid_ses)) = Ncell(iexp_type,valid_ses);
%     
    inx0 = inx0 + length(valid_ses); 
end
M = M(:,:,:,1:inx0);
Nc = Nc(:,:,1:inx0);
Nc0 = Nc0(1:inx0);




mM = mean(M,4)*100;
sM = std(M,0,4)/sqrt(inx0)*100;
sel_thrinx = [1 :2:length(CELLSEL_THR)];
linewidth=2;
CLR={'b','g','r'}
linestyle={'-','--',':'}
figure('Position',[100 100 300 230]); hold on;
for icont1= 1 : ncont_plot
    for icont2 = 1 : ncont_plot
        x = CELLSEL_THR(sel_thrinx);
        errorbar(x',mM(sel_thrinx,icont2,icont1),sM(sel_thrinx,icont2,icont1),'LineWidth',linewidth,'Color',CLR{icont1},'LineStyle',linestyle{icont2});
    end
end
if strcmp(SEL_METHOD,'WEIGHTL1')
    set(gca,'XTick',[0 0.1 0.2 0.3])
    tickv = num2cell([0 0.1 0.2 0.3]);
    tickstr = cellfun(@num2str,tickv, 'UniformOutput', false);
    set(gca,'XTickLabel',tickstr);
    set(gca,'FontSize',22);
    xlim([-0.05 0.35])
else
    set(gca,'XTick',[0 10 20 30])
    tickv = num2cell([0 10 20 30]);
    tickstr = cellfun(@num2str,tickv, 'UniformOutput', false);
    tickstr{1}='ALL';
    set(gca,'XTickLabel',tickstr);
    set(gca,'FontSize',22);
    xlim([-5 35])
end

% -------------------- plot relative change
    
    a = M(sel_thrinx,1,1,:)- M(sel_thrinx,2,1,:);
    a = squeeze(a)
    Ma = mean(a,2);
    Sa = std(a,0,2)/sqrt(size(a,2));
    figure; errorbar(Ma,Sa)
    
    
%----- plot ration of cells selected
if strcmp(SEL_METHOD,'WEIGHTL1')
    a=bsxfun(@rdivide,Nc,reshape(Nc0,[1 1 length(Nc0)]));
    figure('Position',[200 100 300 230]);hold on;
    for icont = 1:ncont_plot
    errorbar(x,mean(a(sel_thrinx,icont,:),3),std(a(sel_thrinx,icont,:),0,3)/sqrt(length(size(a,3))),'LineWidth',2,'Color',CLR{icont})       
    end
    set(gca,'XTick',[0 0.1 0.2 0.3])
    tickv = num2cell([0 0.1 0.2 0.3]);
    tickstr = cellfun(@num2str,tickv, 'UniformOutput', false);
    set(gca,'XTickLabel',tickstr);
    set(gca,'FontSize',22);
    xlim([-0.05 0.35])
    ylim([0 0.4]);
    box off
end
% 
% ------ stat
hs=[];
ps=[];
for icont = 1 : ncont_plot
    Mi = squeeze(M(sel_thrinx(2:end),:,icont,:));
    FWE = 0.05/size(Mi,1);
%     [h,p]=ttest(squeeze(Mi(1,1,:)), squeeze(Mi(1,2,:)),FWE);
    [h,p]=ttest(squeeze(Mi(:,1,:))'- squeeze(Mi(:,2,:))',[],FWE);
    hs(:,icont)=h;
    ps(:,icont)=p*size(Mi,1);
    
    Mi2 = reshape(Mi, [size(Mi,1)*size(Mi,2) size(Mi,3)])';
    [p t st]=anova1(Mi2);
    [c,m,h,nms] = multcompare(st,'alpha',0.05);
end


% 
% Mi = squeeze(M(2:7,:,icont,:));
% FWE = 0.05/size(Mi,1);
% [h,p]=ttest(squeeze(Mi(1,1,:)-Mi(6,2,:))',[],FWE);
% hs(:,icont)=h;
% ps(:,icont)=p*size(Mi,1);