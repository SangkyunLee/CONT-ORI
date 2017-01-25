clear all
close all


ctm = 0.6;
ORIcomp ='-10vs0' %'-10vs0' %  '-10vs0'%    % '0vs90' %

% anestheisa
% explist=[1 2];
%AW
explist=[1:2];

ncont_plot=2;


exp_type={'AN','AN','AN_0TO150','AWAKE'};
fntmp = {'AN1-16_','AN17-22_','',''};
DATA_thr_str = 'thr5';
contrast={'100', '40', '20'};

% SEL_METHOD='WEIGHTL1'
% CELLSEL_THR=0:0.05:0.3;

CELLSEL_THR=0:5:30;
SEL_METHOD='NCELL';

nthr = length(CELLSEL_THR);
SUBCELL_DEC = cell(2,3);
deaccALL = cell(1,3);
Ncell = zeros(3,20);
Ncell_sel = -100*ones(nthr,length(contrast),20,3);

clc


currentFolder = pwd;

for iexp_type = 1:2
    fn1 = sprintf('%sSMLR_L1_ctm%0.2f.mat',fntmp{iexp_type},ctm); %L1norm -regularization
    %fn1 = sprintf('SMLR_L1_ctm%0.2f.mat',ctm); %L1norm -regularization
    fullfn1 = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fn1);
    data0=load(fullfn1);

    fndata = sprintf('%sCELLSEL_%s_CRSCON_SMLR_L2_ctm%0.2f.mat',fntmp{iexp_type},SEL_METHOD,ctm); %L2norm -regularization
    %fndata = sprintf('CELLSEL_%s_CRSCON_SMLR_L2_ctm%0.2f.mat',SEL_METHOD,ctm); %L2norm -regularization
    fullfndata = fullfile(fileparts(currentFolder),'NEW_DECODING',exp_type{iexp_type},DATA_thr_str, fndata);
    data1=load(fullfndata);
    
    
    if iexp_type==2,
        ncont = 2;
        icomp0 = 2; % to get dec from 0vs90
        switch ORIcomp
            case {'0vs30'}
                icomp = 3;
            case {'0vs90'}
                icomp = 2;
            case {'-10vs0'}
                icomp = 6;
            otherwise
        end
    elseif iexp_type ==3,
        ncont = 2; % # contrast
        icomp0 = 4; % to get dec from 0vs90
        switch ORIcomp
            case {'0vs30'}
                icomp = 1;
            case {'0vs90'}
                icomp = 4;
            otherwise
        end

    else
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






    
    for ises = 1 : length(data1.CELLSEL_INX)
        if isempty(data0.deaccALL{ises})
            %deaccALL{iexp_type}(:,ises)=0;
            continue;
        end            
        deaccALL{iexp_type}(:,ises) = data0.deaccALL{ises}(:,icomp0);

        
        for icont = 1:ncont
            for icont2 = 1: ncont
                SUBCELL_DEC{1,iexp_type}(:,icont2,icont,ises) = data1.DEC_SELCELL{icomp,ises}(:,icont2,icont,1);
                SUBCELL_DEC{2,iexp_type}(:,icont2,icont,ises) = data1.DEC_SELCELL{icomp,ises}(:,icont2,icont,2);                
            end            
        end
        CELLSEL_L1DEC = data1.CELLSEL_L1DECORDER{icomp,ises};
        for icont=1:ncont
            for ithr = 1: length(data1.CELLSEL_THRs)
                Ncell_sel(ithr,icont,ises,iexp_type) = length(CELLSEL_L1DEC{ithr,icont});
            end
        end
        Ncell(iexp_type,ises) = length(data1.CELLSEL_INX{ises});
   
        
    end

end
%% -------------- plot

Mtr = Inf*ones(nthr,ncont_plot,ncont_plot,100);
Mte = Inf*ones(nthr,ncont_plot,ncont_plot,100);
Nc = Inf*ones(nthr,ncont_plot,100);
Nc0 = Inf*ones(1,100);

inx0=0;
for iexp_type=explist
    valid_ses = intersect(find(deaccALL{iexp_type}(2,:)>0.7 & deaccALL{iexp_type}(1,:)>0.7), find(Ncell(iexp_type,:)>29));
    if size(deaccALL{iexp_type},1)<ncont_plot
        continue;
    end
    Mtr(:,:,:,inx0+1:inx0+length(valid_ses)) = SUBCELL_DEC{1,iexp_type}(:,1:ncont_plot,1:ncont_plot,valid_ses);
    Mte(:,:,:,inx0+1:inx0+length(valid_ses)) = SUBCELL_DEC{2,iexp_type}(:,1:ncont_plot,1:ncont_plot,valid_ses);
    Nc(1:nthr,1:ncont_plot,inx0+1:inx0+length(valid_ses)) = Ncell_sel(1:nthr,1:ncont_plot,valid_ses,iexp_type);
    Nc0(inx0+1:inx0+length(valid_ses)) = Ncell(iexp_type,valid_ses);
%     
    inx0 = inx0 + length(valid_ses); 
end


Mtr = Mtr(:,:,:,1:inx0);
Mte = Mte(:,:,:,1:inx0);


Nc = Nc(:,:,1:inx0);
Nc0 = Nc0(1:inx0);


mMtr = mean(Mtr,4)*100;
sMtr = std(Mtr,0,4)/sqrt(inx0)*100;
mMte = mean(Mte,4)*100;
sMte = std(Mte,0,4)/sqrt(inx0)*100;

for icont2 = 1 : ncont_plot
    grpstr =cell(1,ncont_plot);
    
    inxk=1;
    mM=zeros(size(mMtr,1),1+(ncont_plot-1)*2);
    sM=zeros(size(mMtr,1),1+(ncont_plot-1)*2);
    
    mM(:,inxk) = mMtr(:,icont2,icont2);
    sM(:,inxk) = sMtr(:,icont2,icont2);
    M(:,inxk,:,icont2)=Mtr(:,icont2,icont2,:);
    grpstr{inxk}= sprintf('%s-%s',contrast{icont2},contrast{icont2});
    inxk = inxk + 1;
    for icont1 = 1 : ncont_plot       

        if icont1~=icont2        
            mM(:,inxk) = mMtr(:,icont2,icont1);
            sM(:,inxk) = sMtr(:,icont2,icont1);
            M(:,inxk,:,icont2)=Mtr(:,icont2,icont1,:);
            grpstr{inxk}= sprintf('S%s-%s',contrast{icont1},contrast{icont2});
            inxk = inxk + 1;
            mM(:,inxk) = mMte(:,icont2,icont1);
            sM(:,inxk) = sMte(:,icont2,icont1);
            M(:,inxk,:,icont2)=Mte(:,icont2,icont1,:);
            grpstr{inxk}= sprintf('%s-S%s',contrast{icont1},contrast{icont2});
            inxk = inxk + 1;

        end
    end

    sel_thrinx = [1 :2:length(CELLSEL_THR)];
    linewidth=2;
    CLR={[0.7 0.7 0.7],'k','b','g'}
    linestyle={':','--','-','--','-'}
    figure('Position',[100 100 300 230]); hold on;
    for icont1= 1 : 1+(ncont_plot-1)*2

        x = CELLSEL_THR(sel_thrinx);
        errorbar(x',mM(sel_thrinx,icont1),sM(sel_thrinx,icont1),'LineWidth',linewidth,'Color',CLR{icont1},'LineStyle',linestyle{icont1});

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
    legend(grpstr);
end

selinx=[3 5 7];
icont=2
a=squeeze(M(selinx,2,:,icont)-M(selinx,3,:,icont))';
[h p]=ttest(a,[],0.05/length(selinx))
pc = p*length(selinx)
%%  ------------
% inx_diag = (0:ncont_plot:ncont_plot*(ncont_plot-1))+(1:ncont_plot);
% inx_offdiag = setdiff(1:ncont_plot^2,inx_diag);
% mM = zeros(size(mMtr,1),length(inx_diag)+2*length(inx_offdiag));
% sM = zeros(size(mMtr,1),length(inx_diag)+2*length(inx_offdiag));
% grpstr =cell(1,length(inx_diag)+2*length(inx_offdiag));
% inx_k=1;
% for inx1 = 1: length(inx_diag)
%     mM(:,inx_k) = mMtr(:,inx1,inx1);
%     sM(:,inx_k) = sMtr(:,inx1,inx1);
%     grpstr{inx_k}= sprintf('%s-%s',contrast{inx1},contrast{inx1});
%     inx_k = inx_k +1;
% end
% for inx1 = inx_offdiag
%     [i,j]=ind2sub([ncont_plot ncont_plot],inx1);
%     mM(:,inx_k) = mMtr(:,j,i);
%     sM(:,inx_k) = sMtr(:,j,i);
%     grpstr{inx_k}= sprintf('%s-%s',contrast{i},contrast{j});
%     mM(:,inx_k+1) = mMte(:,j,i);
%     sM(:,inx_k+1) = sMte(:,j,i);
%     grpstr{inx_k+1}= sprintf('%s-%s',contrast{i},contrast{j});
%     inx_k = inx_k +2;
% end
% 
% sel_thrinx = [1 :2:length(CELLSEL_THR)];
% linewidth=2;
% CLR={'b','g','b','b','g','g'}
% linestyle={'-','-','--',':','--',':'}
% figure('Position',[100 100 300 230]); hold on;
% for icont1= 1 : ncont_plot
% 
%     x = CELLSEL_THR(sel_thrinx);
%     errorbar(x',mM(sel_thrinx,icont1),sM(sel_thrinx,icont1),'LineWidth',linewidth,'Color',CLR{icont1},'LineStyle',linestyle{icont1});
% 
% end
% legend(grpstr);
% if strcmp(SEL_METHOD,'WEIGHTL1')
%     set(gca,'XTick',[0 0.1 0.2 0.3])
%     tickv = num2cell([0 0.1 0.2 0.3]);
%     tickstr = cellfun(@num2str,tickv, 'UniformOutput', false);
%     set(gca,'XTickLabel',tickstr);
%     set(gca,'FontSize',22);
%     xlim([-0.05 0.35])
% else
%     set(gca,'XTick',[0 10 20 30])
%     tickv = num2cell([0 10 20 30]);
%     tickstr = cellfun(@num2str,tickv, 'UniformOutput', false);
%     tickstr{1}='ALL';
%     set(gca,'XTickLabel',tickstr);
%     set(gca,'FontSize',22);
%     xlim([-5 35])
% end
% 
% 
% 
% %----- plot ration of cells selected
% if strcmp(SEL_METHOD,'WEIGHTL1')
%     a=bsxfun(@rdivide,Nc,reshape(Nc0,[1 1 length(Nc0)]));
%     figure('Position',[200 100 300 230]);hold on;
%     for icont = 1:ncont_plot
%     errorbar(x,mean(a(sel_thrinx,icont,:),3),std(a(sel_thrinx,icont,:),0,3)/sqrt(length(size(a,3))),'LineWidth',2,'Color',CLR{icont})       
%     end
%     set(gca,'XTick',[0 0.1 0.2 0.3])
%     tickv = num2cell([0 0.1 0.2 0.3]);
%     tickstr = cellfun(@num2str,tickv, 'UniformOutput', false);
%     set(gca,'XTickLabel',tickstr);
%     set(gca,'FontSize',22);
%     xlim([-0.05 0.35])
%     ylim([0 0.4]);
%     box off
% end
% % 
% % ------ stat
% hs=[];
% ps=[];
% for icont = 1 : ncont_plot
%     Mi = squeeze(M(sel_thrinx(2:end),:,icont,:));
%     FWE = 0.05/size(Mi,1);
% %     [h,p]=ttest(squeeze(Mi(1,1,:)), squeeze(Mi(1,2,:)),FWE);
%     [h,p]=ttest(squeeze(Mi(:,1,:))'- squeeze(Mi(:,2,:))',[],FWE);
%     hs(:,icont)=h;
%     ps(:,icont)=p*size(Mi,1);
%     
%     Mi2 = reshape(Mi, [size(Mi,1)*size(Mi,2) size(Mi,3)])';
%     [p t st]=anova1(Mi2);
%     [c,m,h,nms] = multcompare(st,'alpha',0.05);
% end


% 
% Mi = squeeze(M(2:7,:,icont,:));
% FWE = 0.05/size(Mi,1);
% [h,p]=ttest(squeeze(Mi(1,1,:)-Mi(6,2,:))',[],FWE);
% hs(:,icont)=h;
% ps(:,icont)=p*size(Mi,1);