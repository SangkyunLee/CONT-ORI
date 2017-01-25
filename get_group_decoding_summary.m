
%% analysis of group decoding acc data
clear all
close all

exp_type={'AN','AN_0TO150','AWAKE'};


%--- AN
ORI_list=[-15 0 30 90];
contrasts = [ 100 40 20];
sel_compset = 1:6;
ncellgrp_plot=[3 4 5];
iexp_type=1;
seslist = [(1:8) 11 12 15 16];

%--- AN_0TO150
% ORI_list=[0 30 35 60 90 12 150];
% contrasts = [ 100 40];
% sel_compset = [1 2 3 4 7 12 16 19 21];   
% ncellgrp_plot=[3 4 5];
% iexp_type=2;
% seslist=1:7;


ctm=0.6;
DATA_thr_str = 'thr5';
rThr =0.9;
mThr =0.6;


fnsave = sprintf('SUMMARYTHR_AN1-16_DEC_CELLGRP_ctm%0.2f.mat',ctm);
fullfnsave = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fnsave);




cuts = zeros(length(sel_compset),length(contrasts), length(seslist));
CELL_CBT=cell(length(sel_compset),length(ncellgrp_plot),length(contrasts),length(seslist));
CELL_CBT2 = -1*ones(200,length(sel_compset),length(ncellgrp_plot),length(contrasts),length(seslist));
CELLINX = cell(max(seslist),1);
DEC = cell(max(seslist),1);


bdisp=0;
for ises0 = 1: length(seslist)
    ises = seslist(ises0);
    fndata = sprintf('DEC_CELLGRP_ctm%0.2f_ses%d.mat',ctm,ises);
    fullfndata = fullfile('../NEW_DECODING/',exp_type{iexp_type},DATA_thr_str, fndata);
    fullfndata1 = strrep(fullfndata,'DEC_CELLGRP','CELLGRP');
    
    if ~(exist(fullfndata,'file') && exist(fullfndata1,'file'))
        continue;
    end
    
    load(fullfndata);
    load(fullfndata1);


    for icont = 1 : length(contrasts)
        grpinx={};
        grpACC={};
        if bdisp
            figure;
        end
        
        for icomp2 = 1 : length(sel_compset) 
            icomp =sel_compset(icomp2); 
            kk=1;
            nvals = zeros(61,length(ncellgrp_plot));
            for incell = ncellgrp_plot

                ACC = single(deaccN{incell,icomp,icont});

                a = hist(ACC,40:100); 
                nvals(:,kk) = a./sum(a); 

                len_HACC = length(find(ACC>=deaccALL(icomp,icont)*100));
                if len_HACC>50
                    cellgrp1 = find(ACC>=deaccALL(icomp,icont)*100 & ACC>(mThr*100));
                    if deaccALL(icomp,icont)>mThr
                        cut = deaccALL(icomp,icont)*100;
                    else
                        cut =mThr*100;
                    end
                else
                    cellgrp1 = find(ACC >= rThr*deaccALL(icomp,icont)*100  & ACC>(mThr*100));
                    if rThr * deaccALL(icomp,icont)*100>(mThr*100)
                        cut = deaccALL(icomp,icont)*100*rThr;
                    else
                        cut =mThr*100;
                    end

                end
                if ~isempty(cellgrp1)

                    ACC_sel = double(grplistALL{incell}(cellgrp1,:));    
                    ncase = size(ACC_sel,1);
                    
                    uniquecell = unique(ACC_sel(:));
                    cellnR=zeros(length(uniquecell),2);
                    cellnR2 = zeros(Ncell,1);
                    for icell = 1 : length(uniquecell)            
                        cellnR(icell,:) = [uniquecell(icell) length(find(ACC_sel(:)==uniquecell(icell)))/ncase];
                        cellnR2(uniquecell(icell)) = length(find(ACC_sel(:)==uniquecell(icell)))/ncase;

                    end
                    [a b]=sort(cellnR(:,2),'descend');       


                    CELL_CBT{icomp2,kk,icont,ises0} = cellnR(b,:);
                    CELL_CBT2(1:length(cellnR2),icomp2,kk,icont,ises0)=cellnR2;
                end
                kk = kk+1;
            end
            cuts(icomp2,icont,ises0)=cut;
            if bdisp
                subplot(2,3,icomp2); 
                plot(40:100,nvals);
                hold on; plot(deaccALL(icomp,icont)*100*ones(1,2), [0 0.1],'k');
                plot(cut*ones(1,2), [0 0.1],'k--');
                xlim([50 100]);
                title([condlabel{icomp,icont}]);                                
            end
        end
    end
    CELLINX{ises0}=cellinx_sel;
    DEC{ises0}=deaccALL;
end


Info.ORI_list = ORI_list;
Info.contrasts = contrasts;
Info.sel_compset = sel_compset;
Info.ncellgrp_plot = ncellgrp_plot;
Info.iexp_type = iexp_type;
Info.seslist = seslist;
Info.ctm = ctm;
Info.DATA_thr_str = DATA_thr_str;
Info.rThr = rThr;
Info.mThr = mThr;
save(fullfnsave,'cuts','CELL_CBT','CELL_CBT2','CELLINX','DEC','Info');

%% --------------
% 
% sel_compset = Info.sel_compset;
% contrasts = Info.contrasts;
% seslist = Info.seslist;
% inx_NCG =3; % index of 'n'-cell group
% bins =[0.0 :0.1: 1];
% QuantInfo=zeros(3,3,length(sel_compset),length(contrasts),length(seslist));
% P = zeros(length(bins)-1, 3, length(sel_compset), length(seslist));
% 
% for ises = 1 : length(seslist)      
% 
%     for icomp = 1 :length(sel_compset)        
%         %CELL_CBT2(1:length(cellnR2),icomp2,kk,icont,ises0)
%         M = CELL_CBT2(:,icomp,inx_NCG,:,ises);
%         M = squeeze(M);
%         M = M(M(:,1)>=0,:);
%         QuantInfo(:,:,icomp,ises)=quantile(M,[0.25 0.5 0.75],1);
%         
%         [cnt]=histc(M,bins);
%         P(:,:,icomp,ises)=cnt(1:end-1,:);
%     end
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %-------------- select a session and get summary
% ORI_list =Info.ORI_list;
% ises=16;
% icomp=3;
% ncell=3;
% ncont=2;
% a=zeros(20,ncell*ncont*2);
% condset=combnk(ORI_list,2)'; 
% condsets_disp=[];
% kk=1;
% maxlen=0;
% for icont=1:3    
%     for incell= 2: 3
%         if ~isempty(CELL_CBT{icomp,incell,icont,ises})
%             inxcell =CELL_CBT{icomp,incell,icont,ises}(:,2)>0.1;
%             M = CELL_CBT{icomp,incell,icont,ises};
%             tmp =[CELLINX{ises}(M(inxcell,1))' M(inxcell,2)*100 ];
%             a(1:size(tmp,1),kk:kk+1)=tmp;
%             maxlen = max(maxlen, size(tmp,1));
%         end
%         kk = kk + 2;        
%     end
%     
%     condsets_disp= [ condsets_disp condset(:,icomp)];
% end
% a=a(1:maxlen,:);
% a=uint8(a)
% b=uint8([cuts(icomp,:,ises); DEC{ises}(icomp,:)*100])
% 
% %---- distribution with a given threshold cThr
% 
% nc_inx =3;
% cont_inx=[1 3];
% bins =(0.05:0.1:0.75);
% cThrs=[0.1];
% figure;
% for icthr= 1: length(cThrs)
%     cThr = cThrs(icthr);
%     P = zeros(length(bins),max(sel_compset), length(seslist));
% 
%     for ises = 1: length(seslist)           
%         for icomp =sel_compset
%             M = CELL_CBT2(:,icomp,nc_inx,:,ises);
%             M = squeeze(M);
%             M = M(M(:,1)>=0,:);
%             M1 = M(M(:,cont_inx(1))> cThr | M(:,cont_inx(1))>cThr,:);
%             [cnt]=hist(abs(M1(:,cont_inx(1))-M1(:,cont_inx(2))),bins);
%             P(:,icomp,ises)=cnt/size(M1,1);
%         end
%     end
%     mP = mean(P,3);
%     sP = std(P,0,3)/sqrt(length(seslist));
%     subplot(length(cThrs),1,icthr);
%     errorbar(bins'*ones(1,size(mP,2)),mP,sP);
%     title(['cThr = ' num2str(cThr)]);
%     ylim([0 0.5]);
% end
% 
% 
% %---- distribution with a given threshold cThr bins
% 
% nc_inx =3;
% cont_inx=[1 3];
% bins =[0.0 0.2 0.4 0.8 1];
% %cThrs=[0.0 0.3 0.6 1];
% cThrs=[0.0 0.1 0.2 1];
% QuantInfo=zeros(3,3,max(sel_compset),length(seslist));
% figure;
% for icthr= 1: length(cThrs)-1
%     cThr1 = cThrs(icthr);
%     cThr2 = cThrs(icthr+1);
%     P = zeros(length(bins)-1,max(sel_compset), length(seslist));
% 
%     for ises = 1 : length(seslist)            
%         for icomp =sel_compset
%             M = CELL_CBT2(:,icomp,nc_inx,:,ises);
%             M = squeeze(M);
%             M = M(M(:,1)>=0,:);
%             if icthr ==1,
%                 QuantInfo(:,:,icomp,ises)=quantile(M,[0.25 0.5 0.75],1);
%             end               
%             
%             inx1 = M(:,cont_inx(1))> cThr1 & M(:,cont_inx(1))<= cThr2;
%             inx2 = M(:,cont_inx(2))> cThr1 & M(:,cont_inx(2))<= cThr2;
%             
%             M1 = M(inx1 | inx2,:);
%             [cnt]=histc(abs(M1(:,cont_inx(1))-M1(:,cont_inx(2))),bins);
%             P(:,icomp,ises)=cnt(1:length(bins)-1)/size(M1,1);
%         end
%     end
%     
%     mP = mean(P,3);
%     sP = std(P,0,3)/sqrt(length(seslist));
%     subplot(length(cThrs)-1,1,icthr);
%     bin1 = diff(bins)/2+bins(1:length(bins)-1);    
%     x = bin1'*ones(1,size(mP,2));
%     errorbar(x,mP,sP);
%     title(['cThrbin = [' num2str(cThr1) '-' num2str(cThr2) ']']);
% %     ylim([0 0.5]);
% end

