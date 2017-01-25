%normalized plot to within-session decoding performance
% test Contrast versus Time effect

clear all
close all;



%datapath='/home/slee/data/data_2photon/matlab_2ndLev/NEW_DECODING_NOBIAS_ZMEAN'
datapath = 'Z:/data_2photon/matlab_2ndLev/NEW_DECODING_NOBIAS_ZMEAN';

M{1,1} = load(fullfile(datapath,'/AN/thr5/P1-P2_AN17-22CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'));
M{1,2} = load(fullfile(datapath,'/AN/thr5/P2-P1_AN17-22CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'));
M0{1,1} = load(fullfile(datapath,'/AN/thr5/P1-AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat'));
M0{1,2} = load(fullfile(datapath,'/AN/thr5/P2-AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat'));

M{2,1} = load(fullfile(datapath,'AWAKE_EYE/thr5_eyethr_xy1_p1/P1-P2_AW23-40CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'));
M{2,2} = load(fullfile(datapath,'AWAKE_EYE/thr5_eyethr_xy1_p1/P2-P1_AW23-40CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'));
M0{2,1} = load(fullfile(datapath,'AWAKE_EYE/thr5_eyethr_xy1_p1/P1-SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat'));
M0{2,2} = load(fullfile(datapath,'AWAKE_EYE/thr5_eyethr_xy1_p1/P2-SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat'));





nc= 2;
D1 = cell(6,6);
B1 = cell(6,6);

dispord =[6 3 5 1 2 4];% order with direction difference
ncellord = [ (2:6) 1];
ncell = [NaN 1 3 5 10 20]; % NaN for all cells
explist=[2 5]; % experiment list  5 for AWAKE_EYE
comcontorder=[1 3 4 2]; % original order: 100-100, 40-40, 100-40, 40-100 

for inxcell0 = 1 : 6
    inxcell = ncellord(inxcell0);
    for iori = 1 : 6
        inx_compori = dispord(iori);

        for ises = 1 : size(M,1)
            % collect ORI decoding accs
            nsub = size(M{ises,1}.DEC_SELCELL,2);
            D0 = cell(nsub,1);
            B0 = cell(nsub,1);
            nseg = size(M,2);
            for isub = 1: nsub
                x=zeros(nc*nc,nseg);
                % decoding acc within segment data
                y= zeros(nc*nc,nseg);
                for k = 1 : nseg
                    %%DEC_SELCELL{icomp0,ises}(ithr,icont2,icont)
                    K = M{ises,k}.DEC_SELCELL{inx_compori,isub};                    
                    if isempty(K)
                        continue;
                    end

                    x2 = [];                    
                    for ic1 = 1 : nc                                    
                        x0 = squeeze(K(inxcell,:,ic1));                        
                        if ic1 ==1,
                          x2 = zeros(size(x0,1),nc*nc);
                        end
                        x2(:,(ic1-1)*nc+1 : ic1*nc)=x0;                                    
                    end
                    %%DEC_SELCELL(ithr,icont,icomp0,ises)                         
                    y(:,k) = M0{ises,k}.DEC_SELCELL(inxcell,comcontorder,inx_compori,isub);                    
                    x(:,k)=x2;                    
                end
                

                if sum(abs(x(:)))>0                   
                    D0{isub}=x;
                    B0{isub}=y;                    
                end
            end
            D1{inxcell0,iori} = [D1{inxcell0,iori} cell2mat(D0')];
            B1{inxcell0,iori} = [B1{inxcell0,iori} cell2mat(B0')];
        end % ises
    end % iori
end %inxcell0
% cross-session, Time effect
Y= cell2mat(D1(end,:)');
Y = reshape(Y,[4 6 28]);
Y=permute(Y,[3 2 1]);

% within session, no time effect
Y0= cell2mat(B1(end,:)');
Y0 = reshape(Y0,[4 6 28]);
Y0 = permute(Y0,[3 2 1]);
% %-------------------------






inx1=1;% 100-100
inx2=3; % 40-100 

inx1=4; % 40-40
inx2=2; % 100 -40


A1=Y(1:2:end,:,inx1);
B1=Y(1:2:end,:,inx2);
A2=Y(2:2:end,:,inx1);
B2=Y(2:2:end,:,inx2);

A01=Y0(1:2:end,:,inx1);
B01=Y0(1:2:end,:,inx2);
A02=Y0(2:2:end,:,inx1);
B02=Y0(2:2:end,:,inx2);

% relative performance
C= ((B02)./A02 + (B01)./A01)/2 ; % contrast effect only
T= ((A1)./A02 + (A2)./A01)/2; % time effect only
CT = ((B1)./A02 + (B2)./A01)/2; %contrast-time effect together

signrank(mean(C,2),mean(CT,2))
% time vs contrast-change
figure; hold on; plot(mean(T,2),mean(C,2),'k.','MarkerSize',15);
plot(mean(mean(T,2)),mean(mean(C,2)),'k+','MarkerSize',20);
plot([0.8 1.02],[0.8 1.02],'k','LineWidth',2)
set(gca,'FontSize',20)
set(gca,'XTick',[0.8 0.9 1]);
set(gca,'YTick',[0.8 0.9 1]);
ylim([0.8 1.02])
xlim([0.8 1.02])

% time vs contrast-change,time
figure; hold on; plot(mean(T,2),mean(CT,2),'k.','MarkerSize',15);
plot(mean(mean(T,2)),mean(mean(CT,2)),'k+','MarkerSize',20);
plot([0.8 1.02],[0.8 1.02],'k','LineWidth',2)
set(gca,'FontSize',20)
set(gca,'XTick',[0.8 0.9 1]);
set(gca,'YTick',[0.8 0.9 1]);
ylim([0.8 1.02])
xlim([0.8 1.02])


%%

% Q =                Within-Contrast | Cross-Contrast
%                  ------------------|-------------------
%    Within-ses   |
%    ---------------------------------------------------
%    Cross-ses    |


Q =[mean([A01 A02],2) mean([B01 B02],2);
    mean([A1 A2],2) mean([B1 B2],2)];

 [p,table,stats] = anova2(Q,14);
  [c,m,h,nms] = multcompare(stats,'estimate','row');

% 

 
 
