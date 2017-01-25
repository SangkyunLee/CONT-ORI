
% comcont={{100,100,100},{100,100,40},{40,100,100},{40,100,40}, ...
%     {40,40,100},{40,40,40},{100,40,100},{100,40,40}};

%DEC_SELCELL = zeros(length(CELLSEL_THRs),length(comcont),length(ORI_compindexset),nses);
clear all

% M{1}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM3_AN1-16_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{2}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM3_AN17-22_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{3}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\GAM3_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')

% M{1}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM4_AN1-16_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{2}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM4_AN17-22_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{3}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\GAM4_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')

M{1}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM5_AN1-16_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
M{2}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM5_AN17-22_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
M{3}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\GAM5_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')


% M{1}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM6_AN1-16_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{2}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\GAM6_AN17-22_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{3}=load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\GAM6_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')



nc= size(M{1}.comcont,2);
D1 = cell(6,6,nc);
IX1 = cell(6,6);
Dthr =0;
dispord =[6 3 5 1 2 4];
ncellord = [ (2:6) 1];
ncell = [NaN 1 3 5 10 20];

for inxcell0 = 6 %1 : 6
    inxcell = ncellord(inxcell0);
    for iori = 1 : 6
        inx_compori = dispord(iori);

        D =cell(length(M),nc);
        IX =cell(length(M),nc);
        for ises = 1 : length(M)
            %----------------
            
            
%             nsub = size(M{ises}.CELLSEL_INX,3);
%             IX0 = NaN*ones(nsub,1);
%             
%             for isub = 1 : nsub                
%                 uc =cell(1,nc);
%                 CELLinx = cell(1,nc);
%                 N = cell(1,nc);
%                 for i = 1 : nc
%                     CELLinx0 = M{ises}.CELLSEL_INX{i,inx_compori,isub};
%                     if isempty(CELLinx0)
%                         continue;
%                     end
%                     
%                     CELLinx{i} = cell2mat(CELLinx0(inxcell,:));
%                     uc{i} = unique(CELLinx{i});
%                     N{i} = length(CELLinx{i});
%                 end
%                 if any(cellfun(@isempty,N))
%                     continue;
%                 end
%                 uc = unique(cell2mat(uc));                
%                 [a, x]=cellfun(@hist,CELLinx,repmat({uc},[1 nc]), 'UniformOutput', false);
%                 for i= 2 : nc
%                     if any(x{i}-x{1})
%                         error('not matched in xbin');
%                     end
%                 end
%                 P = cellfun(@rdivide, a, repmat({100},size(a)), 'UniformOutput', false);
%                 NC =  length(M{ises}.CELLSEL_INX{1,1,isub}{1});                
% 
%                
%                 OP=P{1};
%                 for i = 2 : nc
%                     OP = OP.*P{i};
%                 end    
%                 % when not corrected for dispersion: response variability
%                 % over trials
% %                 IX0(isub) = mean(OP);
%                 
%                 % when corrected for dispersion: response variability
%                 % over trials
%                 if isnan(ncell(inxcell))
%                     IX0(isub) = mean(OP);
%                 else
%                     IX0(isub) = sum(OP)/ncell(inxcell);
%                 end
%                 %IX0(isub) = sum(OP)/NC;
%             end
%             
%             IX{ises} = IX0(~isnan(IX0));
            

            % collect ORI decoding accs
            for i = 1 : nc    
                x = squeeze(M{ises}.DEC_SELCELL(inxcell,i,inx_compori,:));               
                D{ises,i} = x(x>0); 
            end 
            
        end % ises
        IX1{inxcell0,iori} = cell2mat(IX);
        for k = 1 : nc
            D1{inxcell0,iori,k} = cell2mat(D(:,k));        
        end
        
        
    end % iori
end %inxcell0

inxcell = 6;
Di = squeeze(D1(inxcell,:,:));

for iori = 1 : 6
    if iori == 1
        Dtmp = cell2mat(Di(iori,:));
    else
        Dtmp = Dtmp + cell2mat(Di(iori,:));
    end
end
Dtmp = Dtmp/6;

% comcont={{100,100,100},{100,100,40},{40,100,100},{40,100,40}, ...
%     {40,40,100},{40,40,40},{100,40,100},{100,40,40}};



figure; hold on;
% {40,40,40} vs {40,100,40}
plot(Dtmp(:,6),Dtmp(:,4),'.','markersize',20);
plot([0.5 1],[0.5 1],'k','linewidth',2)
signrank(Dtmp(:,4),Dtmp(:,6))
set(gca,'FontSize',20)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])
axis equal
xlim([0.5 1]); ylim([0.5 1])




figure; hold on;
% {40,40,40}vs{100,40,40}
plot(Dtmp(:,6),Dtmp(:,8),'.','markersize',20);
plot([0.5 1],[0.5 1],'k','linewidth',2)
signrank(Dtmp(:,5),Dtmp(:,8))
set(gca,'FontSize',20)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])
axis equal
xlim([0.5 1]); ylim([0.5 1])


% comcont={{100,100,100},{100,100,40},{40,100,100},{40,100,40}, ...
%     {40,40,100},{40,40,40},{100,40,100},{100,40,40}};
figure; hold on;

plot(Dtmp(:,4),Dtmp(:,6),'.','markersize',20);
plot([0.5 1],[0.5 1],'k','linewidth',2)
signrank(Dtmp(:,5),Dtmp(:,8))
set(gca,'FontSize',20)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])
axis equal
xlim([0.5 1]); ylim([0.5 1])


%%%%%%%%%%%%%%%%%
figure; hold on;
% {100,100,100}vs{100,100,40}
plot(Dtmp(:,1),Dtmp(:,2),'.','markersize',20);
plot([0.5 1],[0.5 1],'k','linewidth',2)

set(gca,'FontSize',20)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])
axis equal
xlim([0.5 1]); ylim([0.5 1])

figure; hold on;
% {40,40,100}vs{40,40,40}
plot(Dtmp(:,5),Dtmp(:,6),'.','markersize',20);
plot([0.5 1],[0.5 1],'k','linewidth',2)

set(gca,'FontSize',20)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])
axis equal
xlim([0.5 1]); ylim([0.5 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% comcont={{100,100,100},{100,100,40},{40,100,100},{40,100,40}, ...
%     {40,40,100},{40,40,40},{100,40,100},{100,40,40}};
ps=zeros(1,4);
for ii =1:4

figure(ii); hold on;
x= Dtmp(:,(ii-1)*2+1);
y= Dtmp(:,(ii)*2);
plot(x,y,'.','markersize',20);
plot([0.5 1],[0.5 1],'k','linewidth',2)

set(gca,'FontSize',20)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])
axis equal
xlim([0.5 1]); ylim([0.5 1])
ps(ii) = signrank(x,y)

end