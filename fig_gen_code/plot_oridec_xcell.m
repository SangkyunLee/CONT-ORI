% clear all
% close all;




% M{1} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS\AN\thr5\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
% M{2} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS\AN\thr5\AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
% M{3} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS\AWAKE_EYE\thr5_eyethr_xy1_p1\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');

M{1} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN1-16SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
M{2} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
M{3} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')


% M{1} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\SHUFFLE_AN1-16SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
% M{2} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\SHUFFLE_AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{3} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\SHUFFLE_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% 

%-- for common training 
% comcont{1,1} = {100,100}; comcont{1,2} = {100, 40};
% comcont{2,1} = {40,40}; comcont{2,2} = {40, 100};


%---- for common testing
comcont{1,1} = {100,100}; comcont{1,2} = {40, 100};
comcont{2,1} = {40,40}; comcont{2,2} = {100, 40};

%---- for common testing with contrast-independence
% comcont{1,1} = {100,100}; comcont{1,2} = {[100 40], 100};
% comcont{2,1} = {40,40}; comcont{2,2} = {[100 40], 40};


nc= size(comcont,1);
D1 = cell(6,6,nc);
D2 = cell(6,6,nc);
D3 = cell(6,6,nc);
IX1 = cell(6,6);
Dthr =0;
dispord =[6 3 5 1 2 4];
ncellord = [ (2:6) 1];
ncell = [NaN 1 3 5 10 20];
for inxcell0 = 1 : 6
    inxcell = ncellord(inxcell0);
    for iori = 1:6
        inx_compori = dispord(iori);

        D =cell(length(M),nc);
        IX =cell(length(M),nc);
        for ises = 1 : length(M)
            %----------------
            Tcomcont = M{ises}.comcont;

            CSinx= NaN*ones(1,nc); % Contrast-Specific index
            CCinx= NaN*ones(1,nc);%Cross-contrast index

            for i = 1 : nc   
                for j = 1 : size(Tcomcont,2)
                    if all(comcont{i,1}{1}== Tcomcont{j}{1}) &&...
                            all(comcont{i,1}{2}== Tcomcont{j}{2})
                        CSinx(i) = j;
                    end
                end
            end

            for i = 1 : nc   
                for j = 1 : size(Tcomcont,2)
                    if all(comcont{i,2}{1}== Tcomcont{j}{1}) &&...
                            all(comcont{i,2}{2}== Tcomcont{j}{2})
                        CCinx(i) = j;
                    end
                end
            end
            
            %------ calculate overlap ratio of cells selected for ORI
            %decoding between two different contrasts        
            
            
            nsub = size(M{ises}.CELLSEL_INX,3);
            IX0 = NaN*ones(nsub,1);
            
            for isub = 1: nsub                
                uc =cell(1,nc);
                CELLinx = cell(1,nc);
                N = cell(1,nc);
                for i = 1 : nc
                    CELLinx0 = M{ises}.CELLSEL_INX{CSinx(i),inx_compori,isub};
                    if isempty(CELLinx0)
                        continue;
                    end
                    
                    CELLinx{i} = cell2mat(CELLinx0(inxcell,:));
                    uc{i} = unique(CELLinx{i});
                    N{i} = length(CELLinx{i});
                end
                if any(cellfun(@isempty,N))
                    continue;
                end
                uc = unique(cell2mat(uc));                
                [a, x]=cellfun(@hist,CELLinx,repmat({uc},[1 nc]), 'UniformOutput', false);
                for i=2:nc
                    if any(x{i}-x{1})
                        error('not matched in xbin');
                    end
                end
                P = cellfun(@rdivide, a, repmat({100},size(a)), 'UniformOutput', false);
                NC =  length(M{ises}.CELLSEL_INX{1,1,isub}{1});
                

                
               
                OP=P{1};
                for i=2:nc
                    OP=OP.*P{i};
                end    
                % when not corrected for dispersion: response variability
                % over trials
%                 IX0(isub) = mean(OP);
                
                % when corrected for dispersion: response variability
                % over trials
                if isnan(ncell(inxcell))
                    IX0(isub) = mean(OP);
                else
                    IX0(isub) =sum(OP)/ncell(inxcell);
                end
                
                
                
                %IX0(isub) = sum(OP)/NC;
            end
            
            IX{ises} = IX0(~isnan(IX0));
            

            % collect ORI decoding accs
            for i = 1 : nc    
                x = squeeze(M{ises}.DEC_SELCELL(inxcell,CSinx(i),inx_compori,:));
                y = squeeze(M{ises}.DEC_SELCELL(inxcell,CCinx(i),inx_compori,:));  
                D{ises,i} = [x(x>0) y(y>0)]; 
            end 
            
        end % ises
        IX1{inxcell0,iori} = cell2mat(IX);

        
        for k = 1: nc
            Dtmp = cell2mat(D(:,k));
            Dtmp = Dtmp(Dtmp(:,k)>Dthr,:);
            D1{inxcell0,iori,k} = Dtmp;
            Dtmp2 = Dtmp(:,2)- Dtmp(:,1);
            D2{inxcell0,iori,k} = Dtmp2;
            D3{inxcell0,iori,k} = Dtmp2./Dtmp(:,1);
            
            
        end
        
        
    end % iori
end %inxcell0

D4 = cell(size(D3));
for iori =1:6
    for icell = 1:6
        for k = 1:2
        D4{icell,iori,k} = D2{icell,iori,k}./D1{end,iori,k}(:,1);
        end
    end
end

%-------- plot n-cell decoding acc change
inc = 2
D3x = cell2mat(D4(:,:,inc));
nrep = size(D3{1,1,inc},1);
D3x = reshape(D3x,[nrep 6 6]); %nrep x ncell x nori

x = squeeze(mean(D3x,3));


len = size(x,1);                   
                   
SEM = std(x)/sqrt(len);               % Standard Error
ts = tinv([0.025  0.975],len-1);      % T-Score
CI = ones(2,1)*mean(x) + ts'*SEM;                      % Confidence Intervals

fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','w','Facealpha',.5);
figure; 
hold on;
fill_rec([1 3 5 10 20 30],CI(1,:)*100,CI(2,:)*100,[0.3 0.3 0.3])
plot([1 3 5 10 20 30],mean(x)*100,'.-','LineWidth',2,'MarkerSize',20,'Color','k');
set(gca,'XTick',[1 3 5 10 20 30])
set(gca,'XTickLabel',{'1','3','5','10','20','all'})
set(gca,'FontSize',22,'linewidth',2)
xlim([0 31])
ylim([-5 0])
% ylim([-10 5])

% [p,table,stats] = anova1(x);
% [c,m,h,nms] = multcompare(stats)
% 
for iori=1:6
    if iori == 1
        Dec=cell2mat(D1(:,iori,inc));
    else
        Dec=Dec+cell2mat(D1(:,iori,inc));
    end
end
Dec=Dec/6; 
nrep = size(Dec,1)/6;
% [p,table,stats] = anova2(Dec,nrep);
% [c,m,h,nms] = multcompare(stats,'estimate','column');
[p,table,stats] = friedman(Dec,nrep);


Dec1 = reshape(Dec,[nrep 6 2]);

clr = {'b','r'};
clra = {[0 0 0.3],[0.3 0 0]};
fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','none','Facealpha',.15);
figure; 
hp=zeros(1,2);
for i = 1 : 2
    x = Dec1(:,:,i);
    len = size(x,1);                   

    SEM = std(x)/sqrt(len);               % Standard Error
    ts = tinv([0.025  0.975],len-1);      % T-Score
    CI = ones(2,1)*mean(x) + ts'*SEM;                      % Confidence Intervals

    hold on;
    fill_rec([1 3 5 10 20 30],CI(1,:),CI(2,:),clra{i})
    hp(i)=plot([1 3 5 10 20 30],mean(x),'.-','LineWidth',2,'MarkerSize',20,'Color',clr{i});
end
set(gca,'XTick',[1 3 5 10 20 30])
set(gca,'YTick',[0.6 0.8 1])
set(gca,'XTickLabel',{'1','3','5','10','20','all'})
set(gca,'FontSize',22,'linewidth',2)
xlim([0 31])
ylim([0.5 1])

legend(hp,'Within-Cont.','Cross-Cont.')
legend(hp,'Within-Cont.','Indep.-Cont.')
legend('boxoff')

k = 1;
ps =[];
for ii=[1 3 6]
figure; 
plot(Dec1(:,ii,1),Dec1(:,ii,2),'.','LineWidth',2,'MarkerSize',20)
hold on; 
plot([0.45 1],[0.45 1],'k','LineWidth',2)
set(gca,'FontSize',22,'linewidth',2)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])
box off

xlim([0.6 1])
ylim([0.6 1])

ps(k) = signrank(Dec1(:,ii,1),Dec1(:,ii,2));
k = k +1;
end

%% plot commonly used cell ratio between contrasts

nrep = size(IX1{1,1},1);
IXs = reshape(cell2mat(IX1),[nrep 6 6]); %nrep x ncell x nori



x = squeeze(mean(IXs,3));


len = size(x,1);                   
                   
SEM = std(x)/sqrt(len);               % Standard Error
ts = tinv([0.025  0.975],len-1);      % T-Score
CI = ones(2,1)*mean(x) + ts'*SEM;                      % Confidence Intervals

fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','w','Facealpha',.5);
figure; 
hold on;
fill_rec([1 3 5 10 20 30],CI(1,:),CI(2,:),[0.3 0.3 0.3])
plot([1 3 5 10 20 30],mean(x),'.-','LineWidth',2,'MarkerSize',20,'Color','k');
set(gca,'XTick',[1 3 5 10 20 30])
set(gca,'XTickLabel',{'1','3','5','10','20','all'})
set(gca,'FontSize',22,'linewidth',2)
xlim([0 20])

%% plot commonly used cell ratio between contrasts after hard thresholding
%------- without using probability, select n-cell and calculate the ratio

% out of 100 cv tests,when cells are selected at least 70 trials,
% the cells were used for overlapping ratio.
% in addition, the denominator was the minmum cell number between 100% and
% 40% contrasts.
SELTHR=0.7; 
nc=2
IX1 = cell(5,6);

dispord =[6 3 5 1 2 4];
ncellord = 2:6;
ncell = [NaN 1 3 5 10 20];
for inxcell0 = 1 : 5
    inxcell = ncellord(inxcell0);
    for iori = 1 : 6
        inx_compori = dispord(iori);
        IX =cell(length(M),nc);
        for ises = 1 : length(M)


            %------ calculate overlap ratio of cells selected for ORI
            %decoding between two different contrasts    
            nsub = size(M{ises}.CELLSEL_INX,3);
            IX0 = NaN*ones(nsub,1);
            
            for isub = 1: nsub  
                isub
                uc =cell(1,nc);
                CELLinx = cell(1,nc);
                N = cell(1,nc);
                for i = 1 : nc
                    CELLinx0 = M{ises}.CELLSEL_INX{i,inx_compori,isub};
                    if isempty(CELLinx0)
                        continue;
                    end
                    
                    CELLinx{i} = cell2mat(CELLinx0(inxcell,:));
                    uc{i} = unique(CELLinx{i});
                    N{i} = length(CELLinx{i});
                end
                if any(cellfun(@isempty,N))
                    continue;
                end
                uc = unique(cell2mat(uc));                
                [a, x]=cellfun(@hist,CELLinx,repmat({uc},[1 nc]), 'UniformOutput', false);
                for i=2:nc
                    if any(x{i}-x{1})
                        error('not matched in xbin');
                    end
                end
                P = cellfun(@rdivide, a, repmat({100},size(a)), 'UniformOutput', false);
      
               
                minN=Inf;
                for i = 1 : nc
                    Px = P{i};
                    [~,inx] = sort(Px);
                    Px(Px<SELTHR)=0;
                    Px = sign(Px);
                    minN = min(minN, length(find(Px)));
                    
                    if i == 1,
                        OP=Px;
                    end
                    OP=OP.*Px;
                end    

                
                % when corrected for dispersion: response variability
                % over trials
                
                %IX0(isub) =sum(OP)/ncell(inxcell);
                if minN==0
                    minN=1;
                end
                IX0(isub) =sum(OP)/minN;
                
            end
            
            IX{ises} = IX0(~isnan(IX0));
            

            
        end % ises
        IX1{inxcell0,iori} = cell2mat(IX);

        

        
    end % iori
end %inxcell0



nrep = size(IX1{1,1},1);
IXs = reshape(cell2mat(IX1),[nrep 5 6]); %nrep x ncell x nori



x = squeeze(mean(IXs,3));


len = size(x,1);                   
                   
SEM = std(x)/sqrt(len);               % Standard Error
ts = tinv([0.025  0.975],len-1);      % T-Score
CI = ones(2,1)*mean(x) + ts'*SEM;                      % Confidence Intervals

fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','w','Facealpha',.5);
figure; 
hold on;
fill_rec([1 3 5 10 20 ],CI(1,:),CI(2,:),[0.3 0.3 0.3])
plot([1 3 5 10 20 ],mean(x),'.-','LineWidth',2,'MarkerSize',20,'Color','k');
set(gca,'XTick',[1 3 5 10 20 ])
set(gca,'XTickLabel',{'1','3','5','10','20'})
set(gca,'FontSize',22,'linewidth',2)
xlim([0 20])


%% ------------------ plot examples
icont1 = 1 ;% 100
icont2 =2; % 40
inc = 3;
iori= 3
isub =1;

C1 = M{1}.CELLSEL_INX{icont1,iori,isub};
C1 = cell2mat(C1(inc,:));
C2 = M{1}.CELLSEL_INX{icont2,iori,isub};
C2 = cell2mat(C2(inc,:));
uc = unique([C1 C2]);
luc = length(uc);
h1= hist(C1,uc);
h2= hist(C2,uc);

hf=figure;%('Position',[480   478   660   400]);
bar((1:luc)'*ones(1,2),[h2' h1'],'BarWidth',1);
xlim([0.5 luc+.5])
xstr=cell(1,luc);
% for i=1:luc
%     xstr{i}=sprintf('Cell%d',i);
% end
% set(gca,'xticklabel',xstr);
set(gca,'FontSize',22,'linewidth',2)
box off


O{1}=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\AN1-16_ORIsc_ctm0.60_fit_ab4.mat');
O1 = O{1}.ORItun(isub);

ord = O1.evtord(:,1,2);
a = O1.mresp(:,1:2,uc);
ma = max(a(:));
sc = zeros(luc,2);
for i = 1: luc
    sc(i,1) = max(a(:,1,i))/ma;
    sc(i,2) = max(a(:,2,i))/ma;
end


%figure('Position',[480   678   600   450]);
figure('Position',[480   678   1300   210]);
clrs={'r','b'};
cellist =[ 2 3 9 10];
mys = zeros(length(cellist),2);
snrs = zeros(length(cellist),2);
for i0 = 1:length(cellist)
    i = cellist(i0);
    subplot(1,4,i0);
    hold on;
    for j = 1 :2
        my = O1.mresp(:,j,uc(i));
        stdy = O1.stdresp(:,j,uc(i));
        mys(i0,j) = abs(my(2)-my(3));
        snrs(i0,j)= mys(i0,j)/sqrt(mean(stdy(2:3).^2));
        Mmy = max(my);
        my = my/Mmy*sc(i,j);
        sy = O1.semresp(:,j,uc(i));
        sy = sy/Mmy*sc(i,j);
        errorbar(ord,my,sy,'linewidth',2,'color',clrs{j}); 
        
    end 
    set(gca,'ytick',[0 0.5 1]);
    set(gca,'ytick',[0 0.5 1]);
    set(gca,'FontSize',22,'linewidth',2)
    xlim([-20 100])
    ylim([0 1])
    %tltstr=sprintf('Cell%d,%.2f',i, abs(my(2)-my(3)));
    tltstr=sprintf('Cell%d',i);
    title(tltstr);
    if i0==1,
        plot([0 0],[0 1],'k--');
        plot([30 30],[0 1],'k--');
    end
%     axis off
end



ratio = h1/100*h2'/100/3




%% %-------- plot n-cell decoding acc change
% inc = 1
% nrep = size(D3{1,1},1);
% ncell=6; nori=6;
% Dx = zeros(nrep,ncell,nori);
% for iori = 1:6
%     for icell = 1:6        
%         Dx(:,icell,iori) = D4{icell,iori,inc};    
%     end
% end
% 
% fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','w','Facealpha',.5);
% figure('Position', [680 71 500 1000]); 
% for iori = 5
%     x = Dx(:,:,iori);
%     len = size(x,1);                   
%                    
%     SEM = std(x)/sqrt(len);               % Standard Error
%     ts = tinv([0.025  0.975],len-1);      % T-Score
%     CI = ones(2,1)*mean(x) + ts'*SEM;                      % Confidence Intervals
% 
% 
%     subplot(6,1,iori);
%     hold on;
%     fill_rec([1 3 5 10 20 30],CI(1,:)*100,CI(2,:)*100,[0.3 0.3 0.3])
%     plot([1 3 5 10 20 30],mean(x)*100,'.-','LineWidth',1.2,'MarkerSize',10,'Color','k');
%     set(gca,'XTick',[1 3 5 10 20 30])
%     set(gca,'XTickLabel',{'1','3','5','10','20','all'})
%     set(gca,'FontSize',16)
%     xlim([0 31])
%     ylim([-15 0])
% %     
%     
%     [p,table,stats] = anova1(x);
% [c,m,h,nms] = multcompare(stats)
% % [nms num2cell(c)]
% end

%% plot signal-to-noise ratio in contrasts

clear all
D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\BASIC_SUMMARY_AN1-16_ctm0.60.mat')
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\BASIC_SUMMARY_AN17-22_ctm0.60.mat')
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\BASIC_SUMMARY_AW23-40_ctm0.60.mat')


j = 1;
S=cell(28,1);
S1=cell(28,1);
S2=cell(28,1);
selcont = [1 2];
for iexp = 1 : 3
     ns = size(D(iexp).S,2);
    for isub = 1 :ns 
        X1 = squeeze(mean(D(iexp).S(isub).V./D(iexp).S(isub).m*100,2)); % fano factor
        X2 = squeeze(mean(sqrt(D(iexp).S(isub).V),2))*100;% noise only
        X = squeeze(mean(D(iexp).S(isub).m./sqrt(D(iexp).S(isub).V),2)); % snr
        if ~isempty(X)
            S{j} = X(:,selcont);
            S1{j} = X1(:,selcont);
            S2{j} = X2(:,selcont);
            j = j + 1;
        end
    end
end
        

S =cellfun(@mean,S, 'UniformOutput', false);
S = cell2mat(S);

S1 =cellfun(@mean,S1, 'UniformOutput', false);
S1 = cell2mat(S1);

S2 =cellfun(@mean,S2, 'UniformOutput', false);
S2 = cell2mat(S2);


signrank(S(:,1),S(:,2))


hf=figure;%('Position',[680   678   640   420]);

model_series =nanmean(S,1);
model_error = nanstd(S,0,1)/sqrt(sum(~isnan(S(:,1))));
h = bar(model_series);
set(h,'BarWidth',0.4,'FaceColor',[0.91 0.91 0.91]);    % The bars will now touch each other
set(gca,'XTicklabel','100|40')
hold on;
errorbar( model_series, model_error, 'k', 'linestyle', 'none','linewidth',2);
set(gca,'FontSize',22,'linewidth',2);
xlim([0.5 2.5])
box off      

rF1 = bsxfun(@rdivide, bsxfun(@minus,S1,S1(:,1)),S1(:,1));
rF2 = bsxfun(@rdivide, bsxfun(@minus,S2,S2(:,1)),S2(:,1));
rF = [rF1(:,2) rF2(:,2)]*100;

figure;
errorbar(mean(S1,1),std(S1,0,1)/sqrt(size(S1,1)));

figure; hold on;
plot(rF','k.','linewidth',2);
plot(median(rF),'x','MarkerSize',20)
xlim([0.5 2.5])
set(gca,'FontSize',22);
set(gca,'XTick',[1 2])
set(gca,'XTicklabel','|')

plot([0 3],[0 0],'k')

[h p]=ttest(rF(:,2))

signrank(S1(:,1),S1(:,2))
%%
