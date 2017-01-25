clear all
close all;

inxcell = 1;


M{1} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN1-16SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
M{2} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
M{3} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')

M{1} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN1-16SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
M{2} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
M{3} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')

% M{1} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS\AN\thr5\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
% M{2} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS\AN\thr5\AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M{3} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS\AWAKE_EYE\thr5_eyethr_xy0p5_p100\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')



%-- for common training 
comcont{1,1} = {100,100}; comcont{1,2} = {100, 40};
comcont{2,1} = {40,40}; comcont{2,2} = {40, 100};


%---- for common testing
comcont{1,1} = {100,100}; comcont{1,2} = {40, 100};
comcont{2,1} = {40,40}; comcont{2,2} = {100, 40};

%---- for common testing with contrast-independence
% comcont{1,1} = {100,100}; comcont{1,2} = {[100 40], 100};
% comcont{2,1} = {40,40}; comcont{2,2} = {[100 40], 40};


nc= length(comcont);
D1 = cell(nc,6);
D2 = cell(nc,6);
D3 = cell(nc,6);
Dthr =0
for inx_compori = 1 : 6

    D =cell(length(M),nc);
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



        for i = 1 : nc    
            x = squeeze(M{ises}.DEC_SELCELL(inxcell,CSinx(i),inx_compori,:));
            y = squeeze(M{ises}.DEC_SELCELL(inxcell,CCinx(i),inx_compori,:));  
    %         x=x([(1:20) 22]);
    %         y=y([(1:20) 22]);
            D{ises,i} = [x(x>0) y(y>0)]; 
        end       
    end


    % D1 = [ cell2mat(D(:,2))];         
    % D1 = [ cell2mat(D(:,1)); cell2mat(D(:,2))];  

    for k = 1: nc
    Dtmp = cell2mat(D(:,k));
    Dtmp = Dtmp(Dtmp(:,k)>Dthr,:);
    D1{k,inx_compori} = Dtmp;
    D2{k,inx_compori} = Dtmp(:,1)- Dtmp(:,2);
    D3{k,inx_compori} = (Dtmp(:,1)- Dtmp(:,2))./Dtmp(:,1);
    end
end

%--------- 2way-anova

icont = 2;
dispord =[6 3 5 1 2 4];
figure;
X=[];Y=[];
for i=dispord
    x=D1{icont,i}(:,1);
    y=D1{icont,i}(:,2);
    X=[X x];    
    Y=[Y y];
end
plot(X,Y,'.','MarkerSize',20);
lgdstr={'\Delta60','\Delta90','\Delta30','\Delta100(105)','\Delta40(45)','\Delta10(15)'}
[legend_h,object_h,plot_h,text_strings] =legend(lgdstr(dispord),'color','none');
% legend boxoff
hold on;
%plot([min(X(:)) max(X(:))],[min(X(:)) max(X(:))],'k')
plot([0.45 1],[0.45 1],'k','LineWidth',2)
set(gca,'FontSize',20)
set(gca,'XTick',[0.6 0.8 1])
set(gca,'YTick',[0.6 0.8 1])

box off
axis equal
xlim([0.45 1])
ylim([0.45 1])

D1a =cell2mat(D1(icont,:));
D1a=[ D1a(:,1:2:end); D1a(:,2:2:end)]; 

[p,table,stats] = anova2(D1a,size(D1a,1)/2);
% [c,m,h,nms] = multcompare(stats, 'estimate','column')
% 
% 
% [p,table,stats] = anova1((X-Y));
% [c,m,h,nms] = multcompare(stats)
% 
% p=ones(1,6);
% for i=1:6
%     p(i) = signrank(X(:,i),Y(:,i));
% end
% pc = p*6

% X=cell(6,1);
% Y=cell(6,1);
% for i=1:6
%     X{i} = D3{1,i}(:);
%     Y{i}=i*ones(size(X{i},1),1);
% end
% 
% [p,t,st] = anova1(cell2mat(X),cell2mat(Y));
% [c,m,h,nms] = multcompare(st);



%%
% D1 = D1(D1(:,1)>0.7,:);
% 
% (D1(:,1)-D1(:,2))*100
% p = ranksum(D1(:,1),D1(:,2))*12
% figure; hold on;
% plot(D1(:,1),D1(:,2),'.')
% plot([min(D1(:)) max(D1(:))],[min(D1(:)) max(D1(:))]);
% 




%%
Dm =cell(size(D1));
for ii = 1: size(D1,1)
    for jj = 1 : size(D1,2)
        Dm{ii,jj} = D1{ii,jj}(:,1);
    end
end
Dm = mean(cell2mat(Dm),2);
Dm = reshape(Dm,[size(Dm,1)/2 2]);

figure; 
model_series = nanmean(Dm,1);
model_error = nanstd(Dm,0,1)./sqrt(size(Dm,1));
h = bar(model_series);
set(h,'BarWidth',0.4,'FaceColor',[0.91 0.91 0.91]);    % The bars will now touch each other
set(gca,'XTicklabel','100%-100%|40%-40%')
hold on;
errorbar( model_series, model_error, 'k', 'linestyle', 'none','linewidth',2);

set(gca,'FontSize',20);
xlim([0.5 2.5])
ylim([0.5 1])
box off

signrank(Dm(:,1),Dm(:,2));

%%
% clear all
% close all;
% inx_compori = 3;
% inxcell = 1;
% 
% 
% inx_sub=[(1:8) 11 12 15 16];
% M1 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
% M2 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AN\thr5\AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% M3 = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING\AWAKE_EYE\thr5_eyethr_xy1_p1\SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
% 
% 
% 
% 
% D=[];
% inx_comcont= [1 6]
% x = squeeze(M1.DEC_SELCELL(inxcell,inx_comcont(1),inx_compori,inx_sub));
% y = squeeze(M1.DEC_SELCELL(inxcell,inx_comcont(2),inx_compori,inx_sub));
% D =[D; [x(:) y(:)]];
% 
% % 40-40 vs 100-40
% inx_sub2 = 17:22;
% inx_comcont= [2 4]
% x = squeeze(M1.DEC_SELCELL(inxcell,inx_comcont(1),inx_compori,inx_sub));
% y = squeeze(M1.DEC_SELCELL(inxcell,inx_comcont(2),inx_compori,inx_sub));
% D =[D; [x(:) y(:)]];
% 
% %% AN2
% 
% % 100-100 vs 40-100
% inx_comcont=  [1 4]
% x = squeeze(M2.DEC_SELCELL(inxcell,inx_comcont(1),inx_compori,inx_sub2));
% y = squeeze(M2.DEC_SELCELL(inxcell,inx_comcont(2),inx_compori,inx_sub2));
% D =[D; [x(:) y(:)]];
% 
% % 40-40 vs 100-40
% inx_comcont=  [2 3]
% x = squeeze(M2.DEC_SELCELL(inxcell,inx_comcont(1),inx_compori,inx_sub2));
% y = squeeze(M2.DEC_SELCELL(inxcell,inx_comcont(2),inx_compori,inx_sub2));
% D =[D; [x(:) y(:)]];
% 
% % 
% % (D(:,1)-D(:,2))*100
% p = ranksum(D(:,1),D(:,2))
% % 
% % figure; hold on;
% % plot(D(:,1),D(:,2),'.')
% % plot([min(D(:)) max(D(:))],[min(D(:)) max(D(:))]);
% 
% 
% %----- AWAKE
% inx_sub3=[23 25 26 27 29 30 32 33 ] %[23 25 26 27 29 30 32]%
% % inx_compori = 3;
% % inxcell = 3;
% % D=[];
% inx_comcont= [1 4]
% x = squeeze(M3.DEC_SELCELL(inxcell,inx_comcont(1),inx_compori,inx_sub3));
% y = squeeze(M3.DEC_SELCELL(inxcell,inx_comcont(2),inx_compori,inx_sub3));
% D =[D; [x(:) y(:)]];
% 
% inx_comcont= [2 3]
% x = squeeze(M3.DEC_SELCELL(inxcell,inx_comcont(1),inx_compori,inx_sub3));
% y = squeeze(M3.DEC_SELCELL(inxcell,inx_comcont(2),inx_compori,inx_sub3));
% D =[D; [x(:) y(:)]];
% 
% (D(:,1)-D(:,2))*100
% p = ranksum(D(:,1),D(:,2))
% figure; hold on;
% plot(D(:,1),D(:,2),'.')
% plot([min(D(:)) max(D(:))],[min(D(:)) max(D(:))]);
% 
% 
% 
% 
% figure; plot(D(:,1),D(:,2),'.');
% hold on; plot([min(D(:)) max(D(:))],[min(D(:)) max(D(:))],'k');
% 
% 
% 
% 
% 
