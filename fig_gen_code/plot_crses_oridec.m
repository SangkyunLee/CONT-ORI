% test Cross-Contrast,Time versus Contrast-specific,Time effect


clear all
close all;



%datapath='/home/slee/data/data_2photon/matlab_2ndLev/NEW_DECODING_NOBIAS_ZMEAN'
datapath = 'Z:/data_2photon/matlab_2ndLev/NEW_DECODING_NOBIAS_ZMEAN'

M{1,1} = load(fullfile(datapath,'/AN/thr5/P1-P2_AN17-22CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'))
M{2,1} = load(fullfile(datapath,'AWAKE_EYE/thr5_eyethr_xy1_p1/P1-P2_AW23-40CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'));

M{1,2} = load(fullfile(datapath,'/AN/thr5/P1-P2_AN17-22CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'))
M{2,2} = load(fullfile(datapath,'AWAKE_EYE/thr5_eyethr_xy1_p1/P2-P1_AW23-40CELLSEL_SMLRW_CRSCON_SMLR_L2_ctm0.60.mat'));





nc= 2;
D1 = cell(6,6);


dispord =[6 3 5 1 2 4];
ncellord = [ (2:6) 1];
ncell = [NaN 1 3 5 10 20];


for inxcell0 = 1 : 6
    inxcell = ncellord(inxcell0);
    for iori = 1 : 6
        inx_compori = dispord(iori);

        for ises = 1 : size(M,1)

            % collect ORI decoding accs
            nsub = size(M{ises,1}.DEC_SELCELL,2);
            D0 = cell(nsub,1);
            nseg = size(M,2);
            for isub = 1: nsub
                x=[];
                for k = 1:nseg
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
                    if k==1,
                        x = x2;
                    else
                        x = x+x2;
                    end
                end
                x = x/nseg;

                if ~isempty(x)
                    D0{isub}=x;
                end
            end
            D1{inxcell0,iori} = [D1{inxcell0,iori}; cell2mat(D0)];
        end % ises
    end % iori
end %inxcell0

Y= cell2mat(D1(end,:)');
Y = reshape(Y,[14 6 4]);

%-------------------------
%100,100 - 40,100
x=Y(:,:,1);y=Y(:,:,3);


%40,40 - 100,40
x=Y(:,:,4);y=Y(:,:,2);
%-----------------------
%100,100 - 100,40
x=Y(:,:,1);y=Y(:,:,2);

%40,40 - 40,100
x=Y(:,:,4);y=Y(:,:,3);
%----------------------------



figure; plot(x,y,'.','MarkerSize',20);
lgdstr={'60','90','30','100/105','40/45','10/15'}
% [legend_h,object_h,plot_h,text_strings] =legend(lgdstr(dispord),'color','none');
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

p=signrank(x(:),y(:))

 [p,table,stats] = anova2([x;y],14);
 [c,m,h,nms] = multcompare(stats,'estimate','row');
% [p,table,stats] = anova2(,14);
% [c,m,h,nms] = multcompare(stats,'estimate','column');
