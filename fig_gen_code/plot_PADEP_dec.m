clear all
close all;

inxcell = 3;


M{1} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\PADEP_AN1-16SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat');
M{2} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AN\thr5\PADEP_AN17-22SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')
M{3} = load('Z:\data_2photon\matlab_2ndLev\NEW_DECODING_NOBIAS_ZMEAN\AWAKE_EYE\thr5_eyethr_xy1_p1\PADEP_SUBCELL-Xsel_CRSCON_SMLR_L2_ctm0.60.mat')





nc = 16%size(M{1}.comcont,2);
D = cell(6,8);
Do = cell(6,8);



dispord =[6 3 5 1 2 4];
ncellord = [ (2:6) 1];
ncell = [NaN 1 3 5 10 20];
for inxcell0 = 1 : 6
    inxcell = ncellord(inxcell0);
    for ic = 1 : nc    
        D1 =cell(length(M),6);
        IX =cell(length(M),6);
        for ises = 1 : 3%length(M)
            % collect ORI decoding accs
            for iori = 1:6
                inx_compori = dispord(iori);
                x = squeeze(M{ises}.DEC_SELCELL(inxcell,ic,inx_compori,:));                
                D1{ises,iori} = x(x>0); 
            end             
        end % ises
        D{inxcell0, ic}= cell2mat(D1);
        Do{inxcell0, ic}= mean(cell2mat(D1),2);
        
    end % ic % contrast    
end %inxcell0

%% ------------------
% 
% comcont={{'100,H','100,H'},{'100,H','100,L'},{'100,H','40,H'},{'100,H','40,L'},...
%     {'40,H','100,H'},{'40,H','100,L'},{'40,H','40,H'},{'40,H','40,L'},...
%     {'100,L','100,H'},{'100,L','100,L'},{'100,L','40,H'},{'100,L','40,L'},...
%     {'40,L','100,H'},{'40,L','100,L'},{'40,L','40,H'},{'40,L','40,L'}};

%%
% D1=cell2mat(Do);
% 
% inx1=8+2
% inx2=4+2
% incell=6
% x = Do{incell,inx1};
% y = Do{incell,inx2};
% p=signrank(x(:),y(:));
% p=p*8
% figure; 
% plot(x,y,'.'); hold on
% plot([0.45 1],[0.45 1],'k','LineWidth',2)
% [p,table,stats] = anova2([x;y],28);
% [c,m,h,nms] = multcompare(stats);

%------------------------
close all

incell = 5 % incell=6: including all cells
ps = NaN*ones(4,4);
macc = cell(4);
for mode =[1 2 3 4]

    switch mode
        case {1}
            inx_set0 = [0 4 8 12]/4 +1; 
            inx_set = 1+[0 4 8 12]; %100H
        case {2}
            inx_set0 =[8 0 4 12]/4 + 1; %100L
            inx_set = 2 + [8 0 4 12]; %100L
        case {3}
            inx_set0 = [4 0 8 12]/4 + 1; %40H
            inx_set = 3 + [4 0 8 12]; %40H        
        
        case {4}
            inx_set0 = [12 0 4 8]/4 + 1; %40L
            inx_set = 4 + [12 0 4 8]; %40L
    end

    figure; hold on;
    clrs={'r','g','b','y'}
    x=cell2mat(Do(incell,inx_set(1)));
    macc{mode,1}=x;

    bfc =3; % bonferroni correction
    %lgdstr = cell(1,4);
    for i = 2 : 4    
        y = cell2mat(Do(incell,inx_set(i)));
        ps(mode,inx_set0(i)) = signrank(x(:),y(:))*bfc;
        
        plot(x,y,'.','Color',clrs{inx_set0(i)},'MarkerSize',20);         
        macc{mode,i}=y;
        %lgdstr{i}=M{1}.comcont{inx_set(i)}{1};
    end
    plot([0.45 1],[0.45 1],'k','LineWidth',2)

    set(gca,'FontSize',20);
    set(gca,'XTick',[0.6 0.8 1]);
    set(gca,'YTick',[0.6 0.8 1]);
    box off;
    axis equal;
    xlim([0.45 1]);
    ylim([0.45 1]);
    %legend(lgdstr,'Location','northwest')
    % legend boxoff
end



ps1=ps(:,[1 3 2 4]);% to align ps to 100,H|100,L|40,H|40,L
% ps1(ps1(:)>0.05)=Inf;
% 
% psall(:,:,incell)=ps1;


close all

for mode = 1 : 4
    
    switch mode
        case {1}
            
            %inx_set = 1+[0 4 8 12]; %100H
            inx_set = 1+[0 8 4 12]; %100H
            in2 = [1 3 2 4];
            
        case {2}
            %inx_set = 2 + [8 0 4 12]; %100L
            inx_set = 2 + [8 0 4 12]; %100L
            in2 = [1 2 3 4];
        case {3}
            %inx_set = 3 + [4 0 8 12]; %40H
            inx_set = 3 + [4 12 0 8 ]; %40H
            in2 = [1 4 2 3];
        case {4}
            %inx_set = 4 + [12 0 4 8]; %40L
            inx_set = 4 + [12 4 0 8]; %40L
            in2 = [1 3 2 4];
    end
    
    strx = M{1}.comcont(inx_set)
    str = cell(1,4);
    for i=1:4
    str{i}=strx{i}{1};
    end
    
    x = cell2mat(macc(mode,in2));
    figure; hold on;
    bar(mean(x,1),'FaceColor',[0.91 0.91 0.91],'EdgeColor','k','BarWidth',0.5);
    errorbar(mean(x),std(x,0,1)/sqrt(size(x,1)),'k.','LineWidth',2);
    ylim([0.5 1])
    xlim([0.5 4.5])
    set(gca,'XTick',1:4);
    set(gca,'XTickLabel',str);
    set(gca,'YTick',[0.6 0.8 1]);
    set(gca,'FontSize',20);
    
end




%% plot decoding accuracy difference  between contrast-and-PA-dependent models and others in n-cell

close all


for mode =[1 2 3 4]


    switch mode
        case {1}
            inx_set0 = [0 4 8 12]/4 +1; 
            inx_set = 1+[0 4 8 12]; %100H
        case {2}
            inx_set0 =[8 0 4 12]/4 + 1; %100L
            inx_set = 2 + [8 0 4 12]; %100L
        case {3}
            inx_set0 = [4 0 8 12]/4 + 1; %40H
            inx_set = 3 + [4 0 8 12]; %40H        
        
        case {4}
            inx_set0 = [12 0 4 8]/4 + 1; %40L
            inx_set = 4 + [12 0 4 8]; %40L
    end
    a=cell(5,3);
    for incell = 1:5 % incell=6: including all cells
        x=cell2mat(Do(incell,inx_set(1)));
        for i = 2 : 4    
            y = cell2mat(Do(incell,inx_set(i)));
            a{incell,i-1} = (y-x)./x*100;
        end
    end
    strx = M{1}.comcont(inx_set);
    str = cell(1,3);
    for i=2:4
        str{i-1}=strx{i}{1};
    end
    
    
    
%     [p,table,stats]  =anova2(cell2mat(a),28);
%     [c,m,h,nms] = multcompare(stats,'estimate','column','alpha',1e-15);
    
    clrs={'r','g','b','y'}
    sig = cellfun(@mean,a,'UniformOutput',false);
    sig = cell2mat(sig);
    nrep = size(a{1,1},1);
    calerr = @(x) sqrt(std(x,0,1)/nrep);
    err = cellfun(calerr,a);
%     err = cell2mat(err);
    figure; hold on;
    for i=1:3
        errorbar([1 3 5 10 20],sig(:,i),err(:,i),'linewidth',2,'Color',clrs{inx_set0(i+1)});
    end
    set(gca,'FontSize',20);
%     axis equal
    xlim([-1 21])
    ylim([-15 5])
end

%generation of figure legend
% figure; hold on
% for i=1:4
%     plot(i,i,'linewidth',2,'Color',clrs{i});
% end
% legend({'100,H','40,H','100,L','40,L'},'Location','northwest')
% set(gca,'FontSize',20);


%% plot mean population activity depending on contrast and pouplation levels
D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN1-16-Xsel_ctm0.60.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN17-22-Xsel_ctm0.60.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\PADEPDATA_-Xsel_ctm0.60.mat')

%%------ check individual session
% for iexp = 1 : 3
%     ns = size(D(iexp).DPA,3);
%     for isub =1 :ns
%         A=D(iexp).DPA(:,:,isub);
%         if ~isempty(A{1,1})
%             A1=cellfun(@mean,A,'UniformOutput',false);
%             nrep = length(A1{1,1});
%             A1 =cell2mat(A1);
%             [p,table,stats] =anova2(A1',nrep);
%             [c,m,h,nms] = multcompare(stats,'estimate','column');
%              
%         end
%     end
% end


%------ compare the mean response
K = zeros(4,4,28);

j = 1;
for iexp = 1 : 3
    ns = size(D(iexp).DPA,3);
    for isub =1 :ns
        A=D(iexp).DPA(:,:,isub);
        if ~isempty(A{1,1})
            A1 = cellfun(@mean,A,'UniformOutput',false);
            A1 = cellfun(@mean,A1);
            K(:,:,j) = A1;
            j = j + 1;
             
        end
    end
end



%----------  anova2
K = K*100;


figure;
model_series = mean(K,3);
model_error = std(K,0,3)/sqrt(size(K,3));
h = bar(model_series);
set(h,'BarWidth',1);    % The bars will now touch each other
set(gca,'XTicklabel','100,H|100,L|40,H|40,L')
% lh = legend('-15/-10','0','30','90');

hold on;
numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','linewidth',2);
end
set(gca,'FontSize',20);
xlim([0.5 4.5])
box off

%---- get statistics
K1= permute(K,[3 2 1]);
K1 = reshape(K1,[size(K1,1)*size(K1,2) size(K1,3)]);
[p,table,stats]  = anova2(K1,28)
[c,m,h,nms] = multcompare(stats,'estimate','column','alpha',5e-3);




%% plot SNR of response depending on contrast and pouplation levels
clear all;
D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN1-16-Xsel_ctm0.60.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN17-22-Xsel_ctm0.60.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\PADEPDATA_-Xsel_ctm0.60.mat')



%------ compare the mean response
FANO = zeros(4,4,28);
SNR = zeros(4,4,28);
var2 = @(x) var(x,0,2);
mean2 = @(x) mean(x,2);
fano = @(x) var(x,0,2)./mean(x,2);
snr = @(x) mean(x,2)./sqrt(var(x,0,2));

j = 1;
for iexp = 1 : 3
    ns = size(D(iexp).DPA,3);
    for isub =1 :ns
        A=D(iexp).DPA(:,:,isub);
        A = cellfun(@times,A,repmat({100},size(A)),'UniformOutput', false);% corresponding to percentage
        
        if ~isempty(A{1,1})
            FANO1 = cellfun(fano,A,'UniformOutput',false);
            SNR1 = cellfun(snr,A,'UniformOutput',false);

            
            A1 = cellfun(@mean,FANO1);
            FANO(:,:,j) = A1;
            
            A1 = cellfun(@mean,SNR1);
            SNR(:,:,j) = A1;
            
            j = j + 1;
             
        end
    end
end



%----------  anova2

% K=FANO;
K=SNR;

figure;
model_series = mean(K,3);
model_error = std(K,0,3)/sqrt(size(K,3));
h = bar(model_series);
set(h,'BarWidth',1);    % The bars will now touch each other
set(gca,'XTicklabel','100,H|100,L|40,H|40,L')
% lh = legend('-15/-10','0','30','90');

hold on;
numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','linewidth',2);
end
set(gca,'FontSize',20);
xlim([0.5 4.5])
box off
ylim([0.6 1])
%---- get statistics
K1= permute(K,[3 2 1]);
K1 = reshape(K1,[size(K1,1)*size(K1,2) size(K1,3)]);
[p,table,stats]  = anova2(K1,28)
[c,m,h,nms] = multcompare(stats,'estimate','column','alpha',5e-2);



%% measure angle between CPRS(Contrasts-and-population-response-specific) vs non-CPRS
clear all; close all;
D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN1-16-Xsel_ctm0.60.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN17-22-Xsel_ctm0.60.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\PADEPDATA_-Xsel_ctm0.60.mat')
  
compstr = cell(6,1);
mean2 =@(x)mean(x,2);
A = cell(40,1);
k=1;
for iexp = 1 :3
    nsub = size(D(iexp).DPA,3);
    M = cellfun(mean2,D(iexp).DPA,'UniformOutput',false);
    
    
    for isub = 1 : nsub
        if ~isempty(M{1,1,isub})
            dsize = size(M(:,:,isub));
            csize =size(M{1,1,isub},1);
            J = cell2mat(M(:,:,isub));
            J = reshape(J,[csize dsize]);
            
            A0 = zeros(6,4);
            
            for iori = 1 : 4
                j0=1;                
                for iCPR = 1 : 4               
                    J1 = J(:,iCPR,iori);
                    for iCPR2 = iCPR+1 : 4
                        J2 = J(:,iCPR2,iori);
                        A0(j0,iori) = acos(J1'*J2/(sqrt(J1'*J1)*sqrt(J2'*J2)))/pi*180;
                        if iori==1&& k==1
                            compstr{j0}= [D(1).comcont{iCPR} '-' D(1).comcont{iCPR2}];
                        end
                        
                        j0 = j0 +1;
                        
                        
                    end                
                end
            end
            A{k}=A0;
            k = k+1;
        end
    end
end
    

submean = @(x)[mean(x([1 6],:),2); mean(mean(x(2:5,:),1),2)];

B= cellfun(submean,A(1:28),'UniformOutput',false);

B= cell2mat(B')';

%------- get statistics

[p,table,stats]  = anova1(B)
[c,m,h,nms] = multcompare(stats,'alpha',1e-2);



hf=figure('Position',[680   678   640   420]);
model_series =mean(B,1);
model_error = std(B,0,1)/sqrt(size(B,1));
h = bar(model_series);
set(h,'BarWidth',0.4,'FaceColor',[0.91 0.91 0.91]);    % The bars will now touch each other
set(gca,'XTicklabel','100,H-100,L|40,H-40,L|Crs. Cont.')
hold on;
errorbar( model_series, model_error, 'k', 'linestyle', 'none','linewidth',2);

set(gca,'FontSize',20);
xlim([0.5 3.5])
set(gca,'YTick',[0 10 20 30])
box off




%% measure angle between ORIs

D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN1-16-Xsel_ctm0.60.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\PADEPDATA_AN17-22-Xsel_ctm0.60.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\PADEPDATA_-Xsel_ctm0.60.mat')

ori_list = [-10 0 30 90];
oridif = zeros(6,1);

mean2 =@(x)mean(x,2);
A = cell(40,1);
k=1;
for iexp = 1 :3
    nsub = size(D(iexp).DPA,3);
    M = cellfun(mean2,D(iexp).DPA,'UniformOutput',false);
    
    
    for isub = 1 : nsub
        if ~isempty(M{1,1,isub})
            dsize = size(M(:,:,isub));
            csize =size(M{1,1,isub},1);
            J = cell2mat(M(:,:,isub));
            J = reshape(J,[csize dsize]);
            
            A0 = zeros(6,4);
            
            for iCPR  = 1 : 4
                j0=1;                
                for io1 = 1 : 4               
                    J1 = J(:,iCPR,io1);
                    for io2 = io1+1 : 4
                        J2 = J(:,iCPR,io2);
                        A0(j0,iCPR) = acos(J1'*J2/(sqrt(J1'*J1)*sqrt(J2'*J2)))/pi*180;
                        if iCPR==1&& k==1
                            oridif(j0)= abs(ori_list(io1)-ori_list(io2));
                        end
                        
                        j0 = j0 +1;
                        
                        
                    end                
                end
            end
            A{k}=A0;
            k = k+1;
        end
    end
end
    

B=A(1:28,:);
B= reshape(cell2mat(B),[6 28 4]);
Bt = permute(B,[2 1 3]); 
B1= permute(B,[1 3 2]);
[~,newinx] =sort(oridif,'ascend');
Bt =Bt(:,newinx,:);
B1 =B1(newinx,:,:);


%------- get statistics

Bt=mean(Bt,3);
[p,table,stats]  = anova1(Bt)
[c,m,h,nms] = multcompare(stats,'alpha',0.05);

%--- plot
Bt=mean(Bt,3);
hf=figure('Position',[680   678   720  420]);
model_series =mean(Bt,1);
model_error = std(Bt,0,1)/sqrt(size(Bt,1));
h = bar(model_series);
set(h,'BarWidth',0.5,'FaceColor',[0.91 0.91 0.91]);    % The bars will now touch each other
set(gca,'XTicklabel','10(15)|30|40(45)|60|90|100(105)')
hold on;
errorbar( model_series, model_error, 'k', 'linestyle', 'none','linewidth',2);

set(gca,'FontSize',20);
xlim([0.5 6.5])
box off



%------- plot indivdidual cases
% %------- get statistics
% 
% Bt=reshape(Bt, [size(Bt,1)*size(Bt,2) size(Bt,3)]);
% [p,table,stats]  = anova2(Bt,28)
% [c,m,h,nms] = multcompare(stats,'estimate','row','alpha',1e-10);
% 
% %--- plot
% 
% hf=figure('Position',[680   678   860   420]);
% model_series = mean(B1,3);
% model_error = std(B1,0,3)/sqrt(size(B1,3));
% h = bar(model_series);
% set(h,'BarWidth',1);    % The bars will now touch each other
% set(gca,'XTicklabel','10(15)|30|40(45)|60|90|100(105)')
% % lh = legend('100,H','100,L','40,H','40,L');
% 
% hold on;
% numgroups = size(model_series, 1); 
% numbars = size(model_series, 2); 
% groupwidth = min(0.8, numbars/(numbars+1.5));
% for i = 1:numbars
%       % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
%       x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
%       errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','linewidth',2);
% end
% set(gca,'FontSize',20);
% xlim([0.5 6.5])
% box off
% clrs={[0.3 0.3 0.3],[0.5 0.5 0.5],[0.7 0.7 0.7],[0.9 0.9 0.9]}
% for i=1:4
% set(h(i),'facecolor',clrs{i})
% end



%%  ORI scale between different population activity vs between contrast


clear all; close all
D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_AN1-16_PADEP_ORIsc_ctm0.60_fit_ab4.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_AN17-22_PADEP_ORIsc_ctm0.60_fit_ab4.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\L_AW23-40_PADEP_ORIsc_ctm0.60_fit_ab4.mat')
thr=0.5
k=1;
evs=cell(6,28);
A=cell(6,28);
ixcont=[2 2 2 1 2 1]
for iexp = 1 : 3
    M=D(iexp).Out.ORIsc;
    nsub = size(M,2);
    ncx = size(D(iexp).Out.cx,2);
    for isub = 1 : nsub
        for icx = 1 : ncx
            if ~isempty(M{icx,isub})                
                evs{icx,k} = M{icx,isub}.ev;
                b = squeeze(mean(M{icx,isub}.Mresp(:,ixcont(icx),:),1))';
                b = [ones(1,size(b,2)); b];
                A{icx,k} = squeeze(M{icx,isub}.as)./b;
                
            end
        end
        if ~isempty(M{icx,isub})
            k = k + 1;
        end
    end
end

plotsig =1;%scale
plotsig =2;%bias

athr=10
spread = NaN*ones(6,28);
ncell = NaN*ones(1,28);
for isub = 1 : 28
    E1 = cell2mat(evs(:,isub));   
    A1 = cell2mat(A(:,isub));
    a=A1(1:2:end,:);
    
    inx1 = E1(1:2,:)>thr;
    inx2 = ~(sum(a>athr)>0);
    
    inx = (sum(inx1,1)==2) & inx2;
    if sum(inx)>1
        ncell(isub)=sum(inx);
        for j = 1 : 6
            if plotsig==2
                spread(j,isub) = var(A{j,isub}(plotsig,inx));
            else
                spread(j,isub) = var(A{j,isub}(plotsig,inx));
            end
        end
    end
end
spread1 = [spread([1 2],:); mean(spread(3:end,:),1)];
        
hf=figure('Position',[680   578   640   420]);

model_series =nanmean(spread1,2);
model_error = nanstd(spread1,0,2)/sqrt(sum(~isnan(ncell)));
h = bar(model_series);
set(h,'BarWidth',0.4,'FaceColor',[0.91 0.91 0.91]);    % The bars will now touch each other
set(gca,'XTicklabel','100,H-100,L|40,H-40,L|Crs. Cont.')
hold on;
errorbar( model_series, model_error, 'k', 'linestyle', 'none','linewidth',2);
set(gca,'FontSize',20);
xlim([0.5 3.5])
box off        

[p,table,stats]  = anova1(spread1')
[c,m,h,nms] = multcompare(stats,'estimate','row','alpha',1e-4);
% 


%% count franction of cells with ev>0.5
frac=zeros(2,28);
for ii=1:28
    for j=1:2
        ev = evs{j,ii};
        frac(j,ii) = length(find(ev>0.5))/length(ev);
    end
end
    
    


