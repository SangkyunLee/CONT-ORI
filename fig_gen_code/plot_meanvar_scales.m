clear all

D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_WSWC_AN1-16_ORIsc_ctm0.60_fit_ab4.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_WSWC_AN17-22_ORIsc_ctm0.60_fit_ab4.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\L_WSWC_AW23-40_ORIsc_ctm0.60_fit_ab4.mat');

M(1)=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_AN1-16_ORIsc_ctm0.60_fit_ab4.mat');
M(2)=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\L_AN17-22_ORIsc_ctm0.60_fit_ab4.mat');
M(3)=load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\L_AW23-40_ORIsc_ctm0.60_fit_ab4.mat');


parinx =1 % scale
% parinx =2 % bias


thr = 0.5
thr_scale=1

A = cell(28,2);
B = cell(28,1);
varA = NaN*ones(28,2);
varB = NaN*ones(28,1);
for icont = 1 : 2
    j = 1;
    for iexp = 1 : 3
        nsub = size(D(iexp).ORIsc,1);
        for isub = 1 : nsub

            X0 = D(iexp).ORIsc{isub,icont};
            if icont==1,
                X1 = M(iexp).ORIsc{isub};
            end
            if isempty(X0)
                continue;
            end

            %inx = find(X0.ev>thr & squeeze(X0.as(1,1,:))'<thr_scale);
            inx = find(X0.ev>thr);
            p = squeeze(X0.as(parinx,:,inx));
            if parinx==2
                maxR = max(squeeze(D(iexp).ORItun(isub).mresp(:,icont,:)))';
                p= p./maxR(inx);
            end
            A{j,icont} = p;
            varA(j,icont) = var(p);
            
            if icont==1,
                inx = find(X1.ev>thr & squeeze(X1.as(1,1,:))'<thr_scale);
                p = squeeze(X1.as(parinx,:,inx));
                if parinx==2
                    maxR = max(squeeze(M(iexp).ORItun(isub).mresp(:,1,:)))';
                    p= p./maxR(inx);
                end
                
                B{j} = p;
                varB(j) = var(p);
            end
            j = j+1;
        end
    end
end
        

varA(varA(:)==0)=NaN;
varB(varB(:)==0)=NaN;



X =[varA varB];
nX = sum(~isnan(X))

hf=figure('Position',[680   500   560   420]);
model_series = nanmean(X);
model_error = nanstd(X,0,1)./sqrt(nX);
h = bar(model_series);
set(h,'BarWidth',0.4,'FaceColor',[0.91 0.91 0.91]);    % The bars will now touch each other
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTicklabel',{'f_{100}-f_{100}','f_{40}-f_{40}','f_{100}-f_{40}'})
hold on;
errorbar( model_series, model_error, 'k', 'linestyle', 'none','linewidth',2);

set(gca,'FontSize',22,'linewidth',2);
xlim([0.5 3.5])

box off

%[p,table,stats]  = anova1(X)
[p,table,stats] = kruskalwallis(X)
[c,m,h,nms] = multcompare(stats,'alpha',1e-7);


