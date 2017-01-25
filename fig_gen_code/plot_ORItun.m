clear all
close all
D0 = load('Y:\data_2photon\Sedated\thy1_012916B1\04292016_thy1_0129B1_6\matlab\data\DISK_Cont-ORI_ctm0.60_HPF100_tau0.85_F0-10.0sigma.mat');
D0 = D0.data;
D = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN_FULLORI\thr5\DATA_DISK_UNION_CONTRSP_ctm0.60ses1.mat');
load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN_FULLORI\thr5\AN1_DISK_ORIsc_ctm0.60.mat');



%------ plot FOV
mimg = D0(1).meanimg;
mimg(mimg(:)>3000)=3000;
mimg(mimg(:)<300)=3800;
hf=figure; imagesc(mimg); colormap('gray'); axis image;
axis off;

pixel =D0(1).ROI{1}(1).pixelres
 1/(pixel*512/50)
 
 
 A= D0(1).meanimg(201:302,40:142);
 A(A(:)>3300)=3300;
 figure('Position',[100 100 102*1.5 102*1.5]); imagesc(A); colormap('gray'); axis image;
 axis off
%---------- plot tuning function
clr={'r', 'b'}; % red (100% contrast) vs blue (40% contrast)

ises=1
thr = 0.5;
Mf1 = ORIsc{ises,1};
a = squeeze(Mf1.as(1,:,:));  % scale  C40 = a*C100 +b
b = squeeze(Mf1.as(2,:,:));  % bias
Mf1.ev(a(:)<0)=0;
[Mev, mi ]= max(Mf1.ev,[],1);

mi = sub2ind(size(Mf1.ev),mi,1:size(Mf1.ev,2));

a = a(mi);
b = b(mi);


CLIST = find(Mev>thr);
nC= length(CLIST);
sc.a = a(CLIST);
B = squeeze(Mf1.Mresp(:,2,CLIST));
sc.b = b(CLIST)./mean(B,1);
sc.cell = CLIST;


Md1 = ORItun(ises,1);
% Md2 = ORItun(ises,2);
%     Md3 = ORItun(ises,3);

nCell = size(Md1.mresp,3);
for i = 1:nC
    icell = CLIST(i);
    if mod(i,35)==1,
        figure('Position',[680 173 1200 800]);
        kk=1;
    end
    subplot(5,7,kk);hold on
    
    M = Md1.mresp(:,1,icell);
    a = 1/max(M);
    Mn = M*a;
    S = Md1.semresp(:,1,icell);
    Sn = S*a;
    errorbar([0:30:330],Mn,Sn,'Color',clr{1},'linewidth',1.8);
    M = Md1.mresp(:,2,icell);    
    Mn = M*a;
    S = Md1.semresp(:,2,icell);
    Sn = S*a;
    errorbar([0:30:330],Mn,Sn,'Color',clr{2},'linewidth',1.8);
     
    
%     M = Md2.mresp(:,1,icell);    
%     Mn = M*a;
%     S = Md2.semresp(:,1,icell);
%     Sn = S*a;
%     errorbar([0:30:330],Mn,Sn,'Color',clr{1},'LineStyle','--');
%     M = Md2.mresp(:,2,icell);    
%     Mn = M*a;
%     S = Md2.semresp(:,2,icell);
%     Sn = S*a;
%     errorbar([0:30:330],Mn,Sn,'Color',clr{2},'LineStyle','--');
     
%     errorbar([0:30:330],Md2.mresp(:,1,icell),Md2.semresp(:,1,icell),'b--');
%     errorbar([0:30:330],Md2.mresp(:,2,icell),Md2.semresp(:,2,icell),'r--');
%     title(num2str(icell))
%         xlim([-15 110])


    xlim([-5 335])
    axis off
    kk = kk + 1;
end
% axis on
subplot(5,7,kk+2);hold on
plot([0 0],'r','linewidth',2);
plot([1 1],'b','linewidth',2);
legend('100%','40%')
set(gca,'FontSize',20)
%-------------- plot examples in an expanded view
figure('Position',[680 173 300 560]);
kk=1;
for i = [17 20]
    icell = CLIST(i);
 
    subplot(2,1,kk);hold on
    
    M = Md1.mresp(:,1,icell);
    a = 1/max(M);
    Mn = M*a;
    S = Md1.semresp(:,1,icell);
    Sn = S*a;
    errorbar([0:30:330],Mn,Sn,'Color',clr{1});
    M = Md1.mresp(:,2,icell);    
    Mn = M*a;
    S = Md1.semresp(:,2,icell);
    Sn = S*a;
    errorbar([0:30:330],Mn,Sn,'Color',clr{2});
     
    
    M = Md2.mresp(:,1,icell);    
    Mn = M*a;
    S = Md2.semresp(:,1,icell);
    Sn = S*a;
    errorbar(0:30:330,Mn,Sn,'Color',clr{1},'LineStyle','--');
    M = Md2.mresp(:,2,icell);    
    Mn = M*a;
    S = Md2.semresp(:,2,icell);
    Sn = S*a;
    errorbar(0:30:330,Mn,Sn,'Color',clr{2},'LineStyle','--');
     
%     errorbar([0:30:330],Md2.mresp(:,1,icell),Md2.semresp(:,1,icell),'b--');
%     errorbar([0:30:330],Md2.mresp(:,2,icell),Md2.semresp(:,2,icell),'r--');
%     title(num2str(icell))
%         xlim([-15 110])


    xlim([-5 335])
    axis off
    kk = kk + 1;
end
axis on
set(gca,'XTick',[0 180])
set(gca,'FontSize',20)

% old 
% %------------- plot scaling factor
% figure; hist(sc.a,[0.1:0.1:3]);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% set(gca,'FontSize',24,'XTick',[0 1 2]);
% xlim([0 2.5])
% figure; hist(sc.b,[-1:0.1:0.5])
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% set(gca,'FontSize',24,'XTick',[-1 0 1]);
% xlim([-1 1])
%---------------- plot 

b0 = 1; s=4; iscan=3
cellinx = D.cellinx_sel;
%cid = fliplr([2 6 11 14 22]);
cid = fliplr([1 4 15 17  27]);

t = 1:size(D0(iscan).dFF,1);
t = t*D0(iscan).Params.msperframe/1000;
inxt = find(t>200 & t<400);
t = t(inxt);
m0=0;
for i = 1: length(cid)
    ic = cellinx(CLIST(cid(i)));
    a = D0(iscan).Nhat(inxt,ic);
    m0 = max(max(a),m0);
end


figure('Position',[680 173 400 300]); hold on;
for i = 1: length(cid)
    ic = cellinx(CLIST(cid(i)));
    CLIST(cid(i)) 
    dff = D0(iscan).dFF(inxt,ic) + b0 +s;
    a = D0(iscan).Nhat(inxt,ic);
    a = a/m0*s*0.7 + b0; 
    plot(t,dff,'k');
    plot(t,a,'Color',[0.5 0.5 0.5]);
    b0 =  max(dff)+1;
end
xlim([200 400])
axis off


%-------------- mean contrast response
msf = 1000/D0(iscan).Params.msperframe;
spec.frames = -3 : floor(2*msf);
spec.dataType='Nhat';
spec.nCell= size(D0(iscan).Nhat,2);

[Y others]=data_sort(D0(3),spec,1);
evt = others.events{1};

inxC{1}  = find(mod(evt,2)==0); %index for 100% contrast
inxC{2} = find(mod(evt,2)==1); % index for 40% contrast
figure; hold on;
fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','w','Facealpha',.9);
fill_rec([0 0.5],[2.8 2.8],[3 3],'k')
fill_rec = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor','w','Facealpha',.5);
fill_rec([0.18 0.68],[3.3 3.3],[7 7],[1 1 0])

for ii = 1:2
    X= squeeze(mean(Y{1}(inxC{ii},:,cellinx),1));
    mX = mean(X,2);
    sX= (std(X,0,2)/sqrt(spec.nCell));
    errorbar(spec.frames/msf,mX*100,sX*100, 'Color',clr{ii},'LineWidth',2)
end
set(gca,'FontSize',24,'XTick',[-0.5:0.5:2]);
xlim([-0.55 2.05])
ylim([2 8])
% axis off



%% plot scales and bias distribution


load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN_FULLORI\thr5\AN1_DISK_ORIsc_ctm0.60_fit_ab4.mat');

ises=1
thr = 0.5;
Mf1 = ORIsc{ises,1};
a = squeeze(Mf1.as(1,1,:));  % scale  C40 = a*C100 +b
b = squeeze(Mf1.as(2,1,:));  % bias
Mev = Mf1.ev;




CLIST = find(Mev>thr & a'<2);
nC= length(CLIST);
sc.mev = Mev(CLIST);
sc.a = a(CLIST);
B = squeeze(Mf1.Mresp(:,2,CLIST));
sc.b = b(CLIST)./mean(B,1)';
sc.b0= b(CLIST);
sc.cell = CLIST;
sc.M = Mf1.Mresp(:,:,CLIST);

figure; hist(sc.a,[0.1:0.1:2]);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
set(gca,'FontSize',24,'XTick',[0 1 2]);
xlim([0 2])
box off
figure; hist(sc.b,[-1:0.1:1])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
set(gca,'FontSize',24,'XTick',[-1 0 1]);
xlim([-1 1])
box off

%--- case 1
inx =find(sc.a<0.8)
%cellist =1:4;
cellist =[2 4];
%------case 2
inx =find(sc.a>1)
%cellist=[1 3 4 5]
cellist =[1 3]


figure('position',[680   678   280   360]);
for i0=1:2
   i = cellist(i0); 
subplot(2,1,i0);
hold on;
r40 = sc.M(:,1,inx(i));
r100 = sc.M(:,2,inx(i));
n =max(r100);
nr100 = r100/n;
nr40 = r40/n;
a = sc.a(inx(i));
b = sc.b0(inx(i));
plot(0:30:330,nr40,'b.-','LineWidth',2);
plot(0:30:330,nr100,'r.-','LineWidth',2);
plot(0:30:330,(a*r100+b)/n,'k--','LineWidth',1.5);
% tltstr = sprintf('a: %.2f, b: %.2f, EV: %.2f',sc.a(inx(i)),sc.b(inx(i)),sc.mev(inx(i)));
% tltstr = sprintf('a: %.2f, b: %.2f',sc.a(inx(i)),sc.b(inx(i)));
% title(tltstr,'FontSize',15);
xlim([0 330]);ylim([0 1.5]);
set(gca,'FontSize',20);
set(gca,'XTick',[0 180])
set(gca,'YTick',[0 1])
axis off
end

