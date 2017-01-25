clear all
% load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE\thr5\NEW1_ORIsc_ctm0.60.mat')
% seslist=11:17;
load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\DISK_ORIsc_ctm0.60.mat')
seslist=[(1:8) 11 12];


for ises = 3%seslist
Mf1 = ORIsc{ises,2};
Mf2 = ORIsc{ises,2};
Md = ORItun(ises);

a=squeeze(Mf1.as(1,:,:));
Mf1.ev(a(:)<0)=0;
[Mev1, mi1 ]= max(Mf1.ev,[],1);
a=squeeze(Mf2.as(1,:,:));
Mf2.ev(a(:)<0)=0;
[Mev2, mi2 ]= max(Mf2.ev,[],1);

Mev=[Mev1];
ratio(:,ises)=[length(find(max(Mev,[],1)>0.4)) length(find(max(Mev,[],1)>0.4))/length(Mev) ];
end

icell=63
figure;
A=Mf1.as(:,mi1(icell),icell);
X=Mf1.Mresp(:,:,icell);
y = A(1)*X(:,Mf1.conref)+A(2);
subplot(121); plot([X(:,Mf1.conref) X(:,Mf1.concom)]);
hold on;
plot(y,'k--')

A=Mf2.as(:,mi2(icell),icell);
X=Mf2.Mresp(:,:,icell);
y = A(1)*X(:,Mf2.conref)+A(2);
subplot(122); plot([X(:,Mf2.conref) X(:,Mf2.concom)]);
hold on;
plot(y,'k--')

%--------------
hOSI={};
for ises=seslist
    Mf1 = ORIsc{ises,1};
    Mf2 = ORIsc{ises,2};
    Md = ORItun(ises);


    figure;
    for icont=1:3
        Ncell =size(Md.mresp,3)
        OSI1= zeros(Ncell,1);
        ori= Md.evtord(:,1,2);
        ori =[ ori; ori+180];
        for icell =1:Ncell    
            x = Md.mresp(:,icont,icell);    
            K=exp(1i*ori/180*pi);
            x1 = mean(x)*ones(size(ori));
            x1([1 :4])=x;
            OSI1(icell) = abs(sum(K.*x1))/(sum(x1));
        end
        subplot(3,1,icont);
        hist(OSI1,[0:0.01:0.3]);
        
       inx=find(OSI1>0.1);
       hOSI{ises,icont}=inx(:)';
    end
end


%% check SNR
SNR={};
Nr=3; Nc=4;
for ises=seslist(9)
    Md = ORItun(ises);
    Ncell =size(Md.mresp,3);
    selcells = unique(cell2mat(hOSI(ises,:)));
    K=Md.mresp./Md.stdresp;
    K1 = K(:,:,selcells);
    SNR{ises}=K;
    figure; hist(K1(:),100);
    md = median(K1(:));
    stdsnr = std(K1(:));
    title(sprintf('Median SNR:%.2f +/- %.2f',md,stdsnr));
    kk=1;
    for icell =selcells
        if mod(kk,Nr*Nc)==1,
            kk=1;
            figure;
        end
        subplot(Nr,Nc,kk);
        errorbar(Md.evtord(:,:,2),Md.mresp(:,:,icell),Md.stdresp(:,:,icell));
        title(num2str(icell));
        ylim([0 0.1])
        kk = kk +1;
    end
end
%%

%%


 K=exp(1i*ori/180*pi);
% a=abs(sum(K.*mean(x)))/(mean(x)*24)


% icont=1
% Ncell =size(Md.mresp,3)
% OSI= zeros(Ncell,1);
% ori= Md.evtord(:,1,2);
% ori=-15:15:330;
% for icell =1:Ncell    
%     x = Md.mresp(:,icont,icell);    
%     K=exp(1i*ori/180*pi*2);
%     x1 = mean(x)*ones(size(ori));
%     x1([1 2 4 8])=x;
%     OSI(icell) = abs(sum(K.*x1))/(sum(x1));
% end
% 
% 
% 
