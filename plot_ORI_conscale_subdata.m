clear all



for ises = 17:22%seslist
Mf1 = ORIsc{ises,2};

a=squeeze(Mf1.as(1,:,:));
Mf1.ev(a(:)<0)=0;
[Mev1, mi1 ]= max(Mf1.ev,[],1);

Mev=[Mev1];
ratio(:,ises)=[length(find(max(Mev,[],1)>0.4)) length(find(max(Mev,[],1)>0.4))/length(Mev) ];
end


hOSI={};
for ises = 1 %seslist   
    Md1 = ORItun(ises,1);
    Md2 = ORItun(ises,2);
%     Md3 = ORItun(ises,3);

    nCell = size(Md1.mresp,3);
    for icell = 1 :nCell
        if mod(icell,36)==1,
            figure;
            kk=1;
        end
        subplot(6,6,kk);hold on
        errorbar([-10 0 30 90],Md1.mresp(:,1,icell),Md1.semresp(:,1,icell),'b');
        errorbar([-10 0 30 90],Md1.mresp(:,2,icell),Md1.semresp(:,2,icell),'r');
        errorbar([-10 0 30 90],Md2.mresp(:,1,icell),Md2.semresp(:,1,icell),'b--');
        errorbar([-10 0 30 90],Md2.mresp(:,2,icell),Md2.semresp(:,2,icell),'r--');
        
%         errorbar([0:30:330],Md1.mresp(:,1,icell),Md1.semresp(:,1,icell),'b');
%         errorbar([0:30:330],Md1.mresp(:,2,icell),Md1.semresp(:,2,icell),'r');
%         errorbar([0:30:330],Md2.mresp(:,1,icell),Md2.semresp(:,1,icell),'b--');
%         errorbar([0:30:330],Md2.mresp(:,2,icell),Md2.semresp(:,2,icell),'r--');
        title(num2str(icell))
        xlim([-15 110])
%         xlim([-5 335])
%         plot([0:30:330],Md2.mresp(:,1,icell),'bx');
%         plot([0:30:330],Md2.mresp(:,2,icell),'rx');
        kk = kk + 1;
    end
end