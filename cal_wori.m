function Out=cal_wori(X0,evt,calfun)



    inxsample = evt.inxsample;
    cons = evt.cons;
    oris = evt.oris;
    
    
    lcon = unique(cons);    
    conref = find(lcon==evt.conref);
    concom = find(lcon==evt.concom);
    
    N = min(cellfun(@length,inxsample));
    O = length(unique(oris));
    opts.O = O;
    opts.N = N;
    T = N*O;
    D = zeros(T,size(X0,2), length(lcon));
    M=zeros(O,length(lcon),size(X0,2));    
    

    opts.Ncv = 1000;
    if strcmp(calfun,'fit_ab') || strcmp(calfun,'fit_ab2') || strcmp(calfun,'fit_ab4')
        opts.lams =  fliplr(10.^(linspace(-5,1,30)))*T;
%         opts.lams =  (10.^(linspace(-5,1,20)));
    end

    
    for ic = 1 : length(lcon)
        iorg = find(cons==lcon(ic));
        fprintf('CON=%d,ORI=%s\n',lcon(ic),num2str(oris(iorg)));
        for ii = 1: length(iorg)
            J = inxsample{iorg(ii)};
            K = randperm(length(J));
            D((ii-1)*N+(1:N),:,ic) = X0(J(K(1:N)),:);
            M(ii,ic,:)=mean(X0(J(K(1:N)),:));
            %M1(ii,ic,:)=mean(X0(J,:));
        end
    end
    %Out = fit_ab(D,conref,concom,opts);
    %calfun = @fit_ab3
    cal = str2func(calfun);
    Out = cal(D,conref,concom,opts);
    
    Out.Mresp =M;
    Out.conref = conref;
    Out.concom = concom;
    Out.cons = cons;
    Out.oris = oris;
    Out.lcon = lcon;
    if strcmp(calfun,'fit_ab') || strcmp(calfun,'fit_ab2')
        Out.lams = opts.lams;
    end
    Out.Ncv = opts.Ncv;
    Out.calfun = calfun;
    
end



%%

%-----------------------------------------
function Out = fit_ab(D,ir,is,opts)
% function Out = fit_ab(D,ir,is,opts)
%     D: trialxcellxcontrast
%     ir: reference contrast index
%     is: selected contrast index
%     lams = opts.lams;
%     Ncv = opts.Ncv;
%     O = opts.O;
%     N = opts.N;
    lams = opts.lams;
    Ncv = opts.Ncv;
    O = opts.O;
    N = opts.N;
    T=size(D,1);
    assert(N*O==T,'dimension mismatch');
    
    K=[0 0; 0 1];
    Out0(size(D,2))=struct;
    Out=struct;
    
    parfor ix = 1 : size(D,2)
        D1=D; % just for parfor        
        ev = zeros(Ncv,length(lams));
        as = zeros(2,Ncv,length(lams));
        y = D1(:,ix,is);
        x = D1(:,ix,ir);
        
        for il = 1 :length(lams)
            lam = lams(il);
            for icv = 1 : Ncv
                
                y1 = zeros(O,1);
                y2 = zeros(O,1);
                x1 = zeros(O,2);
                x2 = zeros(O,2);
                for io = 1:O
                    N1 = randperm(N);
                    inxtr = N1(1:round(N*0.5));
                    inxte = N1(round(N*0.5)+1:end);
                    y1(io) = mean(y(inxtr+N*(io-1)));
                    y2(io) = mean(y(inxte+N*(io-1)));
                    x1(io,:) = [mean(x(inxtr+N*(io-1))) 1];
                    x2(io,:) = [mean(x(inxte+N*(io-1))) 1];                 
                end
                


                M =x1'*x1+lam*K;
                a =pinv(M)*x1'*y1;
                ev(icv,il) = 1-sum((x2*a-(y2)).^2)/sum((y2-mean(y2)).^2);
                as(:,icv,il) = a;
            end
        end
        ev(ev(:)<0 | ev(:)>1)=NaN;
        a1=squeeze(as(1,:,:));
        b1=squeeze(as(2,:,:));
        a1(ev(:)<0 | ev(:)>1)=NaN;
        b1(ev(:)<0 | ev(:)>1)=NaN;
        
        Out0(ix).ev = nanmean(ev);
        Out0(ix).as(1,:) = nanmean(a1,1);
        Out0(ix).as(2,:) = nanmean(b1,1);
    end
    for ix = 1 : size(D,2)
        Out.ev(:,ix)=Out0(ix).ev;
        Out.as(:,:,ix)=Out0(ix).as;
    end
end

%%

%-----------------------------------------
function Out = fit_ab2(D,ir,is,opts)
% function Out = fit_ab(D,ir,is,opts)
%     D: trialxcellxcontrast
%     ir: reference contrast index
%     is: selected contrast index
%     lams = opts.lams;
%     Ncv = opts.Ncv;
%     O = opts.O;
%     N = opts.N;
    lams = opts.lams;
    Ncv = opts.Ncv;
    O = opts.O;
    N = opts.N;
    T=size(D,1);
    assert(N*O==T,'dimension mismatch');
    
    K=[0 0; 0 1];
    Out0(size(D,2))=struct;
    Out=struct;
    
    for ix = 1 : size(D,2)
        D1=D; % just for parfor        
        ev = zeros(Ncv,length(lams));
        as = zeros(2,Ncv,length(lams));
        y = D1(:,ix,is);
        x = D1(:,ix,ir);
        
        for il = 1 :length(lams)
            lam = lams(il);
            
            y1 = zeros(O*Ncv,1);
            y2 = zeros(O*Ncv,1);
            x1 = zeros(O*Ncv,2);
            x2 = zeros(O*Ncv,2);
            for icv = 1 : Ncv                
                for io = 1:O
                    N1 = randperm(N);
                    inxtr = N1(1:round(N*0.5));
                    inxte = N1(round(N*0.5)+1:end);
                    y1((icv-1)*O+io) = mean(y(inxtr+N*(io-1)));
                    y2((icv-1)*O+io) = mean(y(inxte+N*(io-1)));
                    x1((icv-1)*O+io,:) = [mean(x(inxtr+N*(io-1))) 1];
                    x2((icv-1)*O+io, :) = [mean(x(inxte+N*(io-1))) 1];                 
                end
            end
                scale = mean(x1(:,1),1);

                x1(:,2)=scale*x1(:,2);
                x2(:,2)=scale*x2(:,2);

                M =x1'*x1+lam*K;
                a =pinv(M)*x1'*y1;
                ev(icv,il) = 1-var(y2-x2*a)/var(y2);%1-sum((y2-x2*a).^2)/sum((y2).^2);
                as(:,icv,il) = a;
            
        end
        Out0(ix).ev = mean(ev);
        Out0(ix).as = squeeze(mean(as,2));
    end
    for ix = 1 : size(D,2)
        Out.ev(:,ix)=Out0(ix).ev;
        Out.as(:,:,ix)=Out0(ix).as;
    end
end



%-----------------------------------------
function Out = fit_ab3(D,ir,is,opts)
% function Out = fit_ab3(D,ir,is,opts)%   
%     J = || r40-r100*w||^2 -t*log(w)
%     d
%     D: trialxcellxcontrast
%     ir: reference contrast index
%     is: selected contrast index
%     lams = opts.lams;
%     Ncv = opts.Ncv;
%     O = opts.O;
%     N = opts.N;

    Ncv = opts.Ncv;
    O = opts.O;
    N = opts.N;
    T=size(D,1);
    assert(N*O==T,'dimension mismatch');
    

    Out=struct;
    NC = size(D,2);
    as = zeros(2,1,NC);
    ev = NaN*ones(1,NC);
    parfor ix = 1 : size(D,2)
        D1=D; % just for parfor        

        y = D1(:,ix,is);
        x = D1(:,ix,ir);       

            
        y1 = zeros(O*Ncv,1);        
        x1 = zeros(O*Ncv,2);
        
        for icv = 1 : Ncv                
            for io = 1:O
                N1 = randperm(N);
                inxtr = N1(1:round(N*0.5));                
                y1((icv-1)*O+io) = mean(y(inxtr+N*(io-1)));                
                x1((icv-1)*O+io,:) = [mean(x(inxtr+N*(io-1))) 1];                
            end
        end
              
        [w,status]=l1_ls_nonneg(x1,y1,0,1e-3,true);
        if strcmp(status, 'Failed')
        else
            w = w.*[1 mean(x1(:,2))]';
            ev(1,ix) = 1-var(y1-x1*w)/var(y1);        
            as(:,1,ix)=w;
        end
    end
    Out.ev=ev;
    Out.as=as;

end



%-----------------------------------------
function Out = fit_ab4(D,ir,is,opts)
% function Out = fit_ab4(D,ir,is,opts)%   
%     J = || r40-r100*w||^2 +lambda|w| -t*log(w)
%     w(1): scale, w(2)+(w3): biase
%     D: trialxcellxcontrast
%     ir: reference contrast index
%     is: selected contrast index
%     lams = opts.lams;
%     Ncv = opts.Ncv;
%     O = opts.O;
%     N = opts.N;


    lams = opts.lams;
    O = opts.O;
    N = opts.N;
    T=size(D,1);
    assert(N*O==T,'dimension mismatch');
    

    Out=struct;
    NC = size(D,2);
    as = zeros(2,1,NC);
    ev = NaN*ones(1,NC);
    parfor ix = 1 : size(D,2)
        D1=D; % just for parfor    
        [x1, y1, x2, y2, sc_bias] = get_data4fit_ab4(D1,ix,ir, is,opts);
        
        ev1 = NaN*ones(1,length(lams));
        ws =  NaN*ones(3,length(lams));
        for il = 1: length(lams)
            lam = lams(il);
            %if ir == is
            %    [w,status]=l1_ls(x1,y1-x1(:,1),lam,1e-3,true);
            %    w(1)=w(1)+1;
            %else
                [w,status]=l1_ls_nonneg(x1,y1,lam,1e-3,true);
            %end
            
            
            
            if strcmp(status, 'Failed')
            else
                
                ev1(il) = 1-var(y1-x1*w)/var(y1);                     
                %ev1(il) = 1-sum((y1-x1*w).^2)/sum(y1.^2);                     
                
                ws(:,il)=w;               

            end
        end
        if ~all(isnan(ev1))
            [~,inx]=max(ev1);
            wopt = ws(:,inx);
            ev(1,ix) = 1-var(y2-x2*wopt)/var(y2);        
            %ev(1,ix) = 1-sum((y2-x2*wopt).^2)/sum(y2.^2);        
            
            wsc = wopt.*[1 sc_bias -sc_bias]';            
            as(:,1,ix) = [wsc(1) mean(wsc(2:3))];
            
            
            
        end
    end
    Out.ev=ev;
    Out.as=as;

end




%-----------------------------------------
function [x1, y1, x2, y2, scale] = get_data4fit_ab4(D,icell,ir, is,opts)
    Ncv = opts.Ncv;
    O = opts.O;
    N = opts.N;
    T=size(D,1);
    assert(N*O==T,'dimension mismatch');
    
            
    y1 = zeros(O*Ncv,1);        
    x1 = zeros(O*Ncv,2);
    y2 = zeros(O*Ncv,1);        
    x2 = zeros(O*Ncv,2);
    
  
    
    if is==ir
        d = D(:,icell,is);
        for icv = 1 : Ncv                
            for io = 1:O
%                 N1 = randperm(N);
%                 N2 = randperm(N);
%                 inx1 = N1(1:round(N*0.5));
%                 inx2 = N1(round(N*0.5)+1:end);
%                 inx3 = N2(1:round(N*0.5));
%                 inx4 = N2(round(N*0.5)+1:end);
                
                
                N1 = randperm(round(N/2));
                N0 = N-length(N1);
                N2 = length(N1)+ randperm(N0);
                inx1 = N1(1:round(length(N1)*0.5));
                inx2 = N1(round(length(N1)*0.5)+1:end);
                inx3 = N2(1:round(length(N2)*0.5));
                inx4 = N2(round(length(N2)*0.5)+1:end);



                y1((icv-1)*O+io) = mean(d(inx1+N*(io-1)));                
                x1((icv-1)*O+io,1) = mean(d(inx2+N*(io-1)));  
                y2((icv-1)*O+io) = mean(d(inx3+N*(io-1)));                
                x2((icv-1)*O+io,1) = mean(d(inx4+N*(io-1))); 

            end
        end
    else
        y = D(:,icell,is);
        x = D(:,icell,ir);  

        for icv = 1 : Ncv                
            for io = 1:O
                N1 = randperm(N);
                inxtr = N1(1:round(N*0.5));                
                y1((icv-1)*O+io) = mean(y(inxtr+N*(io-1)));                
                x1((icv-1)*O+io,1) = mean(x(inxtr+N*(io-1)));  

                inxval = N1(round(N*0.5)+1:end);
                y2((icv-1)*O+io) = mean(y(inxval+N*(io-1)));
                x2((icv-1)*O+io, 1) = mean(x(inxval+N*(io-1)));
            end
        end


    end
    scale =  mean(x1(:,1));
    x1(:,2) = scale;
    x1(:,3) = -scale;
    x2(:,2) = scale;
    x2(:,3) = -scale;

end


%-----------------------------------------
% function Out = fit_ab4(D,ir,is,opts)
% % function Out = fit_ab4(D,ir,is,opts)%   
% %     J = || r40-r100*w||^2 +lambda|w| -t*log(w)
% %     w(1): scale, w(2)+(w3): biase
% %     D: trialxcellxcontrast
% %     ir: reference contrast index
% %     is: selected contrast index
% %     lams = opts.lams;
% %     Ncv = opts.Ncv;
% %     O = opts.O;
% %     N = opts.N;
% 
% 
%     lams = opts.lams;
%     Ncv = opts.Ncv;
%     O = opts.O;
%     N = opts.N;
%     T=size(D,1);
%     assert(N*O==T,'dimension mismatch');
%     
% 
%     Out=struct;
%     NC = size(D,2);
%     as = zeros(2,1,NC);
%     ev = NaN*ones(1,NC);
%     parfor ix = 1 : size(D,2)
%         D1=D; % just for parfor        
%         if is==ir
%         else
%             y = D1(:,ix,is);
%             x = D1(:,ix,ir);       
% 
% 
%             y1 = zeros(O*Ncv,1);        
%             x1 = zeros(O*Ncv,3);
%             y3 = zeros(O*Ncv,1);        
%             x3 = zeros(O*Ncv,3);
% 
%             for icv = 1 : Ncv                
%                 for io = 1:O
%                     N1 = randperm(N);
%                     inxtr = N1(1:round(N*0.5));                
%                     y1((icv-1)*O+io) = mean(y(inxtr+N*(io-1)));                
%                     x1((icv-1)*O+io,1) = mean(x(inxtr+N*(io-1)));  
% 
%                     inxval = N1(round(N*0.5)+1:end);
%                     y3((icv-1)*O+io) = mean(y(inxval+N*(io-1)));
%                     x3((icv-1)*O+io, :) = [mean(x(inxval+N*(io-1))) 1 -1];
%                 end
%             end
% 
%             x1(:,2) = mean(x1(:,1));
%             x1(:,3) = -1*mean(x1(:,1));
%         end
%         
%         ev1 = NaN*ones(1,length(lams));
%         ws =  NaN*ones(3,length(lams));
%         for il = 1: length(lams)
%             lam = lams(il);
%             [w,status]=l1_ls_nonneg(x1,y1,lam,1e-3,true);
%             if strcmp(status, 'Failed')
%             else
%                 ev1(il) = 1-var(y1-x1*w)/var(y1);       
%                 wsc = w.*[1 mean(x1(:,2)) mean(x1(:,3))]';
%                 ws(:,il)=wsc;
%             end
%         end
%         if ~all(isnan(ev1))
%             [~,inx]=max(ev1);
%             wopt = ws(:,inx);
%             ev(1,ix) = 1-var(y3-x3*wopt)/var(y3);        
%             as(:,1,ix) = [wopt(1) mean(wopt(2:3))];
%         end
%     end
%     Out.ev=ev;
%     Out.as=as;
% 
% end
