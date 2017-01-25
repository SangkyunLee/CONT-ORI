function Out=cal_wori2(X1, X2, evt,calfun,rejectthr)
% linear fiting between orientation tuning functions from two different
% data sets
    if size(X1,2)~=size(X2,2)
        error('X1 and X2 should have the same number of cells');
    end


    inxsample1 = evt.inxsample1;
    inxsample2 = evt.inxsample2;
    cons = evt.cons;
    oris = evt.oris;
    
    
    lcon = unique(cons);    
    conref = find(lcon==evt.conref);
    concom = find(lcon==evt.concom);
    
    

    opts.Ncv = 1000;
    if  strcmp(calfun,'fit_ab4')
        opts.lams =  10.^(linspace(-10,1,40));
    end

    N1 = min(cellfun(@length,inxsample1));
    N2 = min(cellfun(@length,inxsample2));
    N = min(N1,N2);
    if N < rejectthr
        Out.N=N;
        Out.status ='fail'
        return;
    end
    O = length(unique(oris));
    opts.O = O;
    opts.N = N;
    T = N*O;
    D1 = zeros(T,size(X1,2), length(lcon));
    D2 = zeros(T,size(X2,2), length(lcon));
    
    M1 = zeros(O,length(lcon),size(X1,2));
    M2 = zeros(O,length(lcon),size(X1,2));
    
    for ic = 1 : length(lcon)
        iorg = find(cons==lcon(ic));
        fprintf('CON=%d,ORI=%s\n',lcon(ic),num2str(oris(iorg)));
        for ii = 1: length(iorg)
            J = inxsample1{iorg(ii)};
            K = randperm(length(J));
            D1((ii-1)*N+(1:N),:,ic) = X1(J(K(1:N)),:);
            M1(ii,ic,:)=mean(X1(J(K(1:N)),:));
            
            J = inxsample2{iorg(ii)};
            K = randperm(length(J));
            D2((ii-1)*N+(1:N),:,ic) = X2(J(K(1:N)),:);
            M2(ii,ic,:)=mean(X2(J(K(1:N)),:));
            
           
        end
    end
    cal = str2func(calfun);
    Out = cal(D1,D2,conref,concom,opts);
    
    Out.Mresp1 = M1;
    Out.Mresp2 = M2;
    Out.conref = conref;
    Out.concom = concom;
    Out.cons = cons;
    Out.oris = oris;
    Out.lcon = lcon;
    if strcmp(calfun,'fit_ab4')
        Out.lams = opts.lams;
    end
    Out.Ncv = opts.Ncv;
    Out.calfun = calfun;
    Out.status ='success'
    
end



%%


%-----------------------------------------
function Out = fit_ab4(D1,D2,ir,is,opts)
% function Out = fit_ab4(D1,D2,ir,is,opts)%   
%     J = || r40-r100*w||^2 +lambda|w| -t*log(w)
%     w(1): scale, w(2)+(w3): biase
%     D1, D2: trialxcellxcontrast
%     ir: reference contrast index
%     is: selected contrast index
%     lams = opts.lams;
%     Ncv = opts.Ncv;
%     O = opts.O;
%     N = opts.N;


    lams = opts.lams;
    O = opts.O;
    N = opts.N;
    T1 = size(D1,1);
    T2 = size(D2,1);
    assert(N*O==T1 && N*O==T2,'dimension mismatch');
    

    Out=struct;
    NC = size(D1,2);
    as = zeros(2,1,NC);
    ev = NaN*ones(1,NC);
    parfor ix = 1 : size(D1,2)
           
        [x1, y1, x2, y2] = get_data4fit_ab4(D1,D2,ix,ir, is,opts);
        
        ev1 = NaN*ones(1,length(lams));
        ws =  NaN*ones(3,length(lams));
        for il = 1: length(lams)
            lam = lams(il);
            [w,status]=l1_ls_nonneg(x1,y1,lam,1e-3,true);
            if strcmp(status, 'Failed')
            else
                ev1(il) = 1-var(y1-x1*w)/var(y1);       
                wsc = w.*[1 mean(x1(:,2)) mean(x1(:,3))]';
                ws(:,il)=wsc;
            end
        end
        if ~all(isnan(ev1))
            [~,inx]=max(ev1);
            wopt = ws(:,inx);
            ev(1,ix) = 1-var(y2-x2*wopt)/var(y2);        
            as(:,1,ix) = [wopt(1) mean(wopt(2:3))];
        end
    end
    Out.ev=ev;
    Out.as=as;

end




%-----------------------------------------
function [x1, y1, x2, y2] = get_data4fit_ab4(D1,D2,icell,ir, is,opts)
    Ncv = opts.Ncv;
    O = opts.O;
    N = opts.N;

    
            
    y1 = ones(O*Ncv,1);        
    x1 = ones(O*Ncv,3);
    y2 = ones(O*Ncv,1);        
    x2 = ones(O*Ncv,3);
    
  
    
    if is==ir
        d1 = D1(:,icell,ir);
        d2 = D2(:,icell,is);
        for icv = 1 : Ncv                
            for io = 1:O
                N1 = randperm(N);
                N2 = randperm(N);
                inx1 = N1(1:round(N*0.5));
                inx2 = N1(round(N*0.5)+1:end);
                inx3 = N2(1:round(N*0.5));
                inx4 = N2(round(N*0.5)+1:end);



                y1((icv-1)*O+io) = mean(d1(inx1+N*(io-1)));                
                x1((icv-1)*O+io,1) = mean(d1(inx2+N*(io-1)));  
                y2((icv-1)*O+io) = mean(d2(inx3+N*(io-1)));                
                x2((icv-1)*O+io,1) = mean(d2(inx4+N*(io-1))); 

            end
        end
    else
        
        x1r = D1(:,icell,ir);  
        y1s = D1(:,icell,is);
        x2r = D2(:,icell,ir);  
        y2s = D2(:,icell,is);

        for icv = 1 : Ncv                
            for io = 1:O
                N1 = randperm(N);
                inxtr = N1(1:round(N*0.5));                
                y1((icv-1)*O+io) = mean(y1s(inxtr+N*(io-1)));                
                x1((icv-1)*O+io,1) = mean(x1r(inxtr+N*(io-1)));  

                inxval = N1(round(N*0.5)+1:end);
                y2((icv-1)*O+io) = mean(y2s(inxval+N*(io-1)));
                x2((icv-1)*O+io, 1) = mean(x2r(inxval+N*(io-1)));
            end
        end


    end
    x1(:,2) = mean(x1(:,1));
    x1(:,3) = -1*mean(x1(:,1));

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
