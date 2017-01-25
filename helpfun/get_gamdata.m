
% contcond = comcont{icont}
% contcond{3}='tr'
% ORIsel = ORI_condset(:,icomp);
function [dtr_gam,dte_gam,evttr,evtte,inx_train,inx_test] = get_gamdata(D,EVT, dtype,contcond,ORIsel, getDfun) 
% function [dtr_gam,dte_gam,evttr,evtte,inx_train,inx_test] = get_gamdata(D,EVT, dtype,contcond,ORIsel, getDfun) 
% generate gamma_distribution data from the given data

%---------- select train samples ----------------
    sevts{1} = contcond{1};
    sevts{2} = ORIsel;
    inxsample = TP.select_subdata(EVT,sevts);
    % similar number of samples between two conditions,not exceed 10% more samples               
    inx2 = TP.select_samples_evtratio(D.events_ORI(inxsample),sevts{2}, 1.1);
    inx_train = inxsample(inx2);

    dtr = D.(dtype)(inx_train,:)';
    evttr = D.events_ORI(inx_train);

    

    if all(contcond{1}==contcond{2}) 
        inx_test=[];
        dte=[];
        evtte=[];        
        
    else
        sevts{1} = contcond{2};
        sevts{2} = ORIsel;
        inxsample = TP.select_subdata(EVT,sevts);                        
        inx2 = TP.select_samples_evtratio(D.events_ORI(inxsample),sevts{2}, 1.1);
        inx_test = inxsample(inx2);

        dte = D.(dtype)(inx_test,:)';         
        evtte = D.events_ORI(inx_test);         
    end
    
    if all(contcond{1}==contcond{3})
        dex = dtr;
        evtex = evttr;
    elseif all(contcond{2}==contcond{3})
        dex = dte;
        evtex = evtte;
    else
        sevts{1} = contcond{3};
        sevts{2} = ORIsel;
        inxsample = TP.select_subdata(EVT,sevts);                        
        inx3 = TP.select_samples_evtratio(D.events_ORI(inxsample),sevts{2}, 1.1);
        inx_ex = inxsample(inx3);

        dex = D.(dtype)(inx_ex,:)';         
        evtex = D.events_ORI(inx_ex); 
        
    end
    
    % change response range to percentage
    if mean(dtr(:))<1
        dtr = dtr * 100;
        dte = dte * 100;
        dex = dex * 100;
        bscale = true;
    else
        bscale =false;
    end
    
    
    ori = sevts{2}(:)';
    % mainting the same mean but changing the variance with data 'dex'
    %dtr_gam = get_gamdat1(dtr, dex, evttr, evtex,ori);
    
    get_gamd = str2func(getDfun{1});
    switch getDfun{1}
        case 'get_gamdat1'
            dtr_gam = get_gamd(dtr, dex, evttr, evtex,ori);
        case {'get_gamdat2', 'get_gamdat3', 'get_gamdat4', 'get_gamdat5', 'get_gamdat6'}
            N = getDfun{2};
            [dtr_gam, evttr] = get_gamd(dtr, dex, evttr, evtex,ori,N);
            inx_train=[];
        otherwise
    end
    

    if isempty(inx_test)        
        dte_gam = [];        
    else     
        %dte_gam = get_gamdat1(dte, dex, evtte, evtex,ori);        
        
        switch getDfun{1}
            case 'get_gamdat1'
                dte_gam = get_gamd(dte, dex, evtte, evtex,ori);
            case {'get_gamdat2', 'get_gamdat3', 'get_gamdat4', 'get_gamdat5', 'get_gamdat6'}
                [dte_gam, evtte] = get_gamd(dte, dex, evtte, evtex,ori,N);
                inx_test=[];
            otherwise
        end
        
        
        % old code when contcond{3} ='tr' or 'te'
        % mode = contcond{3};
        % [dtr_gam, dte_gam] = get_gamdat3(dtr, dte, evttr, evtte,ori, mode);
        
    end
    
    
    if bscale
        dtr_gam = dtr_gam / 100;
        dte_gam = dte_gam / 100;
    end
end




function d1_gam = get_gamdat1(d1, d2, evt1, evt2,ori)
% function d1_gam = get_gamdat1(d1, d2, evt1, evt2,ori)
% estimate a new gamma distribution that match with
% the mean of d1 and the variance of d2
    d1_gam = NaN*ones(size(d1));


    
    for io = 1 : length(ori)
        
        
        inx1 = evt1==ori(io);        
        inx2 = evt2==ori(io);
        d1_sub = d1(:,inx1);
        d2_sub = d2(:,inx2);
        

        nc = size(d1_sub,1);
        x1 = NaN*ones(size(d1_sub));

        for ic = 1 : nc
            d1x = d1_sub(ic,:);
            ah1 = gamfit(d1x);
            d2x = d2_sub(ic,:);
            ah2 = gamfit(d2x);
            
            ah = NaN*ones(1,2);            
            ah(1) = (ah1(1)*ah1(2))^2/(ah2(1)*ah2(2)^2);
            ah(2) = ah1(1)*ah1(2)/ah(1);
            x1(ic,:) = gamrnd(ah(1),ah(2),size(d1x));    
        
        end
        d1_gam(:,inx1)=x1;
      
    end 
 
end


function [d1_gam, evtnew] = get_gamdat2(d1, d2, evt1, evt2,ori, N)
% function d1_gam = get_gamdat3(d1, d2, evt1, evt2,ori,N)
% estimate a new gamma distribution that match with
% the mean of d1 and the variance of d2 but generate m*N samples
% m is the number of conditions
    d1_gam = NaN*ones(size(d1,1),2*N);
    evtnew = NaN*ones(1,2*N);


    
    for io = 1 : length(ori)
        
        
        inx1 = evt1==ori(io);        
        inx2 = evt2==ori(io);
        d1_sub = d1(:,inx1);
        d2_sub = d2(:,inx2);
        

        nc = size(d1_sub,1);
        x1 = NaN*ones(size(d1_sub,1),N);

        for ic = 1 : nc
            d1x = d1_sub(ic,:);
            ah1 = gamfit(d1x);          
            d2x = d2_sub(ic,:);
            ah2 = gamfit(d2x);
            
            ah = NaN*ones(1,2);           
            ah(1) = (ah1(1)*ah1(2))^2/(ah2(1)*ah2(2)^2);
            ah(2) = ah1(1)*ah1(2)/ah(1);
            x1(ic,:) = gamrnd(ah(1),ah(2),1,N);    
        
        end
        d1_gam(:,(io-1)*N +(1:N)) = x1;
        evtnew((io-1)*N +(1:N)) = ori(io);
      
    end 
 
end

% in contrast get_gamdata1 and get_gamdat2, get_gamdat3 does not gamfit to raw data
% but using mean and variance of the raw data to construct a new gamma
% distribution. Therefore, the original mean and variances are preserved 
function [d1_gam, evtnew] = get_gamdat3(d1, d2, evt1, evt2,ori, N)
% function d1_gam = get_gamdat3(d1, d2, evt1, evt2,ori,N)
% estimate a new gamma distribution that match with
% the mean of d1 and the variance of d2 but generate m*N samples
% m is the number of conditions
    d1_gam = NaN*ones(size(d1,1),2*N);
    evtnew = NaN*ones(1,2*N);


    
    for io = 1 : length(ori)
        
        
        inx1 = evt1==ori(io);        
        inx2 = evt2==ori(io);
        d1_sub = d1(:,inx1);
        d2_sub = d2(:,inx2);
        

        nc = size(d1_sub,1);
        x1 = NaN*ones(size(d1_sub,1),N);

        for ic = 1 : nc
            d1x = d1_sub(ic,:);
            [a1, b1] = est_gampar(d1x);            
            d2x = d2_sub(ic,:);
            [a2, b2] = est_gampar(d2x);
            
            ah = NaN*ones(1,2);            
            ah(1) = (a1*b1)^2/(a2*b2^2);
            ah(2) = a1*b1/ah(1);
            x1(ic,:) = gamrnd(ah(1),ah(2),1,N);    
        
        end
        d1_gam(:,(io-1)*N +(1:N)) = x1;
        evtnew((io-1)*N +(1:N)) = ori(io);
      
    end 
 
end

function [a, b] = est_gampar(x)
% estimate gamma parameters from mean and variance
    M = mean(x); 
    V = var(x,0);
    a = M*M/V;
    b = V/M;
end





% in contrast get_gamdata1 and get_gamdat2, get_gamdat3 does not gamfit to raw data
% but using mean and variance of the raw data to construct a new gamma
% distribution. Therefore, the original mean and variances are preserved 
function [d1_gam, evtnew] = get_gamdat4(d1, d2, evt1, evt2,ori, N)
% function d1_gam = get_gamdat4(d1, d2, evt1, evt2,ori,N)
% estimate a new gamma distribution that match with
% the mean of d1 and the variance of entad2 but generate m*N samples
% m is the number of conditions
    d1_gam = NaN*ones(size(d1,1),2*N);
    evtnew = NaN*ones(1,2*N);


    
    for io = 1 : length(ori)
        
        
        inx1 = evt1==ori(io);        
        inx2 = evt2==ori(io);
        d1_sub = d1(:,inx1);
        d2_sub = d2(:,inx2);
        

        nc = size(d1_sub,1);
        x1 = NaN*ones(size(d1_sub,1),N);
        
        a2= zeros(nc,1); b2 = zeros(nc,1);
        for ic = 1 : nc
            d2x = d2_sub(ic,:);
            [a2i, b2i] = est_gampar(d2x);
            a2(ic) = a2i;
            b2(ic) = b2i;
        end
        a2 = mean(a2); b2 = mean(b2);
        
        
        for ic = 1 : nc
            d1x = d1_sub(ic,:);
            [a1, b1] = est_gampar(d1x);            
            
            
            ah = NaN*ones(1,2);            
            ah(1) = (a1*b1)^2/(a2*b2^2);
            ah(2) = a1*b1/ah(1);
            x1(ic,:) = gamrnd(ah(1),ah(2),1,N);    
        
        end
        d1_gam(:,(io-1)*N +(1:N)) = x1;
        evtnew((io-1)*N +(1:N)) = ori(io);
      
    end 
 
end


% in contrast get_gamdata1 and get_gamdat2, get_gamdat5 does not gamfit to raw data
% but using mean and variance of the raw data to construct a new gamma
% distribution. Therefore, the original mean and variances are preserved 
function [d1_gam, evtnew] = get_gamdat5(d1, d2, evt1, evt2,ori, N)
% function d1_gam = get_gamdat5(d1, d2, evt1, evt2,ori,N)
% estimate a new gamma distribution that match with
% the mean of d1 and the snr of d2 but generate m*N samples
% m is the number of conditions
    d1_gam = NaN*ones(size(d1,1),2*N);
    evtnew = NaN*ones(1,2*N);


    
    for io = 1 : length(ori)
        
        
        inx1 = evt1==ori(io);        
        inx2 = evt2==ori(io);
        d1_sub = d1(:,inx1);
        d2_sub = d2(:,inx2);
        

        nc = size(d1_sub,1);
        x1 = NaN*ones(size(d1_sub,1),N);

        for ic = 1 : nc
            d1x = d1_sub(ic,:);
            [a1, b1] = est_gampar(d1x);    
            
            
            
            d2x = d2_sub(ic,:);
            [a2, b2] = est_gampar(d2x);
            snr =sqrt(a2);% a2*b2/std(a2*b2^2);
            
            ah = NaN*ones(1,2);            
            ah(1) = snr^2;
            ah(2) = a1*b1/ah(1);
            x1(ic,:) = gamrnd(ah(1),ah(2),1,N);    
        
        end
        d1_gam(:,(io-1)*N +(1:N)) = x1;
        evtnew((io-1)*N +(1:N)) = ori(io);
      
    end 
 
end


% in contrast get_gamdata1 and get_gamdat2, get_gamdat6 does not gamfit to raw data
% but using mean and variance of the raw data to construct a new gamma
% distribution. Therefore, the original mean and variances are preserved 
function [d1_gam, evtnew] = get_gamdat6(d1, d2, evt1, evt2,ori, N)
% function d1_gam = get_gamdat6(d1, d2, evt1, evt2,ori,N)
% estimate a new gamma distribution that match with
% the mean of d1 and the snr of d2 but generate m*N samples
% m is the number of conditions
% identical snr across cells
    d1_gam = NaN*ones(size(d1,1),2*N);
    evtnew = NaN*ones(1,2*N);


    
    for io = 1 : length(ori)
        
        
        inx1 = evt1==ori(io);        
        inx2 = evt2==ori(io);
        d1_sub = d1(:,inx1);
        d2_sub = d2(:,inx2);
        

        nc = size(d1_sub,1);
        x1 = NaN*ones(size(d1_sub,1),N);
        
        a2= zeros(nc,1); b2 = zeros(nc,1);
        for ic = 1 : nc
            d2x = d2_sub(ic,:);
            [a2i, b2i] = est_gampar(d2x);
            a2(ic) = a2i;
            b2(ic) = b2i;
        end
        a2 = mean(a2); b2 = mean(b2);
        
        snr = sqrt(a2);%a2*b2/std(a2*b2^2);
        
        
        for ic = 1 : nc
            d1x = d1_sub(ic,:);
            [a1, b1] = est_gampar(d1x);            
            
            
            ah = NaN*ones(1,2);            
            ah(1) = snr^2;
            ah(2) = a1*b1/ah(1);
            x1(ic,:) = gamrnd(ah(1),ah(2),1,N);    
        
        end
        d1_gam(:,(io-1)*N +(1:N)) = x1;
        evtnew((io-1)*N +(1:N)) = ori(io);
      
    end 
 
end

%%
% 
% % contcond = comcont{icont}
% % contcond{3}='tr'
% % ORIsel = ORI_condset(:,icomp);
% function [dtr_gam,dte_gam,evttr,evtte,inx_train,inx_test] = get_gamdata(D,EVT, dtype,contcond,ORIsel) 
% % function [dtr_gam,dte_gam,evttr,evtte,inx_train,inx_test] = get_gamdata(D,EVT, dtype,contcond,ORIsel) 
% % generate gamma_distribution data from the given data
% 
% %---------- select train samples ----------------
%     sevts{1} = contcond{1};
%     sevts{2} = ORIsel;
%     inxsample = TP.select_subdata(EVT,sevts);
%     % similar number of samples between two conditions,not exceed 10% more samples               
%     inx2 = TP.select_samples_evtratio(D.events_ORI(inxsample),sevts{2}, 1.1);
%     inx_train = inxsample(inx2);
% 
%     dtr = D.(dtype)(inx_train,:)';
%     evttr = D.events_ORI(inx_train);
% 
%     
% 
%     if all(contcond{1}==contcond{2})        
%         inx_test=[];
%         dte=[];
%         evtte=[];
%     else
%         sevts{1} = contcond{2};
%         sevts{2} = ORIsel;
%         inxsample = TP.select_subdata(EVT,sevts);                        
%         inx2 = TP.select_samples_evtratio(D.events_ORI(inxsample),sevts{2}, 1.1);
%         inx_test = inxsample(inx2);
% 
%         dte = D.(dtype)(inx_test,:)';         
%         evtte = D.events_ORI(inx_test);        
%     end
%     
%     % change response range to percentage
%     if mean(dtr(:))<1
%         dtr = dtr * 100;
%         dte = dte * 100;
%     end
%     
%     
%     ori = sevts{2}(:)';
%     if isempty(inx_test)
%         dtr_gam = get_gamdat1(dtr,evttr,ori); 
%         dte_gam = [];
%         
%         
%     else
%         mode = contcond{3};
%         [dtr_gam, dte_gam] = get_gamdat2(dtr, dte, evttr, evtte,ori, mode);
%         
%     end
%     
% end
% 
% 
% function dtr_gam = get_gamdat1(dtr,evttr,ori)
%     dtr_gam = NaN*ones(size(dtr));
%     for io = 1 : length(ori)
%         inx1 = evttr==ori(io);
%         d1 = dtr(:,inx1);
% 
%         nc = size(d1,1);
%         x = NaN*ones(size(d1));
%         for ic = 1 : nc
%             dx = d1(ic,:);
%             ahat = gamfit(dx);
%             x(ic,:) = gamrnd(ahat(1),ahat(2),size(dx));
%         end
%         dtr_gam(:,inx1)=x;
%     end   
% end
% 
% 
% 
% function [dtr_gam, dte_gam] = get_gamdat2(dtr, dte, evttr, evtte,ori, mode)
% % function [dtr_gam, dte_gam] = get_gamdat2(dtr, dte, evttr, evtte,ori, mode)
% % when mode='tr', mean, var of dtr_gam = mean(dtr), var(dtr), dte_gam =mean(dte), var(dtr)
% % when mode='te', mean, var of dtr_gam = mean(dtr), var(dte), dte_gam =mean(dte), var(dte)
% 
% 
%     dtr_gam = NaN*ones(size(dtr));
%     dte_gam = NaN*ones(size(dte));
%     
%     if strcmp(mode,'tr')
%         d1_gam = NaN*ones(size(dtr));
%         d2_gam = NaN*ones(size(dte));
%     elseif strcmp(mode,'te')
%         d1_gam = NaN*ones(size(dte));
%         d2_gam = NaN*ones(size(dtr));
%     end
%     
%     for io = 1 : length(ori)
%         
%         if strcmp(mode,'tr')
%             inx1 = evttr==ori(io);        
%             inx2 = evtte==ori(io);
%             d1 = dtr(:,inx1);
%             d2 = dte(:,inx2);
%         elseif strcmp(mode,'te')
%             inx1 = evtte==ori(io);        
%             inx2 = evttr==ori(io);
%             d1 = dte(:,inx1);
%             d2 = dtr(:,inx2);            
%         end
% 
%         nc = size(d1,1);
%         x1 = NaN*ones(size(d1));
%         x2 = NaN*ones(size(d2));
%         for ic = 1 : nc
%             d1x = d1(ic,:);
%             ah1 = gamfit(d1x);
%             d2x = d2(ic,:);
%             ah2 = gamfit(d2x);
%             
%             x1(ic,:) = gamrnd(ah1(1),ah1(2),size(d1x));
%             % estimate a new gamma distribution that match with
%             % the mean of d2 and the variance of d1
%             
%             ah = NaN*ones(1,2);            
%             ah(1) = (ah2(1)*ah2(2))^2/(ah1(1)*ah1(2)^2);
%             ah(2) = ah2(1)*ah2(2)/ah(1);
%             x2(ic,:) = gamrnd(ah(1),ah(2),size(d2x));    
%         
%         end
%         d1_gam(:,inx1)=x1;
%         d2_gam(:,inx2)=x2;        
%     end 
%     if strcmp(mode,'tr')
%         dtr_gam = d1_gam;
%         dte_gam = d2_gam;
%     elseif strcmp(mode,'te')
%         dtr_gam = d2_gam;
%         dte_gam = d1_gam;
%     end
% end
% 
% % x = [0.1:0.1:20];
% % Y = gampdf(x,ahat(1),ahat(2));            
% % ahat = expfit(a);            
% % Y2 = exppdf(x,ahat);
% % P = Y/sum(Y);
% % P2 = Y2/sum(Y2);
% % figure; plot(x,P);
% % hold on; 
% % plot(x,P2,'k');
% % m=hist(a,x);
% % plot(x,m/sum(m))
% 
% % figure; hold on;
% % m=hist(a,x);
% % plot(x,m);
% % m=hist(b,x);
% % plot(x,m,'r');

    