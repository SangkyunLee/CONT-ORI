% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5_eyethr_xy1_p1'
% iexp_type=5;
% 
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt)

% 
% 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5_eyethr1'
% iexp_type=5;
% 
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt)
% 
% pause(10);
% 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{20,20},...
% {100, 40},{100,20},...
% {40,100},{40,20},...
% {20,100},{20,40},...
% {[100 40],100},{[100 40],40},...
% {[100 40 20],100},{[100 40 20],40},{[100 40 20],20},...
% };
% 
% DATA_thr_str = 'thr5'
% iexp_type=1;
% 
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt)
% 
%% --------------------------------
% ctm = 0.6;
% for iexp_type = [1 2 5]
%     if iexp_type==5,
%         DATA_thr_str = 'thr5_eyethr_xy1_p1';
%     else
%         DATA_thr_str = 'thr5';
%     end
%         
%     get_ORI_conscale2(iexp_type,DATA_thr_str,ctm);
%         
% end
% 
%% ----------------------



% marker=15380912
% cellselopt.CELLSEL_THRs=[0 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% 
% DATA_thr_str = 'thr5'
% iexp_type=2;
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt, false,0.6)


% DATA_thr_str = 'thr5_eyethr_xy1_p1'
% iexp_type=5;
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt, false,0.6)



%-----------------------
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% DATA_thr_str = 'thr5_eyethr_xy1_p1'
% iexp_type=5;
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt, false,0.6,[],true)




%------------------
% 
% marker=1836; %for 1 NOBIAS_ZMEAN
% 
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% for i=[1 ]
% fnpf1='P1'
% fnpf2='P2'
% get_subcell_crsdecoding(i,fnpf1,fnpf2,DATA_thr_str{i},false,true)
% fnpf1='P2'
% fnpf2='P1'
% get_subcell_crsdecoding(i,fnpf1,fnpf2,DATA_thr_str{i},false,true)
% end
% 
% 
% % 
% % %-------------
% % 
% % DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% % for i=[1 2 5 ]   
% % get_ORI_conscale2(i,DATA_thr_str{i},0.6,'fit_ab3')
% % end
% 
% %-----------------
% % % 
% % for iexp_type=[1 2 5]
% % get_decscore(iexp_type,[],0:0.1:1);
% % end
% % %
% % for iexp_type=[5]
% % get_decscore(iexp_type,[],0);
% % end
% 
% %----------------------------
% %for 1 NOBIAS_ZMEAN shuffle
% mk = [1722 627 3456 0 1]
% 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% 
% 
% 
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% for i=[1]
% marker = mk(i)
% % ses = sespart{i}
% get_subcellcrcon_decoding(i, DATA_thr_str{i},...
%    comcont, cellselopt,true,0.6,[],true);
% end
% 
% 
% 
% %------ population activity
% comcont={{'100,H','100,H'},{'100,H','100,L'},{'100,H','40,H'},{'100,H','40,L'},...
%     {'40,H','100,H'},{'40,H','100,L'},{'40,H','40,H'},{'40,H','40,L'},...
%     {'100,L','100,H'},{'100,L','100,L'},{'100,L','40,H'},{'100,L','40,L'},...
%     {'40,L','100,H'},{'40,L','100,L'},{'40,L','40,H'},{'40,L','40,L'}};
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% mk = [1 2 3 4 5]
% 
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% fail = cell(1,5);
% for i=[1]
% marker = mk(i)
% 
% fail0 = get_subcellPA_decoding(i, DATA_thr_str{i},...
%    comcont, cellselopt,false,0.6,[],true);
% 
% fail{i}=fail0;
% end
% 
% 
% 
% %%
% %------ get population activity
% comcont={'100,H','100,L','40,H','40,L'};
% ctm=0.6
% DATA_thr_str ={'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% 
% for iexp_type= [1 2 5]
% 
% get_datawithPA(iexp_type, DATA_thr_str{iexp_type}, comcont, ctm)
% 
% 
% end
% 
% 
% %---------- 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% 
% dpart='P1'
% DATA_thr_str = {'','thr5','','','thr5_eyethr_xy1_p1'}
% for i=[2]
% 
% % ses = sespart{i}
% get_subcellcrcon_subdata_decoding(i, DATA_thr_str{i},dpart,...
%    comcont, cellselopt,false,0.6,[],true);
% end
% 
% 
% 
% %
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% for i=[1 2 5 ]   
% get_ORI_conscale2(i,DATA_thr_str{i},0.6,'fit_ab4')
% end
% 
% 
% 
% 
%   
% get_ORI_conscale_PADEP('fit_ab4')
% 
% 
% 
% 
% %% ctm=0.7
% 
% clear all
% 
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% DATA_thr_str = 'thr5'
% iexp_type=1;
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt, false,0.7,[],true)
% 
% 
% DATA_thr_str = 'thr5'
% iexp_type=2;
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt, false,0.7,[],true)
% 
% 
% DATA_thr_str = 'thr5_eyethr_xy1_p1'
% iexp_type=5;
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt, false,0.7,[],true)
% 
% 
% 
% %%
% 
% iexp_type=1,DATA_thr_str='thr5',ctm=0.6,calfunstr='fit_ab4'
% complist=[[100 100]' [40 40]' [100 40]' [40 100]'];
% 
% get_ORI_conscale4(iexp_type,DATA_thr_str,ctm,calfunstr,complist)

%% to match noise levels of between-contrast decoders
% clear all
% 
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'};
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100,100},{100,100,40},{40,100,100},{40,100,40}, ...
%     {40,40,100},{40,40,40},{100,40,100},{100,40,40}};
% cellselopt.CELLSEL_THRs=[0 1 3 5 10];
% 
% nori=6;
% %getDfun={'get_gamdat2',1000}
% % getDfun={'get_gamdat3',1000}
% % getDfun={'get_gamdat4',1000}
%  getDfun={'get_gamdat5',1000}  
% %  getDfun={'get_gamdat6',1000}     
% 
% seslist{1} = [(1:8) 11 12 15 16];
% seslist{2} = 17:22;
% seslist{3} =[23 25 26 27 29 30 32 33 36 40];
% 
% 
% % console
% sublist = 5%7:9%7:9%4:6%1:3;
% for iexp = [2]    
%     parpoolid = set_env(true,10);    
%     
%     get_subcell_gamnoise_crcon_decoding(iexp, DATA_thr_str{iexp},...
%         comcont, cellselopt, false,0.6,seslist{iexp}(sublist),true,getDfun, parpoolid);
%     if ~isempty(parpoolid),
%         delete(parpoolid)   
%     end
% end


%% re-validate fit_ORIscale 
% 1) because old code the smallest lambda with which explained variance of fitting is maximized.
% 2) In this test, I chose the largest lambda with which explained variance of fitting is maximized.
% % particularly, I change the code in "cal_wori.m"
% opts.lams =  10.^(linspace(-3,1,20));
% into
% opts.lams =  fliplr(10.^(linspace(-3,1,20)));

DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}


parpoolid = set_env(true);

DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
for i=[1 2 5 ]   
get_ORI_conscale2(i,DATA_thr_str{i},0.6,'fit_ab4',parpoolid)
end
for i=[1 2 5 ]   
get_ORI_conscale3(i,DATA_thr_str{i},0.6,'fit_ab4',parpoolid)
end

if ~isempty(parpoolid),
    delete(parpoolid)
    pause(3);
end
