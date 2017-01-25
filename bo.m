

%--------------------------------------


%%
%---------------------------------------
% marker=1818;
% configCluster;
% ClusterInfo.setQueueName('mpi');
% % ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('72:00'); 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5_eyethr_xy0p5_p100'
% iexp_type=5;
% 
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt)

% 
% % %---------------
% marker=1823;
% configCluster;
% % ClusterInfo.setQueueName('mpi');
% ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('72:00'); 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5_eyethr_xy1_p1'
% iexp_type=5;
% 
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%     comcont, cellselopt)

% 




%--------------------
%  %suffle
%   marker=15380916; %for 2
% % 
%  configCluster;
%  % ClusterInfo.setQueueName('mpi');
%  ClusterInfo.setQueueName('parallel');
%  ClusterInfo.setWallTime('96:00'); 
%  cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
%  cellselopt.SEL_METHOD='ncell';
%  comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
%  DATA_thr_str = 'thr5'
%  
%  
%  for iexp_type=1:2
%  
%  get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%      comcont, cellselopt,true,0.6);
%  end



%----------------------------
 %marker=1836; %for 1 NOBIAS_ZMEAN
% mk = [1722 627]
% configCluster;
% %ClusterInfo.setQueueName('mpi');
% ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('96:00'); 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% DATA_thr_str = 'thr5'
% 
% 
% % sespart ={[18 19 17],[20 21]}
% % iexp_type=2
% 
% for iexp_type=2
% marker = mk(i)
% % ses = sespart{i}
% get_subcellcrcon_decoding(iexp_type, DATA_thr_str,...
%    comcont, cellselopt,false,0.6);
% end


%%
%  marker=1836; %for 1 NOBIAS_ZMEAN
% 
% configCluster;
% %ClusterInfo.setQueueName('mpi');
% ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('96:00'); 
% 
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% for i=[1 2 5]
% fnpf1='P1'
% fnpf2='P2'
% get_subcell_crsdecoding(i,fnpf1,fnpf2,DATA_thr_str{i},false,true)
% fnpf1='P2'
% fnpf2='P1'
% get_subcell_crsdecoding(i,fnpf1,fnpf2,DATA_thr_str{i},false,true)
% end


% %----------------------------
% %for 1 NOBIAS_ZMEAN
% mk = [1722 627 3456 0 1] 
% configCluster;
% %ClusterInfo.setQueueName('mpi');
% ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('96:00'); 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% 
% 
% 
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% for i=2%[1 2 5]
% marker = mk(i)
% % ses = sespart{i}
% get_subcellcrcon_decoding(i, DATA_thr_str{i},...
%    comcont, cellselopt,true,0.6,[18],true);
% end


% ------ population activity dependence
% configCluster;
% ClusterInfo.setQueueName('mpi');
% %ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('96:00'); 
% 
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
% for i=[2]
% marker = mk(i)
% 
% fail0 = get_subcellPA_decoding(i, DATA_thr_str{i},...
%    comcont, cellselopt,false,0.6,[],true);
% 
% fail{i}=fail0;
% end



%---------- 

% configCluster;
% ClusterInfo.setQueueName('mpi');
% %ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('96:00'); 
% cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
% cellselopt.SEL_METHOD='ncell';
% comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
% 
% dpart='P2'
% DATA_thr_str = {'','thr5','','','thr5_eyethr_xy1_p1'}
% ses0=[25 26  29 30 32 33 36 40];
% for i=[2 5]
%     i
%     if i==2,
%         ses=[];
%     elseif i==5,
%         ses=ses0; 
%     end
% 
%     get_subcellcrcon_subdata_decoding(i, DATA_thr_str{i},dpart,...
%        comcont, cellselopt,true,0.6,ses,true);
% end



%%
% 
% configCluster;
% % ClusterInfo.setQueueName('mpi');
% ClusterInfo.setQueueName('parallel');
% ClusterInfo.setWallTime('96:00'); 
% 
% 
% 
% DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'}
% for i=[1 2 5 ]   
% get_ORI_conscale3(i,DATA_thr_str{i},0.6,'fit_ab4')
% end

%%

configCluster;
%ClusterInfo.setQueueName('mpi');
ClusterInfo.setQueueName('parallel');
ClusterInfo.setWallTime('96:00'); 


DATA_thr_str = {'thr5','thr5','','','thr5_eyethr_xy1_p1'};
cellselopt.SEL_METHOD='ncell';
comcont={{100,100,100},{100,100,40},{40,100,100},{40,100,40}, ...
    {40,40,100},{40,40,40},{100,40,100},{100,40,40}};
cellselopt.CELLSEL_THRs=[0 1 3 5 10];

nori=6;
% getDfun={'get_gamdat2',1000}
% getDfun={'get_gamdat3',1000}
getDfun={'get_gamdat5',1000}
% getDfun={'get_gamdat6',1000}   
   

seslist{1} = [(1:8) 11 12 15 16];
seslist{2} = 17:22;
seslist{5} =[23 25 26 27 29 30 32 33 36 40];


% console
sublist = 1:2 %1:3;
% sublist = [1] %4:6%1:3 %
for iexp = [2]    
    parpoolid = set_env(true,10);    
    
    get_subcell_gamnoise_crcon_decoding(iexp, DATA_thr_str{iexp},...
        comcont, cellselopt, false,0.6,seslist{iexp}(sublist),true,getDfun, parpoolid);
    if ~isempty(parpoolid),
        delete(parpoolid)   
    end
end

%%
% ses{1}= {(1:8), [11 12 15 16]};
% ses{2}= {17:22,[]};
% ses{5}= {[23 25 26 27 29],[30 32 33 36 40]};
% for iexp = 5
%     for js= 2 
%         parpoolid = set_env(true,nori*length(comcont));
%         ses{iexp}{js}
%         get_subcell_gamnoise_crcon_decoding1(iexp, DATA_thr_str{iexp},...
%             comcont, cellselopt, false,0.6,ses{iexp}{js},true,getDfun, parpoolid);
%         if ~isempty(parpoolid),
%             delete(parpoolid)   
%         end
%     end
% end

