% This script works for eyetracking and movemovement detected by rotaroad.
clear all
close all

% if ispc
%     addpath(genpath('Z:\codes2P'));
% elseif isunix
%     addpath(genpath('/home/slee/data/codes2P'));
% end


% clear all
% close all

if isunix
    mainpath4 ='/media/sdb_WD4T/data_2photon/Sedated'; 
    mainpath3 ='/media/sdb_WD4T/data_2photon/AWAKE';    
    mainpath2='/home/slee/data/data_2photon';
    mainpath1='/media/slee-SS';
else    
    mainpath4 = 'Y:\data_2photon\Sedated'; 
    mainpath3 = 'Y:\data_2photon\AWAKE';    
    mainpath2='Z:\data_2photon';   
    mainpath1='D:\data_2photon';
    
end
exp_type={'AN','AN_0TO150','AWAKE','AN_FULLORI','AWAKE_EYE'};

idtype =5;
if idtype==1 || idtype==4,
    path_list =METADATA_CONT_ORI (exp_type{idtype},mainpath1, mainpath2,mainpath3,mainpath4);
elseif idtype==3
    [path_list, eyeinfo] =METADATA_CONT_ORI (exp_type{idtype},mainpath1, mainpath2,mainpath3);
else
    path_list =METADATA_CONT_ORI (exp_type{idtype},mainpath1, mainpath2,mainpath3);
end
ctm = 0.7;

if idtype==1,
    fndata=sprintf('DISK_Cont-ORI_ctm%0.2f_HPF100_tau0.85_F0-10.0sigma.mat',ctm) ;
else
    fndata=sprintf('DISK_Cont-ORI_ctm%0.2f_HPF100_tau0.85_F0-10.0sigma.mat',ctm) ;   
end


cell_sel_method = 'UNION_CONTRSP';
datatype=['DATA_DISK_' cell_sel_method];
% datatype=['DATA_RIM_' cell_sel_method];
thr_str ='thr5_eyethr_xy1_p1';


fnsave_str =sprintf('%s_ctm%0.2f',datatype,ctm); 
fnsave = fullfile('../GRP_data/', exp_type{idtype} ,thr_str,fnsave_str);




pre_stimtime_ms=500;
motionthr.twin = 1;
motionthr.IMG_motion_thr =2;
motionthr.dist_twin_thr = 1;
motionthr.channel = [3 5]; %DAQ channel for rotaroad


eyepar.lmovthr=2;
eyepar.dxythr=1;
eyepar.puppar.ls = 1;
eyepar.puppar.us = 1;

switch exp_type{idtype}
    case { 'AN_0TO150'}
        selcontrast=[100 40];
        fndats = cell(1,7);
        fndats(:)={fndata};
        fndata = fndats;
        fndata([1 3]) ={['navg2' fndata{1}]};
    case{'AN'}
        selcontrast=[100 40 20];
        %%- when idata=17:22
%         selcontrast=[100 40];
    case {'AWAKE','AWAKE_EYE'}        
        selcontrast=[100 40];        
        
    case{'AN_FULLORI'}
        selcontrast=[100 40];
end
        
% if strcmp(exp_type{idtype},'AN_0TO150')
%     selcontrast=[100 40];
%     fndats = cell(1,7);
%     fndats(:)={fndata};
%     fndata = fndats;
%     fndata([1 3]) ={['navg2' fndata{1}]};
% else
%     selcontrast=[100 40 20];
% end
visual_resp_thr=1.05;
evtsort={'Contrast','Orientation'};


switch exp_type{idtype}
    case {'AN','AN_0TO150','AN_FULLORI','AN2'}
        info.animal_state='AN';
    case {'AWAKE','AWAKE_EYE'}
        info.animal_state='AW';
end
param4cellselect.selcontrast = selcontrast;
param4cellselect.pre_stimtime_ms = pre_stimtime_ms; 
param4cellselect.visual_resp_thr = visual_resp_thr;
switch upper(info.animal_state)
    case {'AW','AWAKE','AWAKE_EYE'}        
        param4cellselect.motionthr = motionthr;  
        param4cellselect.eyethr = eyepar;
    otherwise
        param4cellselect.motionthr=struct([]);

end
param4cellselect.evtsort=evtsort;

%an
%seslist = [23 (25:30) 32 33 (35:37) 40];
% seslist =[11 12 15 16] %1:8;%[(1:8) 11 12 15 16];


seslist= [23 25 26 27 29 30 32 33 36 40]
sepdata = get_scanspersession_eye(exp_type{idtype});

% seslist=[25 26 29 30 32 33 36 40]
% sepdata = get_scanspersession_eye2(exp_type{idtype});

collect_data_BY_contunion(seslist, path_list, info, fndata, fnsave, param4cellselect,sepdata);


% collect_data_BY_contunion([14:length(path_list)], path_list, info, fndata, fnsave, param4cellselect);


