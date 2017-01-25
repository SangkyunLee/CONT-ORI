% clc, clear all

if isunix
    addpath('/home/slee/data/data_2photon/matlab_2ndLev/contrast-ori-code2');
    mainpath='/home/slee/data/data_2photon';
    mainpath1 = '/media/sdb_WD4T/data_2photon/Sedated';
else
    %addpath('Z:\data_2photon\matlab_2ndLev\contrast-ori-code2');
    mainpath='Z:\data_2photon\';
end

fn_list{1} = fullfile(mainpath,'110315BN','02252016_110315BN','matlab');%spiral
fn_list{2} = fullfile(mainpath,'110315BN','02252016_110315BN_2','matlab');%sprial


fn_list{3} = fullfile(mainpath1,'thy1_012916B1','04172016_thy1_0129','matlab');%spiral512, High lum - bad light shielding
fn_list{4} = fullfile(mainpath1,'thy1_012916B1','04172016_thy1_0129_2','matlab');%sprial256x256, high lum - bad light shielding

fn_list{5} = fullfile(mainpath1,'thy1_012916B2','04192016_thy1_0129B2','matlab');%sprial-lefteyeblock
fn_list{6} = fullfile(mainpath1,'thy1_012916B2','04192016_thy1_0129B2_2','matlab');%sprial -lefteyeblock

fn_list{7} = fullfile(mainpath1,'thy1_012916B2','05072016_thy1_0129B2_7','matlab');%sprial -lefteyeblock, two-sessions
fn_list{8} = fullfile(mainpath1,'thy1_012916B1','05062016_thy1_0129B1_8','matlab');%sprial -lefteyeblock, two-sessions
fn_list{9} = fullfile(mainpath1,'thy1_030416BN','05052016_thy1_030416BN','matlab');%sprial -lefteyeblock, two-sessions
fn_list{10} = fullfile(mainpath1,'thy1_032016BR','05102016_thy1_032016BR','matlab');%sprial -lefteyeblock, two-sessions
fn_list{11} = fullfile(mainpath1,'022216BL','05112016_022216BL_1','matlab');%sprial -lefteyeblock, two-sessions
fn_list{12} = fullfile(mainpath1,'022216BN','05122016_022216BN_1','matlab');%sprial -lefteyeblock, two-sessions



% full ori experiment
fn_list{21} = fullfile(mainpath1,'thy1_012916B1','04292016_thy1_0129B1_6','matlab');%sprial-lefteyeblock
fn_list{22} = fullfile(mainpath1,'thy1_012916B1','05032016_thy1_0129B1_7','matlab');%sprial-lefteyeblock-mainBCM
fn_list{23} = fullfile(mainpath1,'thy1_012916B2','05112016_thy1_0129B2_8','matlab');%sprial -lefteyeblock

p = mfilename('fullpath');
p = fileparts(p);
not_performed_set=[];
for data_list=[1 2 (5:12)]%[(7:12) 23]
    pathstr = fullfile(fn_list{data_list});
    cd (pathstr);
%     open('infer_spk.m')
%     try
        
        %eval('newinfer_spk(3,true,true,100,''Cont-ORI'')' ); 
        %eval('newinfer_spk(3,false,false,100,''Cont-ORI'')' ); 
        
        eval('newinfer_spk(10,false,true,100,''Cont-ORI'')' ); 
%         if data_list<6
%             eval('infer_spk(5,0,100,''Con'')' );
%         else
%             eval('infer_spk(5,0,100,''Cont-ORI'')' );
%         end
%     catch err
%         not_performed_set=[not_performed_set data_list];
%         continue;
%     end
end
cd(p)