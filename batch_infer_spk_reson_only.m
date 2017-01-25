clc, clear all
addpath('/home/slee/data/data_2photon/matlab_2ndLev/contrast-ori-code2');
addpath('Z:\data_2photon\matlab_2ndLev\contrast-ori-code2');
if isunix
    mainpath='/home/slee/data/data_2photon'
    mainpath1='/media/slee-SS'
    mainpath2='/media/sdb_WD4T/data_2photon/AWAKE'
else
    mainpath='Z:/data_2photon'
    mainpath1='D:/data_2photon'
    mainpath2='Y:\data_2photon\AWAKE'
    addpath(genpath('Z:\codes2P'));
end
    
    

fn_list{1} = fullfile(mainpath,'Thy1_4P3_05042015BN','07102015_thy1_05042015BN','matlab')
fn_list{2} = fullfile(mainpath,'Thy1_4P3_05042015BN','07102015_thy1_05042015BN_2','matlab')
fn_list{3} = fullfile(mainpath,'Thy1_4P3_05042015BR','07092015_thy1_05042015BR_1','matlab')
fn_list{4} = fullfile(mainpath,'Thy1_4P3_05042015BR','07122015_thy1_05042015BR_2','matlab')
fn_list{5} = fullfile(mainpath,'Thy1_4P3_05102015BL','07142015_thy1_0510BL','matlab')
fn_list{6} = fullfile(mainpath,'Thy1_4P3_05102015BL','07142015_thy1_0510BL_2','matlab')
fn_list{7}  = fullfile(mainpath1,'150507BRL','07222015_0507BRL','matlab')

fn_list{8} = fullfile(mainpath2,'0705BR','09182015_0705BR_2','matlab');%reso
fn_list{9} = fullfile(mainpath2,'0705BR','09212015_0705BR_Nat_Grating','matlab'); %reso


p = mfilename('fullpath')
p = fileparts(p) ;

for data_list=5:7%1:4%[(8:9) 2 (4:7) ] %ubuntu2[7:9]
    pathstr = fullfile(fn_list{data_list})
    cd (pathstr)
    if data_list==1 || data_list==3,
        eval('newinfer_spk(10,true,true,100,''Cont-ORI'',2)' );
    else
        eval('newinfer_spk(10,true,true,100,''Cont-ORI'')' );
    end
end
cd(p)