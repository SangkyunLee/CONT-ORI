% clc, clear all

if isunix
addpath('/home/slee/data/data_2photon/matlab_2ndLev/contrast-ori-code2');
    mainpath='/media/sdb_WD4T/data_2photon/AWAKE';
else
    addpath('Z:\data_2photon\matlab_2ndLev\contrast-ori-code2');
    mainpath='Y:\data_2photon\AWAKE';
end
    
fn_list{1} = fullfile(mainpath,'0705BL','09142015_0705BL','matlab');%spiral
fn_list{2} = fullfile(mainpath,'0705BL','09142015_0705BL_2','matlab');%sprial
fn_list{4} = fullfile(mainpath,'0705BN','09152015_0705BN','matlab');%spiral
fn_list{6} = fullfile(mainpath,'0705BL','10072015_0705BL_Nat_Grating','matlab');%spiral
fn_list{7} = fullfile(mainpath,'0705BR','10042015_0705BR_Nat_Grating','matlab');%spiral
fn_list{8} = fullfile(mainpath,'0705BR','10122015_0705BR_Nat_Grating','matlab');%spiral
fn_list{9} = fullfile(mainpath,'0705BRL','10032015_0705BRL_Nat_Grating','matlab');%spiral
fn_list{10} = fullfile(mainpath,'0705BRL','10112015_0705BRL_Nat_Grating','matlab');%spiral

%----------------------
fn_list{11} = fullfile(mainpath,'0822BRLF','11082015_0822BRLF','matlab');
fn_list{12} = fullfile(mainpath,'0822BRLF','11082015_0822BRLF_2','matlab');
fn_list{13} = fullfile(mainpath,'0828BLF','11082015_0828BLF','matlab');
fn_list{14} = fullfile(mainpath,'0828BLF','11082015_0828BLF_2','matlab');

fn_list{15} = fullfile(mainpath,'thy1_0822BL','11042015_thy1_0822BL_2','matlab');   
fn_list{16} = fullfile(mainpath,'thy1_0822BN','11052015_thy1_0822BN','matlab');
fn_list{17} = fullfile(mainpath,'thy1_0822BN','11052015_thy1_0822BN_2','matlab');
%-------------------------------------

fn_list{21} = fullfile(mainpath,'thy1_0129B2','04192016_thy1_0129B2_3','matlab');
fn_list{22} = fullfile(mainpath,'thy1_0129B1','04192016_thy1_0129B1_3','matlab');

%------------------------------------------
fn_list{23} = fullfile(mainpath,'thy1_0129B1','04272016_thy1_0129B1_4','matlab');
fn_list{24} = fullfile(mainpath,'thy1_0129B1','04272016_thy1_0129B1_5','matlab'); % two sessions, cell not identified yet
fn_list{25} = fullfile(mainpath,'thy1_0129B1','05122016_thy1_0129B1_9','matlab'); % two sessions
fn_list{26} = fullfile(mainpath,'thy1_0129B2','04282016_thy1_0129B2_4','matlab');
fn_list{27} = fullfile(mainpath,'thy1_0129B2','04282016_thy1_0129B2_5','matlab'); % two sessions, cell not identifed
fn_list{28} = fullfile(mainpath,'thy1_0129B2','05062016_thy1_0129B2_6','matlab'); % two sessions

fn_list{29} = fullfile(mainpath,'thy1_0129B3','05102016_thy1_0129B3_1','matlab');% two sessions
fn_list{30} = fullfile(mainpath,'thy1_0129B3','05112016_thy1_0129B3_2','matlab'); % two sessions, cell not identifed
fn_list{31} = fullfile(mainpath,'thy1_0129B3','05182016_thy1_0129B3_3','matlab'); % one session(0:30:330)

fn_list{32} = fullfile(mainpath,'022216BL','05132016_022216BL_2','matlab');% two sessions
fn_list{33} = fullfile(mainpath,'022216BL','05152016_022216BL_3','matlab');% two sessions
fn_list{34} = fullfile(mainpath,'022216BL','05182016_022216BL_4','matlab');% two sessions (0:30:330)

fn_list{35} = fullfile(mainpath,'022216BN','05132016_022216BN_2','matlab');% two sessions
fn_list{36} = fullfile(mainpath,'022216BN','05152016_022216BN_3','matlab');% two sessions 
%---------
fn_list{37} = fullfile(mainpath,'022016BN','05162016_022016BN_1','matlab');% two sessions
fn_list{38} = fullfile(mainpath,'022016BN','05182016_022016BN_2','matlab');% two sessions 
fn_list{39} = fullfile(mainpath,'022016BRF','05172016_022016BRF','matlab');% two sessions
fn_list{40} = fullfile(mainpath,'thy1_032016BR','05172016_Thy1_032016BR_2','matlab');% two sessions

p = mfilename('fullpath');
p = fileparts(p);
not_performed_set=[];
for data_list= [23 25 26 27 29 30 32 33 36 40]%32:34%35:36%29:30 %[26 27 28]%[1 2 4 (6:10)]%11:17%[(11:17) 1 2 4 (6:10)]
    pathstr = fullfile(fn_list{data_list});
    fprintf('%s\n',pathstr);
    cd (pathstr);
    pause(1);
%     open('infer_spk.m')
%     try
%         if data_list<6
%             eval('newinfer_spk(10,true,true,100,''Cont'')' );
%         else
            %eval('newinfer_spk(3,true,true,100,''Cont-ORI'')' ); %RIM-spatial     
            eval('newinfer_spk(10,false,true,100,''Cont-ORI'')' ); %DISK-spatial
            
%             eval('newinfer_spk(0,true,false,100,''Cont-ORI'')' ); %avg-           
%         end
%     catch err
%         not_performed_set=[not_performed_set data_list];
%         continue;
%     end
end
cd(p)