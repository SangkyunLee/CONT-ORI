
function [contrasts, ORI_list, ORI_compindexset, nses, seslist] =get_expinfo(iexp_type)

if iexp_type==1
    contrasts=[100 40 20];
    ORI_list=[-15 0 30 90];
    ORI_compindexset = 1:6;
    nses=22;
    seslist =[(1:8) 11 12 15 16];
elseif iexp_type==2
    
    contrasts=[100 40];
    ORI_list=[-10 0 30 90];
    ORI_compindexset = 1:6;
    nses=22;
    seslist =17:22;

elseif iexp_type==3
    
    contrasts=[100 40];    
    ORI_list=[0 30 35 60 90 120 150];
    ORI_compindexset = [1 2 3 4 7 12 16 19 21];   
    nses=7;
    seslist =1:7;
elseif iexp_type==4
    contrasts=[100 40 20];    
    ORI_list=[-10 0 30 90];
    ORI_compindexset = 1:6;
    nses=17;
    seslist =1:17;
elseif iexp_type==5
    contrasts=[100 40];    
    ORI_list=[-10 0 30 90];
    ORI_compindexset = 1:6;
    nses=40;
    %seslist =[23 (25:30) 32 33 (35:37) 40];  
    %seslist =[25 28 29 30 32 33 (35:37) 40];  
    seslist=[23 25 26 27 29 30 32 33 36 40];
end
