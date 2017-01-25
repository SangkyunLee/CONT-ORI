clear all

D(1) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\BASIC_SUMMARY_AN1-16_ctm0.60.mat');
D(2) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AN\thr5\BASIC_SUMMARY_AN17-22_ctm0.60.mat');
D(3) = load('Z:\data_2photon\matlab_2ndLev\GRP_data\AWAKE_EYE\thr5_eyethr_xy1_p1\BASIC_SUMMARY_AW23-40_ctm0.60.mat');

thr = 0.5
fs = plot_orisc(thr,false);
 
 nses = length(fs);
 a=[]; b=[]; Y=[];
 for i = 1 : nses
     iexp = fs(i).iexp;
     ises = fs(i).ises;
     X = D(iexp).S(ises);
     icell = fs(i).cell;     
     %V =[V mean(mean(X.V(icell,:,1:2),3),2)'];
     %M =[M mean(mean(X.m(icell,:,1:2),3),2)'];
     %Y =[Y mean(abs(X.V(icell,:,1)-X.V(icell,:,2)),2)'];
  
     
     a = [a fs(i).a];
     b = [b fs(i).b];
 end
     
[beta r handles] = plot_linearreg([Y' a'])


%%
 mean2 = @(x)mean(x,2);
 for i = 1 : nses
     iexp = fs(i).iexp;
     ises = fs(i).ises;
     X = D(iexp).S(ises);
     icell = fs(i).cell;     
     
     P = cellfun(mean2,X.Raw,'UniformOutput', false);
     Ps = cellfun(@sepdata,P,repmat({[[0 0.5]' [0.5 1]']},size(P)),'UniformOutput', false);
     
 end