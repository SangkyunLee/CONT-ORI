% new data analysis pipeline

% specify the scan lists to be used 
% actually, all valid scans (i.e., successful eyetracked scans) belong to a
% single session
get_scanspersession_eye.m

% specify the scan lists to be used 
% actually, all valid scans (i.e., successful eyetracked scans) belong 
% to two separate session with 20-30 mins break
% This is commensurate with get_scanspersession.m.
% but get_scanspersession.m function arbitrary separate two sessions 
% regardless of physical time separation.
% but get_scanspersession_eye2.m function allow only physically time
% seperated exeperiments such as AN17-22
get_scanspersession_eye2.m

% loading scan data to session data with eyemovements
loading_data_grouping_decoding5.m

get_subcellcrcon_decoding
% test cross-contrast or contrast-independent models in selected cells
% in addition, it can be used to test the same-contrast data.
% therefore, with cellselopt.CELLSEL_THRs = 0, will result in the same results of get_crcon_decoding.m 

get_subcellcrcon_subdata_decoding
% simply same as get_subcellcrcon_decoding
% but getting decoding accuracy from e.g.,
% DATA_DISK_UNION_CONTRSP_ctm0.60ses16-P1.mat


get_crcon_decoding.m
% test cross-contrast or contrast-independent models

get_respstat.m
% collect all data together to present the basic statistics such as mean and variance 

get_datawithPA.m
% collect data based on mean population response levels

get_ORI_conscale_PADEP.m
% fitting ORI tuning function at one contrast and PADEP (e.g., 100,H) to
% those at other contrasts and PA (e.g., 100,L)


get_ORI_conscale2.m get_ORI_conscale_subdata2.m--- function
% fitting ORI tuning function at 100% contrast to those at other contrasts
% when imaging scans are acquired at different time points, based on the
% timeline of acquisition, the fitting can be separated with
% get_ORI_conscale_subdata.m

get_ORI_conscale3
% data fitting within contrast within session
get_ORI_conscale4
%data fitting within/across contrast across sessions -- fitting between time-separate session data 



%% re-validate fit_ORIscale with get_ORI_conscale2, get_ORI_conscale3, and get_ORI_conscale_PADEP
% 1) because old code the smallest lambda with which explained variance of fitting is maximized.
% 2) In this test, I chose the largest lambda with which explained variance of fitting is maximized.
% % particularly, I change the code in "cal_wori.m"
% opts.lams =  10.^(linspace(-3,1,20));
% into
% opts.lams =  fliplr(10.^(linspace(-5,1,20)))*T;
% These results were saved into
% L_AN1-16_ORIsc_ctm0.60_fit_ab4.mat with get_ORI_conscale2.m
% L_WSWC_AN1-16_ORIsc_ctm0.60_fit_ab4.mat with get_ORI_conscale3.m





%---------------- used before 08/29/16
% loading scan data to session data
loading_data_grouping_decoding4.m

get_decoding_cellsel.m
% perform SMLR decoding L1 or L2 regularization
% and select cell based on weight
% Due to this reason, data are normalized only to use weight for the
% selection
% when sub-data (i.e., parts of all sessions) are used, cells commonly
% selected across sessions are used for future consistency
% of cell selection



get_crossCONdecoding_cellsel.m
% perform orientation decoding across contrast
% to check contrast-independent decoder


get_SUBCELL_crossCONdecoding.m
% perform orientation decoding across contrast
% to check contrast-independent decoder
% additionally, perform orientation decoding on cells selected by smlr
% particularly, the output from "get_crossCONdecoding_cellsel.m" were used.

plot_subcell_decoding
% generate summary tables from the outputs generated from
% "get_SUBCELL_crossCONdecoding.m"


get_subcell_decoding_multises.m
% cross-validation between two-consecutive sessions.
% 1. select cells from one sessions 
% 2. use the selected cells for cross-validation of decoding accuracies in
% the other session.
% prior to call this function, the output from get_decoding_cellsel.m is
% required. 

get_group_decoding
% generate all possible combinations in a given "n"-cell group
% Then, calculate decoding accuracies



get_group_decoding_summary.m
% after performing get_group_decoding, collect the results across sessions
% and make summary statistics.

get_decscore.m
% after performing get_group_decoding, score the contribution of individual
% cells in n-cell-group decoding.
% briefly, the scores of individual groups are calculated with  (ACC -50)/(maxdec-50) 
% and then averaged over all groups, where a cell always participates.

get_subcell_crsdecoding.m
% perform classification cross sessions by using selected cells and
% training classifiers in a session and then applying the classifiers into
% the remaining session.
% This code will check the preservance of coding model across sessions and
% across contrasts




plot_ori_grpdec.m, plot_hist_grpdec.m
% plot orientation tunings of cells selected from get_group_decoding and
% get_group_decoding_summary



get_ORI_conscale.m, get_ORI_conscale_subdata.m -- old files

%% manual eyetracking memo

11 right eye not blocked

    0.7800    0.7700    0.7240    0.8267    0.7950    0.5046
    0.7421    0.7239    0.6875    0.6206    0.7551    0.5211
    0.7279    0.7806    0.6105    0.7296    0.6755    0.6188
    0.7465    0.6667    0.6702    0.7007    0.7205    0.5229    

23 eyebrow clipped, small pupil in some scans,

    0.8900    0.8600    0.8425    0.8712    0.8475    0.6511
    0.7266    0.6265    0.7016    0.7628    0.8333    0.6000
    0.8920    0.8388    0.8263    0.8120    0.7820    0.7197
    0.6786    0.6803    0.5562    0.7200    0.7219    0.6333

    
25 eyebrow clipped, not very obvious for pupil border, eye trembling occasionally

    0.9392    0.9350    0.9100    0.9479    0.9842    0.8343
    0.8214    0.8043    0.6970    0.8557    0.7757    0.6234
    0.8950    0.8917    0.8350    0.8567    0.9083    0.7583
    0.7317    0.6557    0.6529    0.6803    0.6622    0.5198
    
There was bug in function setinvalid of reest_gui.m.
This bug was found in 8/8/2016 and thus 
the manual corrections made before 10:29AM of 08/08/2016 
have not correct info in fitout.out.sPIX( but all others are correct).
Therefore, the length of fitout.out.sPIX is smaller than other info
    
26 good eyetracking
    0.7640    0.8056    0.7758    0.7105    0.7679    0.6127
    0.4441    0.5761    0.7035    0.6074    0.5889    0.5666
    0.9350    0.9010    0.9430    0.7950    0.9310    0.6694
    0.6374    0.6720    0.6732    0.6515    0.6313    0.4600
    
    
27 good eyetracking, but pupil size is too small
    0.9100    0.8850    0.8409    0.9000    0.8600    0.5799
    0.6235    0.7733    0.7200    0.8571    0.6995    0.5479
    0.8298    0.7347    0.5479    0.7267    0.7356    0.5797
    0.6628    0.7000    0.6923    0.6605    0.5897    0.5602
    
28 substantially small pupil size, probably needed to cut data with small pupil, say r=10
I will abandon this data set because the pupil size is significantly smaller since the 2nd scan
    0.7357    0.8036    0.6850    0.7957    0.7057    0.5661
    0.5812    0.6257    0.5256    0.6457    0.5730    0.4762
    0.7412    0.6637    0.7225    0.8037    0.8150    0.7816
    0.6037    0.6101    0.5440    0.6046    0.6279    0.5071    
    
29-best decoding, good eyetracking, but pupil suface is not tangential to eyecamera


30- good eyetracking, 3rd scan missing eyetracking, 8th scan has very small pupil and difficult to identify
% this session might need to have pupil-size-cut-off criterion
% The 8th scan can be potentially completely-abandoned
    0.7857    0.7300    0.7193    0.8271    0.7914    0.6719
    0.8364    0.8450    0.7000    0.8393    0.8264    0.6530
    0.8093    0.7975    0.5997    0.7612    0.6698    0.5859
    0.7600    0.7475    0.5491    0.7222    0.6139    0.5240

    

32- good eyetracking, 
difficulty on automatic separation of 7th scan eyetracking (008.mj2) due to shade,    
Let's decide whether I include 7th scan data in the analysis.
Currently, it is out in the data selected.
    0.9306    0.9350    0.8413    0.9413    0.8344    0.7275
    0.8063    0.8294    0.7013    0.8444    0.6675    0.5874    
    0.8467    0.8983    0.8636    0.8550    0.8043    0.7500
    0.7900    0.8057    0.6905    0.7821    0.6443    0.6622
    
    
    
33   
difficulty on automatic separation of 2nd scan eyetracking  due to shade, 
I excluded 8th scan due to small pupil and difficulty of eyetracking
    0.9600    0.9400    0.8208    0.9236    0.8467    0.5922
    0.8536    0.8650    0.6406    0.8736    0.7371    0.5389
    0.9575    0.9419    0.8294    0.9306    0.8550    0.6333
    0.8137    0.8406    0.6484    0.8106    0.6556    0.4679
    
    
35 pupil size is large, but eyetracking is hard, I have not processsed yet
    0.7329    0.8700    0.6888    0.8536    0.7547    0.5325
    0.6523    0.7237    0.5208    0.7561    0.5637    0.4863
    0.8580    0.9310    0.8458    0.9260    0.8290    0.5158
    0.7258    0.8583    0.7092    0.8190    0.6457    0.5488
    
    
36
    0.6354    0.7310    0.6917    0.7702    0.7408    0.5418
    0.6046    0.7370    0.6341    0.7575    0.6199    0.5556
    0.6812    0.7856    0.6925    0.7750    0.6572    0.5589
    0.5320    0.5669    0.5511    0.5361    0.5332    0.5288 
    
    
40 - good eyetracking
    0.8333    0.8850    0.8783    0.7985    0.7803    0.7005
    0.6997    0.7298    0.6701    0.7717    0.6111    0.5679
    0.7394    0.8512    0.8325    0.8930    0.8650    0.8037
    0.7423    0.7230    0.7227    0.7580    0.6415    0.6054    


