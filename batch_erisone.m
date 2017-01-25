cellselopt.CELLSEL_THRs=[0 1 3 5 10 20];
cellselopt.SEL_METHOD='ncell';
comcont={{100,100},{40,40},{100, 40},{40,100},{[100 40],100},{[100 40],40}};
DATA_thr_str = 'thr5_eyethr0p3'
iexp_type=5;

get_subcellcrcon_decoding_erisone(iexp_type, DATA_thr_str,...
    comcont, cellselopt)
