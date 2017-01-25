    %--------- collecting data
function [inxsample, cons, oris] = collect(data)


        inxs_valtrial = data.events(:)>0 & ~isinf(data.events(:));
        unique_evt = unique(data.events(inxs_valtrial)');
        unique_evt = setdiff(unique_evt,0);
        nevt = length(unique_evt);
        evts = data.events(:);
        inxsample = cell(1,nevt);
        cons = zeros(1,nevt);
        oris = zeros(1,nevt);
        for ievt0 = 1:nevt         
            ievt = unique_evt(ievt0);
            inxsample{ievt0} = TP.select_subdata(evts,{ievt});
            cons(ievt0)=data.events_cont(inxsample{ievt0}(1));
            oris(ievt0)=data.events_ORI(inxsample{ievt0}(1));
        end
end