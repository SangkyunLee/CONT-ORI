function [D1,D2] = subdata(subinx,D1,D2,dtyps,dr)
% function [D1,D2] = subdata(M,D1,D2,dtyps)
% retrieve subdata with subinx in given dtyps
% NOTE: it does not retrieve all data automatically, 
% it retrive data only with dyps (cell type variable)
% dr specifies dimension inx of Dx.(dyps{ityp}) 
    
    if nargin<5
        dr = 2;
    end
      
    inx1 = D1.cellinx_sel;
    inx2 = D2.cellinx_sel;

    [~,~,ni1]=intersect(subinx,inx1);
    [~,~,ni2]=intersect(subinx,inx2);

    D1.cellinx_sel =D1.cellinx_sel(ni1);
    D2.cellinx_sel =D2.cellinx_sel(ni2);
    for ityp = 1 : length(dtyps)
        dtyp = dtyps{ityp};
        if dr==1,
            D1.(dtyp) = D1.(dtyp)(ni1,:);
            D2.(dtyp) = D2.(dtyp)(ni2,:);
        elseif dr==2,
            D1.(dtyp) = D1.(dtyp)(:,ni1);
            D2.(dtyp) = D2.(dtyp)(:,ni2);
        end
    end
end