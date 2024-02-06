function [dd] = fncHank(Xdr)
% dat es la matriz de la venta.
dd = cell(1,numel(Xdr));
for s = 1:numel(Xdr)
    sample = numel(Xdr{s});
    dd{s} = cell(sample,1);
    for j = 1:sample
        xv = Xdr{s}{j};
        xh = hankel(xv,[xv(1),xv(1:end-1)]);
        ddat = hankel(Xdr{s}{j},[Xdr{s}{j}(1),Xdr{s}{j}(1:end-1)]); % realiza con el vector la hankelizacion.
        num = size(ddat,1);
        dd{s}{j} = ddat(1:num,1:num);
    end
end

 ddat = hankel(Xdr{s}{j}',[Xdr{s}{j}(end)' Xdr{s}{j}']); %