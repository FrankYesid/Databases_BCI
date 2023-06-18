function [dd] = fncHankeling_(Xdr,e)

tr = numel(Xdr);
%     dd{s} = cell(sample,1);
for j = 1:tr
    x = Xdr{j};
    dd{j} = hankel(x(1:e,1),x(e:end,1)); % realiza con el vector la hankelizacion.
    %dd{j} = hankel(x(1:e),[x(e:end);x(1:e-1)]); % realiza con el vector la hankelizacion.
    %         num = size(ddat,2);
    %dd{s}{j} = ddat(:,1:end-1);
end
