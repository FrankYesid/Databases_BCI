function [X T] = FilterWave(x,Wave,lv)
% Esta funcion devuelve los canales filtrados por cada uno de los 2^lv filtros
% Wave: wavelet madre 'sym2';
% lv: nivel de descomposicion  3;

if size(x,1)==1 | size(x,2)==1
    x = x(:);
end

[N N_canal]= size(x);
% Inicializar matriz
% X = cell(N_canal,1);
% N_fil = sum(2.^(1:lv));%2^lv;
% X(:) = {zeros(N,N_fil)};
% X = zeros(N_fil,N_canal,N);
X = cell(lv,1);

% aplicar la desocposiccion de wavelet
for l = 1:lv
    % filtrar se?al
    X{l} = zeros(2^l,N_canal,N);
    for cnl = 1:N_canal        
        T = wpdec(x(:,cnl),l,Wave);
        % devuelve los nodos en orden de frecuencia (ten_freq) y en orden paley % (tn_pal)
        [~,tn_freq] = otnodes(T);        
        for fil = 1:2^l            
            % reconstruye el nodo fil
            X{l}(fil,cnl,:) = wprcoef(T,tn_freq(fil));
        end
    end
end
X=cell2mat(X);
