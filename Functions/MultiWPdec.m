function [X X_coef]= MultiWPdec(x1,Wave,lv)
% Esta funcion devuelve los canales filtrados por cada uno de los 2^lv filtros
% Wave: wavelet madre 'sym2';
% lv: nivel de descomposicion  3;
% X: contiene la informacion de los filtros de menor a mayor frecuencia asi:
% X{NivelDescomposicion}(tiempo,filtro,canal)

N_tri = length(x1);
[N N_canal]= size(x1{1});
% Inicializar matriz
% X = cell(N_canal,1);
% N_fil = 2^lv;
% X(:) = {zeros(N,N_fil)};

X = cell(lv,1);
X_coef = cell(lv,1);

% filtrar se?al
for tri = 1:N_tri
    x=x1{tri};
    for lev = lv
        N_fil = 2^lev;
        for cnl = 1:N_canal
            fprintf(['trial: ' num2str(tri) ' de ' num2str(N_tri) '...lv: ' num2str(lev) ' de ' num2str(lv) ' ...channel: ' num2str(cnl) ' de ' num2str(N_canal) '\n'])
            % aplicar la desocposiccion de wavelet
            T = wpdec(x(:,cnl),lev,Wave);
            % devuelve los nodos en orden de frecuencia (ten_freq) y en orden paley % (tn_pal)
            [tn_pal,tn_freq] = otnodes(T);
            for fil = 1:N_fil
                
                % reconstruye el nodo fil
                X{lev}(:,fil,cnl,tri) = wprcoef(T,tn_freq(fil));
                X_coef{lev}(:,fil,cnl,tri) = wpcoef(T,tn_freq(fil));
            end
        end
    end
end
