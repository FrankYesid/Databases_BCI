% Esta funcion selecciona el mejor arbol de descomposicion seguna la
% energia de cada nodo
function [BestRec,BestT] = f_BestTree(Xrec,Coef)
% Xrec: es una celda que contiene la se?al filtrada en cada nivel Xrec{lv}(time,fil,canal)
% Coef: es una celda que contiene los coeficientes wavelets de cada nivel Coef{lv}(time,fil,canal)

N_time = size(Xrec{1},1);
[~,lv,cnl,N_tri] = size(Coef{end});
lv = log2(lv);

% calcular la energia de cada nodo
E = cellfun(@(x) sum(sum(sum(x.^2,1),3),4), Coef,'UniformOutput',false);

% inicializar matrices
EMatrix = zeros(lv,2^lv);
for i_lv = 1:lv
    EMatrix(i_lv,:) = [E{i_lv} zeros(1,2^lv-2^i_lv)];
end
EMatrix(1,:)= EMatrix(1,:)/sum(EMatrix(1,:));
Mat = zeros(lv,2^lv); 
Mat(1,1:2) = 1;

% Construir el BestTree, mediante el siguiente criterio:
% si el nodo que se analiza tiene una energia mayor a un umbral se divide,
% en otro caso se elimina.
for i_lv = 1:lv-1
    index = 1;
    for node = 1:2^lv
        if EMatrix(i_lv,node) > 0.1
            Mat(i_lv,node) = 0;
            Mat(i_lv+1,index) = 1; index = index + 1;
            Mat(i_lv+1,index) = 1; index = index + 1;
        else
            index = index + 2;
        end
    end
    % Normalizar la Enegia del nodo
    temp = EMatrix(i_lv+1,:).*Mat(i_lv+1,:);
    EMatrix(i_lv+1,:) = temp/sum(temp);
end

% BestTree por niveles: BestT(lv)(time,filtro,canal,trial)
BestT = cell(lv,1);
[N_time,N_fil,~,N_trial] = size(Xrec{1});
BestT(:) = {zeros(N_time,N_fil,cnl,N_trial)};

for i = 1:lv
    ind = find(EMatrix(i,:)>0);
    BestT{i,1} = Xrec{i}(:,ind,:,:);
    fprintf(['lv: ' num2str(i) ' de ' num2str(lv)   '\n'])
end

% BestTree por niveles: BestT(trial,filtro,canal,tiempo)

temp = cell(size(BestT,1),1);

for j = 1:size(BestT,1)
    Temp = zeros(size(BestT{j},2),cnl,N_time,N_trial);
    for y = 1:N_trial
        for i = 1:size(BestT{j},2)
            tt = squeeze(BestT{j,1}(:,i,:,y))';
            Temp(i,:,:,y) = tt;
        end
        temp{j,1} = Temp;
    end
end
BestT = cell2mat(temp);


% Nodos para la reconstruccion ordenados de menor a mayor frecuencia
[row, col] = find(Mat>0);
[~,ind] = sort(row,'descend');
coor = [row(ind) col(ind)];
BestRec = zeros(N_time,length(row),cnl,N_tri);
for i = 1:length(row)
%     fprintf(['lv: ' num2str(i) ' de ' num2str(length(row))   '\n'])
    BestRec(:,i,:,:) = Xrec{coor(i,1)}(:,coor(i,2),:,:);
end