function [Cpos,Cneg,Npos,Nneg,Kpos,Kneg] = fnc_density(W)
% DENSITY_UND        Density
%
%   kden = density_und(W);
%   [kden,N,K] = density_und(W);
%
%   Density is the fraction of present connections to possible connections.
%
%   Input:      W,    undirected (weighted/binary) connection matrix
%
%   Output:     Cpos,   density positives
%               Cneg,   density negatives
%               Npos,      number of vertices pos
%               Kpos,      number of edges pos
%               Nneg,      number of vertices neg
%               Kneg,      number of edges neg
%
%   Notes:  Assumes CIJ is undirected and has no self-connections.
%           Weight information is discarded.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

W_pos    = W.*(W>0);
[Cpos,Npos,Kpos] = density_und(W_pos);

W_neg    = -W.*(W<0);
[Cneg,Nneg,Kneg] = density_und(W_neg);             