function [SI] = calcSymmetryIndex(L,R)
% Fonction qui calcule le Symmetry Index entre 2 jambes, sans postulat de
% prédominance de jambe, etc..
%   L = jambe gauche, R = jambe droite
SI = (abs(L - R) / (0.5 * (L + R))) * 100;
end