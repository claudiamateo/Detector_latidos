% Asignatura Tratamiento de Señales
% Autoras: Lucia Herrador Dominguez
%          Claudia Mateo Burillo

%% Cargamos la señal escogida

% Leemos la señal 118e_6 de 30K muestras de Physionet
[sign,Fs,tm] = rdsamp('nstdb/119e06',[],129999,100000);