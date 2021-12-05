% Asignatura Tratamiento de Señales
% Autoras: Lucía Herrador Domínguez
%          Claudia Mateo Burillo

%% Cargamos la señal escogida

clc, close all, clear all;
% Leemos la señal 119e06 de Physionet
[sign,Fs,tm] = rdsamp('nstdb/119e06',[],118800,108000);

%Elegimos el tamaño de la muestra. Se sabe que el ruido comienza en el
%minuto 5
frecuencia_de_muestreo = 1/Fs;
conversion_a_segundos = 5*60;
primera_muestra_con_ruido = conversion_a_segundos/frecuencia_de_muestreo;

% Como 1/Fs es 0.028 y el ruido comienza pasados 5 minutos segun la BBDD
% de Noise Stress Test, multiplicamos ese tiempo por 60 para pasarlo a
% segundos y lo dividimos entre la frecuencia de muestre obtenida. De
% esta manera, obtenemos un valor (108000) que representa en qué segundo 
% se produce la primera muestra con ruido.

%Se quieren 30 segundos de la señal. Miramos lo que equivale en muestras de
%la señal.
conversion_a_segundos_final = 5.5*60;
ultima_muestra = conversion_a_segundos_final/frecuencia_de_muestreo;

%Se trabajará solo con una derivación de la señal
soriginal= sign(:,1);

%% Actividad 2. Pre-procesar la señal
%2.1. Visualizamos la señal en tiempo y frecuencia. Visualizamos en el
%dominio del tiempo.
figure('Name', 'Señal original en el dominio del tiempo'),
plot(tm, soriginal), title('Señal original');

%Visualizamos en  el dominio de la frecuencia
TF_sign = fftshift(fft(sign));
n = length(TF_sign);
fshiftaxis = linspace(-Fs/2,Fs/2,n);
figure("Name",'Señal original en el dominio de la frecuencia');
plot(fshiftaxis,abs(TF_sign),'b');
title('TF de la señal original');


%2.2. Quitamos la componente continua a la señal
SSC= soriginal(:,1)-mean(soriginal(:,1));
figure('Name','Señal original sin componente continua'), plot(tm,SSC),
title('Señal original sin componente continua');

%Visualizamos en el dominio de la frecuencia la señal ahora sin componente
%continua
TF_SSC = fftshift(fft(SSC));
n2 = length(TF_SSC);
fshiftaxis = linspace(-Fs/2,Fs/2,n2);
figure("Name",'Señal original en el dominio de la frecuencia');
plot(fshiftaxis,abs(TF_SSC),'b');
title('TF de la señal sin continua');

%% 3. Filtrado de artefactos en la señal ECG

%3.1. Filtro paso bajo para eliminar altas frecuencias
%Probamos con el filtro paso bajo implementado en la presentación Prezi
a = [1 -2 1];
b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
SSC_filtrada_paso_bajo = filter(b,a,SSC);

figure, subplot(2,1,1), plot(tm,SSC,'b'),title('Señal ECG original'),
subplot(2,1,2),plot(tm,SSC_filtrada_paso_bajo,'r');
title('Señal ECG con filtro paso bajo presentación de frec de corte de 15Hz');

%Frecuencia de los filtros
Fcb = 15; % Frecuencia de corte (en Hz)
Wcb = Fcb/(Fs/2); % Frecuencia de corte normalizada (en radianes)

%Probamos con un filtro cheby2 de orden 6
[bc2,ac2] = cheby2(6,80,Wcb); % 80 son los decibelios en la banda de rechazo

figure('Name','Respuesta del sistema en frecuencia filtro Chebyshev 2'),
freqz(bc2,ac2),title('Respuesta en frecuencia del filtro Chebyshev 2');

SSC_filtrada_cheby2=filter(bc2,ac2,SSC);
figure('Name','Señal vs Señal filtrada paso bajo con Chebyshev 2'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,SSC_filtrada_cheby2),
title('Señal filtrada paso bajo con Chebyshev 2');

%Probamos con un filtro cheby2 de orden X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mostramos la diferencia entre los filtros Chebyshev 2 del mismo orden 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Probamos con el filtro elíptico de orden 6
[be,ae] = ellip(6,5,80,Wcb);

figure('Name','Respuesta del sistema en frecuencia filtro elíptico'),
freqz(be,ae),title('Respuesta en frecuencia del filtro Elíptico');

FPB_E = filter(be,ae,SSC); 
figure('Name','Señal vs Señal filtrada paso bajo con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPB_E),
title('Señal filtrada paso bajo con filtro elíptico');

%Probamos con un filtro elíptico de orden X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mostramos la diferencia entre los filtros elípticos 2 del mismo orden 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%3.2. Filtro paso alto para eliminar la componente de continua y para eliminar las
% bajas frecuencias 

%Probamos con el filtro paso alto implementado en la presentación Prezi
a1 = [1 - 1];
b1 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1];
SSC_filtrada_paso_alto = filter(a1,b1,SSC);
SSC_filtrada_paso_alto = SSC_filtrada_paso_alto/32;
SSC_filtrada_paso_alto = SSC-SSC_filtrada_paso_alto;

figure, subplot(2,1,1), plot(tm,SSC,'b'),title('Señal ECG original'),
subplot(2,1,2),plot(tm,SSC_filtrada_paso_alto,'r');
title('Señal ECG con filtro paso alto presentación de frec de corte de 5Hz');

%Calculamos la frecuencia de corte normalizada de los filtros
Fca = 5; % Frecuencia de corte (en Hz)
Wca = Fca/(Fs/2); % Frecuencia de corte normalizada (en radianes)

%Probamos con filtro Chebyshev 2 de orden 6.

%Probamos con un filtro cheby2 de orden X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mostramos la diferencia entre los filtros Chebyshev 2 del mismo orden 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Probamos con Filtro elíptico de orden 6.
[be2,ae2] = ellip(6,5,80,Wca,'high');
figure();
freqz(be2,ae2), title('Respuesta en frecuencia del filtro Elíptico');
FPA_E = filter(be2,ae2,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_E),
title('Señal filtrada paso alto con filtro elíptico');

%Probamos con un filtro elíptico de orden X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mostramos la diferencia entre los filtros elípticos 2 del mismo orden 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Actividad 4: Derivación

%% Actividad 5: Detectar los QRS y marcarlos

