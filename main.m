% Asignatura Tratamiento de Señales
% Autoras: Lucía Herrador Domínguez
%          Claudia Mateo Burillo

%% Cargamos la señal escogida

clc, close all, clear all;
% Leemos la señal 119e06 de Physionet
[sign,Fs,tm] = rdsamp('nsrdb/16265',[],3840,0);


frecuencia_de_muestreo = 1/Fs;

%Se quieren 30 segundos de la señal. Miramos lo que equivale en muestras de
%la señal.
conversion_a_segundos_final = 30;
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

%% 3.1. Filtro paso bajo para eliminar altas frecuencias


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

%Probamos con un filtro cheby2 de orden 2
[bc2_2,ac2_2] = cheby2(2,80,Wcb); % 80 son los decibelios en la banda de rechazo

figure('Name','Respuesta del sistema en frecuencia filtro Chebyshev 2'),
freqz(bc2_2,ac2_2),title('Respuesta en frecuencia del filtro Chebyshev 2');

SSC_filtrada_2_cheby2=filter(bc2_2,ac2_2,SSC);
figure('Name','Señal vs Señal filtrada paso bajo con Chebyshev 2'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,SSC_filtrada_2_cheby2),
title('Señal filtrada paso bajo con filtro Chebyshev 2 de orden 2');

%Probamos con un filtro cheby 2 de orden 12
[bc2_12,ac2_12] = cheby2(12,80,Wcb); % 80 son los decibelios en la banda de rechazo

figure('Name','Respuesta del sistema en frecuencia filtro Chebyshev 2'),
freqz(bc2_12,ac2_12),title('Respuesta en frecuencia del filtro Chebyshev 2');

SSC_filtrada_12_cheby2=filter(bc2_12,ac2_12,SSC);
figure('Name','Señal vs Señal filtrada paso bajo con Chebyshev 12'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,SSC_filtrada_12_cheby2),
title('Señal filtrada paso bajo con filtro Chebyshev 2 de orden 12');

%Mostramos la diferencia entre los filtros Chebyshev 2 de diferente orden 
figure('Name','Señal vs Señal filtrada paso bajo con filtro Chebyshev 2'),
subplot(411),plot(tm,SSC),title('Señal original')
subplot(412),plot(tm,SSC_filtrada_cheby2),title('Señal filtrada paso bajo con filtro elíptico orden 6'),
subplot(413),plot(tm,SSC_filtrada_2_cheby2),
title('Señal filtrada paso bajo con filtro Chebyshev 2 orden 2'),
subplot(414),plot(tm,SSC_filtrada_12_cheby2),
title('Señal filtrada paso bajo con filtro Chebyshev 2 orden 12');




%Probamos con el filtro elíptico de orden 6
[be,ae] = ellip(6,5,80,Wcb);

figure('Name','Respuesta del sistema en frecuencia filtro elíptico orden 6'),
freqz(be,ae),title('Respuesta en frecuencia del filtro Elíptico orden 6');

FPB_E = filter(be,ae,SSC); 
figure('Name','Señal vs Señal filtrada paso bajo con filtro elíptico orden 6'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPB_E),
title('Señal filtrada paso bajo con filtro elíptico orden 6');

%Probamos con un filtro elíptico de orden 2
[be2,ae2] = ellip(2,5,80,Wcb);

figure('Name','Respuesta del sistema en frecuencia filtro elíptico orden 2'),
freqz(be2,ae2),title('Respuesta en frecuencia del filtro Elíptico orden 2');

FPB2_E = filter(be2,ae2,SSC); 
figure('Name','Señal vs Señal filtrada paso bajo con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPB2_E),
title('Señal filtrada paso bajo con filtro elíptico orden 2');

%Probamos con Filtro elíptico de orden 12.
[be3,ae3] = ellip(12,5,80,Wcb);
figure();
freqz(be3,ae3), title('Respuesta en frecuencia del filtro Elíptico');

FPB3_E = filter(be3,ae3,SSC); 
figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico de orden 12'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPB3_E),
title('Señal filtrada paso bajo con filtro elíptico de orden 12');

%Mostramos la diferencia entre los filtros elípticos 2 del mismo orden 
figure('Name','Señal vs Señal filtrada paso bajo con filtro elíptico'),
subplot(411),plot(tm,SSC),title('Señal original')
subplot(412),plot(tm,FPB_E),title('Señal filtrada paso bajo con filtro elíptico orden 6'),
subplot(413),plot(tm,FPB2_E),title('Señal filtrada paso bajo con filtro elíptico orden 2');
subplot(414),plot(tm,FPB3_E),title('Señal filtrada paso bajo con filtro elíptico orden 12');


%% 3.2. Filtro paso alto para eliminar la componente de continua y para eliminar las
% bajas frecuencias 


%Probamos con el filtro paso alto implementado en la presentación Prezi
a1 = [1 -1];
b1 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1];
SSC_filtrada_paso_alto = filter(b1,a1,SSC);
SSC_filtrada_paso_alto = SSC_filtrada_paso_alto/32;
SSC_filtrada_paso_alto = SSC-SSC_filtrada_paso_alto;

figure, subplot(2,1,1), plot(tm,SSC,'b'),title('Señal ECG original'),
subplot(2,1,2),plot(tm,SSC_filtrada_paso_alto,'r');
title('Señal ECG con filtro paso alto presentación de frec de corte de 5Hz');

%Calculamos la frecuencia de corte normalizada de los filtros
Fca = 5; % Frecuencia de corte (en Hz)
Wca = Fca/(Fs/2); % Frecuencia de corte normalizada (en radianes)



%Probamos con filtro Chebyshev 2 de orden 6.
[bc2_PA,ac2_PA] = cheby2(6,80,Wca,'high'); % 80 son los decibelos en la banda de rechazo
figure; 
freqz(bc2_PA,ac2_PA), title('Respuesta en frecuencia del filtro PA Chebysehv 2 orden 6');
title('Filtro Chebyshev 2');
FPA_C2 = filter(bc2_PA,ac2_PA,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_C2),
title('Señal filtrada paso alto con filtro Chebyshev 2 de orden 6');

%Probamos con un filtro cheby2 de orden 2
[bc2_PA_2,ac2_PA_2] = cheby2(2,80,Wca,'high'); % 80 son los decibelos en la banda de rechazo
figure; 
freqz(bc2_PA_2,ac2_PA_2), title('Respuesta en frecuencia del filtro PA Chebysehv 2 orden 2');
title('Filtro Chebyshev 2');
FPA_C2_2 = filter(bc2_PA_2,ac2_PA_2,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_C2_2),
title('Señal filtrada paso alto con filtro Chebyshev 2 de orden 2');

%Probamos con un filtro cheby2 de orden 12
[bc2_PA_12,ac2_PA_12] = cheby2(12,80,Wca,'high'); % 80 son los decibelos en la banda de rechazo
figure; 
freqz(bc2_PA_12,ac2_PA_12), title('Respuesta en frecuencia del filtro PA Chebysehv 2 orden 12');
title('Filtro Chebyshev 2');
FPA_C2_12 = filter(bc2_PA_12,ac2_PA_12,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_C2_12),
title('Señal filtrada paso alto con filtro Chebyshev 2 de orden 12');

%Mostramos la diferencia entre los filtros Chebyshev 2 de diferente orden 
figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
subplot(411),plot(tm,SSC),title('Señal original')
subplot(412),plot(tm,FPA_C2),title('Señal filtrada paso alto con filtro Chebyshev 2 orden 6'),
subplot(413),plot(tm,FPA_C2_2),
title('Señal filtrada paso alto con filtro Chebyshev 2 orden 12'),
subplot(414),plot(tm,FPA_C2_12),
title('Señal filtrada paso alto con filtro Chebyshev 2 orden 2');




%Probamos con Filtro elíptico de orden 6.
[be_PA,ae_PA] = ellip(6,5,80,Wca,'high');
figure();
freqz(be_PA,ae_PA), title('Respuesta en frecuencia del filtro Elíptico');
FPA_E = filter(be_PA,ae_PA,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_E),
title('Señal filtrada paso alto con filtro elíptico de orden 6');

%Probamos con un filtro elíptico de orden 2
[be_PA2,ae_PA2] = ellip(2,5,80,Wca,'high');
figure();
freqz(be_PA2,ae_PA2), title('Respuesta en frecuencia del filtro Elíptico');
FPA_E2 = filter(be_PA2,ae_PA2,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_E2),
title('Señal filtrada paso alto con filtro elíptico de orden 2');

%Probamos con Filtro elíptico de orden 12.
[be_PA3,ae_PA3] = ellip(12,5,80,Wca,'high');
figure();
freqz(be_PA3,ae_PA3), title('Respuesta en frecuencia del filtro Elíptico');
FPA_E3 = filter(be_PA3,ae_PA3,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_E3),
title('Señal filtrada paso alto con filtro elíptico de orden 12');

%Mostramos la diferencia entre los filtros elípticos 2 de diferente orden 
figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
subplot(411),plot(tm,SSC),title('Señal original')
subplot(412),plot(tm,FPA_E),title('Señal filtrada paso alto con filtro elíptico orden 6'),
subplot(413),plot(tm,FPA_E2),title('Señal filtrada paso alto con filtro elíptico orden 2'),
subplot(414),plot(tm,FPA_E3),title('Señal filtrada paso alto con filtro elíptico orden 12');



%Señal final con los filtros seleccionados
%Elegimos el FP presentación Prezi
soriginal_filtrada_PB=SSC_filtrada_paso_bajo;
figure, plot(tm,soriginal_filtrada_PB);
%Filtramos esa señal con el filtro elíptico PA de orden 2
sfinal=filter(b1,a1,soriginal_filtrada_PB); 
figure, plot(tm,sfinal);

%Para la presentación de Prezi añadimos esto
sfinal = sfinal/32;
sfinal = soriginal_filtrada_PB-sfinal;
figure, plot(tm,sfinal);

%% Actividad 4: Derivación
%Derivamos sfinal con diff
D_sfinal=diff(sfinal);
figure(),plot(D_sfinal);
D_sfinal_cuadrado=(D_sfinal).^2;
figure(),plot(D_sfinal_cuadrado);


%% Actividad 5: Detectar los QRS y marcarlos

