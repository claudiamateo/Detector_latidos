% Asignatura Tratamiento de Señales
% Autoras: Lucía Herrador Domínguez
%          Claudia Mateo Burillo

%% Cargamos la señal escogida

clc, close all, clear all;

% Leemos la señal de Physionet
% Registro 1
%[sign,Fs,tm] = rdsamp('nsrdb/16265',[],3840,0);
% Registro 2
%[sign,Fs,tm] = rdsamp('nsrdb/16272',[],3840,0);
% Registro 3
%[sign,Fs,tm] = rdsamp('nsrdb/16273',[],3840,0);
% Registro 4
%[sign,Fs,tm] = rdsamp('nsrdb/16420',[],3840,0);
% Registro 5
%[sign,Fs,tm] = rdsamp('nsrdb/16483',[],3840,0);
% Registro 6
%[sign,Fs,tm] = rdsamp('nsrdb/16539',[],3840,0);
% Registro 7
[sign,Fs,tm] = rdsamp('nsrdb/16773',[],3840,0);
% Registro 8
%[sign,Fs,tm] = rdsamp('nsrdb/16786',[],3840,0);
% Registro 9
%[sign,Fs,tm] = rdsamp('nsrdb/16795',[],3840,0);
% Registro 10
%[sign,Fs,tm] = rdsamp('nsrdb/17052',[],3840,0);


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

%Normalizamos la señal
mins=min(SSC);
SSC=SSC-mins; %El punto mínimo será 0
maxs=max(SSC);
SSC=SSC/maxs; % El valor máximo ahora valdrá 1 y los demás valores estarán entre 0 y 1.

%% 3. Filtrado de artefactos en la señal ECG

%% 3.1. Filtro paso bajo para eliminar altas frecuencias


%Probamos con el filtro paso bajo implementado en la presentación Prezi
% a = [1 -2 1];
% b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
% SSC_filtrada_paso_bajo = filter(b,a,SSC);
% 
% figure, subplot(2,1,1), plot(tm,SSC,'b'),title('Señal ECG original'),
% subplot(2,1,2),plot(tm,SSC_filtrada_paso_bajo,'r');
% title('Señal ECG con filtro paso bajo presentación de frec de corte de 15Hz');



%Frecuencia de los filtros
Fcb = 15; % Frecuencia de corte (en Hz)
Wcb = Fcb/(Fs/2); % Frecuencia de corte normalizada (en radianes)

% %Probamos con un filtro cheby2 de orden 6
% [bc2,ac2] = cheby2(6,80,Wcb); % 80 son los decibelios en la banda de rechazo
% 
% figure('Name','Respuesta del sistema en frecuencia filtro Chebyshev 2'),
% freqz(bc2,ac2),title('Respuesta en frecuencia del filtro Chebyshev 2');
% 
% SSC_filtrada_cheby2=filtfilt(bc2,ac2,SSC);
% figure('Name','Señal vs Señal filtrada paso bajo con Chebyshev 2'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,SSC_filtrada_cheby2),
% title('Señal filtrada paso bajo con Chebyshev 2');

% %Probamos con un filtro cheby2 de orden 2
% [bc2_2,ac2_2] = cheby2(2,80,Wcb); % 80 son los decibelios en la banda de rechazo
% 
% figure('Name','Respuesta del sistema en frecuencia filtro Chebyshev 2'),
% freqz(bc2_2,ac2_2),title('Respuesta en frecuencia del filtro Chebyshev 2');
% 
% SSC_filtrada_2_cheby2=filter(bc2_2,ac2_2,SSC);
% figure('Name','Señal vs Señal filtrada paso bajo con Chebyshev 2'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,SSC_filtrada_2_cheby2),
% title('Señal filtrada paso bajo con filtro Chebyshev 2 de orden 2');

% %Probamos con un filtro cheby 2 de orden 12
% [bc2_12,ac2_12] = cheby2(12,80,Wcb); % 80 son los decibelios en la banda de rechazo
% 
% figure('Name','Respuesta del sistema en frecuencia filtro Chebyshev 2'),
% freqz(bc2_12,ac2_12),title('Respuesta en frecuencia del filtro Chebyshev 2');
% 
% SSC_filtrada_12_cheby2=filter(bc2_12,ac2_12,SSC);
% figure('Name','Señal vs Señal filtrada paso bajo con Chebyshev 12'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,SSC_filtrada_12_cheby2),
% title('Señal filtrada paso bajo con filtro Chebyshev 2 de orden 12');
% 
% %Mostramos la diferencia entre los filtros Chebyshev 2 de diferente orden 
% figure('Name','Señal vs Señal filtrada paso bajo con filtro Chebyshev 2'),
% subplot(411),plot(tm,SSC),title('Señal original')
% subplot(412),plot(tm,SSC_filtrada_cheby2),title('Señal filtrada paso bajo con filtro elíptico orden 6'),
% subplot(413),plot(tm,SSC_filtrada_2_cheby2),
% title('Señal filtrada paso bajo con filtro Chebyshev 2 orden 2'),
% subplot(414),plot(tm,SSC_filtrada_12_cheby2),
% title('Señal filtrada paso bajo con filtro Chebyshev 2 orden 12');



% 
% %Probamos con el filtro elíptico de orden 6
% [be,ae] = ellip(6,5,80,Wcb);
% 
% figure('Name','Respuesta del sistema en frecuencia filtro elíptico orden 6'),
% freqz(be,ae),title('Respuesta en frecuencia del filtro Elíptico orden 6');
% 
% FPB_E = filter(be,ae,SSC); 
% figure('Name','Señal vs Señal filtrada paso bajo con filtro elíptico orden 6'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,FPB_E),
% title('Señal filtrada paso bajo con filtro elíptico orden 6');
% 
%Probamos con un filtro elíptico de orden 2
[be2,ae2] = ellip(2,5,80,Wcb);

figure('Name','Respuesta del sistema en frecuencia filtro elíptico orden 2'),
freqz(be2,ae2),title('Respuesta en frecuencia del filtro Elíptico orden 2');

FPB2_E = filtfilt(be2,ae2,SSC); 
figure('Name','Señal vs Señal filtrada paso bajo con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPB2_E),
title('Señal filtrada paso bajo con filtro elíptico orden 2');

% %Probamos con Filtro elíptico de orden 12.
% [be3,ae3] = ellip(12,5,80,Wcb);
% figure();
% freqz(be3,ae3), title('Respuesta en frecuencia del filtro Elíptico');
% 
% FPB3_E = filter(be3,ae3,SSC); 
% figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico de orden 12'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,FPB3_E),
% title('Señal filtrada paso bajo con filtro elíptico de orden 12');
% 
% %Mostramos la diferencia entre los filtros elípticos 2 del mismo orden 
% figure('Name','Señal vs Señal filtrada paso bajo con filtro elíptico'),
% subplot(411),plot(tm,SSC),title('Señal original')
% subplot(412),plot(tm,FPB_E),title('Señal filtrada paso bajo con filtro elíptico orden 6'),
% subplot(413),plot(tm,FPB2_E),title('Señal filtrada paso bajo con filtro elíptico orden 2');
% subplot(414),plot(tm,FPB3_E),title('Señal filtrada paso bajo con filtro elíptico orden 12');


%% 3.2. Filtro paso alto para eliminar la componente de continua y para eliminar las
% bajas frecuencias 


% %Probamos con el filtro paso alto implementado en la presentación Prezi
% a1 = [1 -1];
% b1 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1];
% SSC_filtrada_paso_alto = filter(b1,a1,SSC);
% SSC_filtrada_paso_alto = SSC_filtrada_paso_alto/32;
% SSC_filtrada_paso_alto = SSC-SSC_filtrada_paso_alto;
% 
% figure, subplot(2,1,1), plot(tm,SSC,'b'),title('Señal ECG original'),
% subplot(2,1,2),plot(tm,SSC_filtrada_paso_alto,'r');
% title('Señal ECG con filtro paso alto presentación de frec de corte de 5Hz');

%Calculamos la frecuencia de corte normalizada de los filtros
Fca = 5; % Frecuencia de corte (en Hz)
Wca = Fca/(Fs/2); % Frecuencia de corte normalizada (en radianes)



%Probamos con filtro Chebyshev 2 de orden 6.
% [bc2_PA,ac2_PA] = cheby2(6,80,Wca,'high'); % 80 son los decibelos en la banda de rechazo
% figure; 
% freqz(bc2_PA,ac2_PA), title('Respuesta en frecuencia del filtro PA Chebysehv 2 orden 6');
% title('Filtro Chebyshev 2');
% FPA_C2 = filter(bc2_PA,ac2_PA,SSC); 
% 
% figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,FPA_C2),
% title('Señal filtrada paso alto con filtro Chebyshev 2 de orden 6');
% 
% %Probamos con un filtro cheby2 de orden 2
% [bc2_PA_2,ac2_PA_2] = cheby2(2,80,Wca,'high'); % 80 son los decibelos en la banda de rechazo
% figure; 
% freqz(bc2_PA_2,ac2_PA_2), title('Respuesta en frecuencia del filtro PA Chebysehv 2 orden 2');
% title('Filtro Chebyshev 2');
% FPA_C2_2 = filter(bc2_PA_2,ac2_PA_2,SSC); 
% 
% figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,FPA_C2_2),
% title('Señal filtrada paso alto con filtro Chebyshev 2 de orden 2');
% 
% %Probamos con un filtro cheby2 de orden 12
% [bc2_PA_12,ac2_PA_12] = cheby2(12,80,Wca,'high'); % 80 son los decibelos en la banda de rechazo
% figure; 
% freqz(bc2_PA_12,ac2_PA_12), title('Respuesta en frecuencia del filtro PA Chebysehv 2 orden 12');
% title('Filtro Chebyshev 2');
% FPA_C2_12 = filter(bc2_PA_12,ac2_PA_12,SSC); 
% 
% figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,FPA_C2_12),
% title('Señal filtrada paso alto con filtro Chebyshev 2 de orden 12');
% 
% %Mostramos la diferencia entre los filtros Chebyshev 2 de diferente orden 
% figure('Name','Señal vs Señal filtrada paso alto con filtro Chebyshev 2'),
% subplot(411),plot(tm,SSC),title('Señal original')
% subplot(412),plot(tm,FPA_C2),title('Señal filtrada paso alto con filtro Chebyshev 2 orden 6'),
% subplot(413),plot(tm,FPA_C2_2),
% title('Señal filtrada paso alto con filtro Chebyshev 2 orden 12'),
% subplot(414),plot(tm,FPA_C2_12),
% title('Señal filtrada paso alto con filtro Chebyshev 2 orden 2');



% 
% %Probamos con Filtro elíptico de orden 6.
% [be_PA,ae_PA] = ellip(6,5,80,Wca,'high');
% figure();
% freqz(be_PA,ae_PA), title('Respuesta en frecuencia del filtro Elíptico');
% FPA_E = filter(be_PA,ae_PA,SSC); 
% 
% figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,FPA_E),
% title('Señal filtrada paso alto con filtro elíptico de orden 6');

%Probamos con un filtro elíptico de orden 2
[be_PA2,ae_PA2] = ellip(2,5,80,Wca,'high');
figure();
freqz(be_PA2,ae_PA2), title('Respuesta en frecuencia del filtro Elíptico');
FPA_E2 = filter(be_PA2,ae_PA2,SSC); 

figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
subplot(212),plot(tm,FPA_E2),
title('Señal filtrada paso alto con filtro elíptico de orden 2');
% 
% %Probamos con Filtro elíptico de orden 12.
% [be_PA3,ae_PA3] = ellip(12,5,80,Wca,'high');
% figure();
% freqz(be_PA3,ae_PA3), title('Respuesta en frecuencia del filtro Elíptico');
% FPA_E3 = filter(be_PA3,ae_PA3,SSC); 
% 
% figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
% subplot(211),plot(tm,SSC),title('Señal original sin filtrar'),
% subplot(212),plot(tm,FPA_E3),
% title('Señal filtrada paso alto con filtro elíptico de orden 12');

% %Mostramos la diferencia entre los filtros elípticos 2 de diferente orden 
% figure('Name','Señal vs Señal filtrada paso alto con filtro elíptico'),
% subplot(411),plot(tm,SSC),title('Señal original')
% subplot(412),plot(tm,FPA_E),title('Señal filtrada paso alto con filtro elíptico orden 6'),
% subplot(413),plot(tm,FPA_E2),title('Señal filtrada paso alto con filtro elíptico orden 2'),
% % subplot(414),plot(tm,FPA_E3),title('Señal filtrada paso alto con filtro elíptico orden 12');



%Señal final con los filtros seleccionados
%Elegimos el FP cheby Orden 6
soriginal_filtrada_PB=FPB2_E;
figure, plot(tm,soriginal_filtrada_PB);
%Filtramos esa señal con el filtro elíptico PA de orden 2
sfinal=filtfilt(be_PA2,ae_PA2,soriginal_filtrada_PB); 
figure, plot(tm,sfinal);
% % 
%Para la presentación de Prezi añadimos esto
sfinal = sfinal/32;
sfinal = soriginal_filtrada_PB-sfinal;
figure, plot(tm,sfinal);

%% Actividad 4: Derivación
% Derivacion codigo de la presentacion
a3 = 8;
b3 = [2 1 0 -1 -2];
s_derivada = (filter(b3,a3,sfinal).^2);
figure, subplot(2,1,1), plot(tm,sfinal,'b'), title('Señal filtrada');
subplot(2,1,2), plot(tm,s_derivada,'r'), title('Señal filtrada derivada');

%Normalizamos la señal
mins=min(s_derivada);
s_derivada=s_derivada-mins;

maxs=max(s_derivada(50:end-50));

s_derivada=s_derivada/maxs;


% %Derivamos sfinal con diff
% D_sfinal=diff(sfinal);
% figure(),plot(D_sfinal);
% D_sfinal_cuadrado=(D_sfinal);
% figure(),plot(D_sfinal_cuadrado);

%% Actividad 5: Detectar los QRS y marcarlos
% Integracion codigo de la presentacion
vit = 0.150*Fs; % Ventana integradora
vit = round(vit);
b4 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
a4 = vit;
s_integrada = filter(b4,a4,s_derivada);

%Normalizamos la señal
mins=min(s_integrada);
s_integrada=s_integrada-mins;

maxs=max(s_integrada);
s_integrada=s_integrada/maxs;

figure, subplot(2,1,1), plot(tm,s_derivada,'b'), title('Señal derivada');
subplot(2,1,2), plot(tm,s_integrada,'r'), title('Señal integrada');


%Marcamos el umbral 
u=0.3;
ene=length(s_integrada);
detector_umbral=s_integrada;
for i=1:ene
    if detector_umbral(i)<u %Si es menor que el umbral igualamos a 0
        detector_umbral(i)=0;
    end
    if detector_umbral(i)>=u %Si es mayor que el umbral igualamos a 1
        detector_umbral(i)=1;
    end
end
figure(),subplot(2,1,1), plot(s_integrada),
title('Señal integrada'),
subplot(2,1,2), plot(detector_umbral)
title('Señal con umbral');


%% Actividad 5.b: Contador de latidos
z=0;
cont=1;
ind=[];
for i=2:(ene-1) %desde segunda muestra hasta penultima muestra
    j=i;
    if detector_umbral(i)>detector_umbral(i-1) && j>=2
        while z<5 && j>=2
            m=detector_umbral(j)-detector_umbral(j-1);
            if m<.001 %Si los dos puntos distintos tinen una diferencia de valor de 0.001
                z=z+1; %Añado uno al contador z
            else
                z=0; %Sino contador a 0
            end
            j=j-1;
        end
    end
end


%% Act 6. Calcular la performance del algoritmo

%Leemos las anotaciones para saber dónde hay latidos
%[ann,anntype] = rdann('nsrdb/16265','atr',[],3840,0,'N');
%[ann,anntype] = rdann('nsrdb/16272','atr',[],3840,0,'N');
%[ann,anntype] = rdann('nsrdb/16273','atr',[],3840,0,'N');
%[ann,anntype] = rdann('nsrdb/16420','atr',[],3840,0,'N');
%[ann,anntype] = rdann('nsrdb/16483','atr',[],3840,0,'N');
%[ann,anntype] = rdann('nsrdb/16539','atr',[],3840,0,'N');
[ann,anntype] = rdann('nsrdb/16773','atr',[],3840,0,'N');
%[ann,anntype] = rdann('nsrdb/16786','atr',[],3840,0,'N');
%[ann,anntype] = rdann('nsrdb/16795','atr',[],3840,0,'N');
% [ann,anntype] = rdann('nsrdb/17052','atr',[],3840,0,'N');

%Cuento los latidos detectados por mi umbral
contador_latidos=0;
detector_umbral=detector_umbral';
detector=[]; %Para guardar el principio y el final
for i=2:numel(detector_umbral)-1
    if detector_umbral(i)==1 && detector_umbral(i-1)==0
        contador_latidos=contador_latidos+1;
        detector=[detector i];
    elseif detector_umbral(i)==1 && detector_umbral(i+1)==0
        detector=[detector i];
    end
end

%Pintamos la señal con principio y final
figure(),plot(tm,sfinal), hold on,
for i=1:numel(detector)
    plot(tm(detector(i)),sfinal(detector(i)),'x- ')
end

%Contamos verdaderos positivos y falsos positivos
verdadero_positivo=0;
falso_negativo=0;
j=1;
for i=1:numel(ann)
    if detector_umbral(ann(j))==1 || detector_umbral(ann(j)+10)==1
        verdadero_positivo=verdadero_positivo+1;
    else
        falso_negativo=falso_negativo+1;
    end
    j=j+1;
end

%Pintamos las anotaciones sobre la señal umbral para ver que esto es cierto
figure("Name","Anotaciones sobre señal umbral");
plot(tm,detector_umbral);
hold on;
grid on;
plot(tm(ann),detector_umbral(ann),'ro','MarkerSize',4);

% Para pintar el tipo de anotación al graficar
text(tm(ann),detector_umbral(ann),anntype);

%Calculamos sensibilidad
sensibilidad= verdadero_positivo/(verdadero_positivo+falso_negativo);
