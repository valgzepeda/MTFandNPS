%% Lectura de archivos DICOM

clear;
clc;
close all;

base_path = 'C:\Users\valgz\OneDrive\Documentos\MAESTRÍA\Tesis_maestría\RTPAbdom_1700\';
file_numbers1 = 26632:26641; % Números de archivo
num_files = length(file_numbers1);
file_numbers2 = 26642:26651;

% Inicializar el array para cada Kernel 

K1 = cell(1, num_files);
K2 = cell(1, num_files);

for i = 1:num_files
    % Leer volumenes DICOM Kernel 1 y Kernel 2
    K1{1,i} = squeeze(dicomreadVolume([base_path 'Abdomen3_' num2str(file_numbers1(i)) '\']));
    K2{1,i} = squeeze(dicomreadVolume([base_path 'Abdomen3_' num2str(file_numbers2(i)) '\']));
end


%% Selección de cortes a usar para medir NPS 

% range = 17:28; %Rango que cortes que usaremos para medir el ruido 
% num_slices = length(range);
% 
% % Inicializamos el espacio para las matrices de cortes
% K1_noise = cell(1, num_files);
% K2_noise = cell(1, num_files);
% 
% for i = 1:num_files
%     K1_noise{1,i} = K1{1,i}(:,:,range);
%     K2_noise{1,i} = K2{1,i}(:,:,range);
% end

%% Seleccionar el corte central NPS 
% central_slice = round(size(range,2)/2);
% 
% %Inicializamos arrays
% K1_central = cell(1, num_files);
% K2_central = cell(1, num_files);
% 
% for i = 1:num_files
%     K1_central{1,i} = K1_noise{1,i}(:,:,central_slice);
%     K2_central{1,i} = K2_noise{1,i}(:,:,central_slice);
% end
%% Función Selección de cortes 

%Rango de cortes a usar 
range = 17:28;

function [K_noise] = SeleccionarCortes(K,range)
    %Definimos el número de archivos que habrá en la celda 
    num_files = size(K,2);
    % Inicializamos el espacio para las matrices de cortes
    K_noise = cell(1, num_files);

    % Selección de los cortes para cada archivo llenando K_noise
    for i = 1:num_files
        K_noise{1,i} = K{1,i}(:,:,range);
    end

    %Inicializamos la celda 
    % K_central = cell(1, num_files);
    % 
    % for i = 1:num_files
    % K_central{1,i} = K_noise{1,i}(:,:,central_slice);
    % end
end

[K1_noise] = SeleccionarCortes(K1,range); 
[K2_noise] = SeleccionarCortes(K2,range); 

%% Seleccionar el corte central NPS 

% %Inicializamos arrays
% K1_central = cell(1, num_files);
% K2_central = cell(1, num_files);
% 
% for i = 1:num_files
%     K1_central{1,i} = K1_noise{1,i}(:,:,central_slice);
%     K2_central{1,i} = K2_noise{1,i}(:,:,central_slice);
% end

%% Función Covertir a Array, corrección de escala y promedio 

function [K_array, K_array_correc, avg] = ArrayEscalaAvg(K,K_noise)

range = 17:28;
num_files = size(K_noise,2); %10
num_slices = length(range); %12

K_array = zeros(size(K{1,1},1), size(K{1,1},2), num_slices*num_files);

    for i = 1:num_files

        start_slice = (i-1)*num_slices + 1;
        end_slice = i*num_slices;

        K_array(:,:,start_slice:end_slice) = K_noise{1,i};

    end

RescaleSlope = 1; 
RescaleIntercept = -8192; 

K_array_correc = RescaleSlope.*K_array + RescaleIntercept; 
avg = mean(K_array_correc,3);

end

[K1_array, K1_array_correc, avg1] = ArrayEscalaAvg(K1,K1_noise); 
[K2_array, K2_array_correc, avg2] = ArrayEscalaAvg(K2,K2_noise); 

%% Convertir a arrays 3D
%Inicializamos arrays del tamaño original X,Y y 
%size(K1_central,2) = num_files
% 
%  [dimX,dimY,dimZ] = size(K1{1});
% 
% K1_array = zeros(dimX, dimY, num_files);
% K2_array = zeros(dimX, dimY, num_files);
% 
% for i = 1:num_files
%     K1_array(:,:,i) = K1_central{i};
%     K2_array(:,:,i) = K2_central{i};
% end
% 
% %RescaleSlope* pixel value+rescaleIntercept
% 
% RescaleSlope = 1; 
% RescaleIntercept = -8192; 
% 
% K1_array_correc = RescaleSlope.*K1_array + RescaleIntercept; 
% K2_array_correc = RescaleSlope.*K2_array + RescaleIntercept; 
% 
% avgcorrec1 = mean(K1_array_correc,3);
% avgcorrec2 = mean(K2_array_correc,3);

%% Figuritas 

figure(1) 

subplot(2,2,1), imshow(K1_array(:,:,3),[]);
axis image
title('corte central del Kernel 1 '); 
colorbar
subplot(2,2,2), imshow(K2_array(:,:,3),[]);
axis image
title('corte central del Kernel 2 '); 
colorbar

subplot(2,2,3), imshow(avg1,[]);
axis image
title('promedio de los slices centrales para el kernel 1 (Sa36) corregido');
colorbar
subplot(2,2,4), imshow(avg2,[]);
axis image
title('promedio de los slices centrales para el kernel 2 (Hn 44) corregido');
colorbar

%% Función que da el ensamble de ruido (Corte - avg) 

function resta = Ensamble(K_array_correc,avg)

    for i = 1:size(K_array_correc,3)
        %Ensamble de ruido para el kernel 
        resta(:,:,i) = K_array_correc(:,:,i)-avg;
    end

end

    resta1 = Ensamble(K1_array_correc,avg1);
    resta2 = Ensamble(K2_array_correc,avg2); 

%% Creación de ROIs y Padding 

function Padtotal = ROIsAndPadding(resta)
    ROI1 = zeros(50,50,size(resta,3));
    ROI2 = zeros(50,50,size(resta,3));
    ROI3 = zeros(50,50,size(resta,3));
    ROI4 = zeros(50,50,size(resta,3));
    ROI5 = zeros(50,50,size(resta,3));

        for i = 1:size(resta,3)
            ROI1(:,:,i) = resta(231:280,231:280,i);

            ROI2(:,:,i) = resta(231:280,321:370,i);

            ROI3(:,:,i) = resta(321:370,231:280,i);

            ROI4(:,:,i) = resta(231:280,141:190,i);

            ROI5(:,:,i) = resta(141:190,231:280,i);

            bigROI(:,:,i)=resta(148:348,148:348,i);
        end
    %Concatenamos todas las ROIs obtenidas 
    ROItotales = cat(3,ROI1,ROI2,ROI3,ROI4,ROI5);
    Padtotal = padarray(ROItotales,[500 500],0,"both");
end

  Padtotal1 = ROIsAndPadding(resta1);
  Padtotal2 = ROIsAndPadding(resta2); 

  function bigPad = bigROIandPadding(resta)
    bigROI = zeros(200,200,size(resta,3));
    for i = 1:size(resta,3)
        bigROI(:,:,i)=resta(158:357,158:357,i);
    end 
    bigPad = padarray(bigROI,[500 500],'both');
  end

    bigPad1 = bigROIandPadding(resta1);
    bigPad2 = bigROIandPadding(resta2);

%% Pruebas de ROIs

ROI1 = zeros(50,50,size(resta1,3));
for i = 1:size(resta1,3)
            ROI1(:,:,i) = resta1(231:280,231:280,i);
end
PadTROI1 = padarray(ROI1,[500 500],'both');

ROI2 = zeros(50,50,size(resta1,3));
for i = 1:size(resta1,3)
            ROI2(:,:,i) = resta1(231:280,331:380,i);
end
PadTROI2 = padarray(ROI2,[500 500],'both');

ROI3 = zeros(50,50,size(resta1,3));
for i = 1:size(resta1,3)
            ROI3(:,:,i) = resta1(331:380,231:280,i);
end
PadTROI3 = padarray(ROI3,[500 500],'both');

ROI4 = zeros(50,50,size(resta1,3));
for i = 1:size(resta1,3)
            ROI4(:,:,i) = resta1(231:280,95:144,i);
end
PadTROI4 = padarray(ROI4,[500 500],'both');

ROI5 = zeros(50,50,size(resta1,3));
for i = 1:size(resta1,3)
            ROI5(:,:,i) = resta1(95:144,231:280,i);
end
PadTROI5 = padarray(ROI5,[500 500],'both');


%% Prueba NPS para cada ROI 

Dft_prueba_1 = zeros(size(PadTROI1,1),size(PadTROI1,2));
Dft_prueba_2 = zeros(size(PadTROI2,1),size(PadTROI2,2));
Dft_prueba_3 = zeros(size(PadTROI3,1),size(PadTROI3,2));
Dft_prueba_4 = zeros(size(PadTROI4,1),size(PadTROI4,2));
Dft_prueba_5 = zeros(size(PadTROI5,1),size(PadTROI5,2));

for i = 1:size(PadTROI1,3)
Dft_prueba_1(:,:,i) = (abs(fftshift(fft2(PadTROI1(:,:,i))))).^2;
Dft_prueba_2(:,:,i) = (abs(fftshift(fft2(PadTROI2(:,:,i))))).^2;
Dft_prueba_3(:,:,i) = (abs(fftshift(fft2(PadTROI3(:,:,i))))).^2;
Dft_prueba_4(:,:,i) = (abs(fftshift(fft2(PadTROI4(:,:,i))))).^2;
Dft_prueba_5(:,:,i) = (abs(fftshift(fft2(PadTROI5(:,:,i))))).^2;

end 

PixelSize = 0.4454;
 
NPStotalPrueba1 = ((PixelSize/50)^2).*mean(Dft_prueba_1,3);
NPStotalPrueba2 = ((PixelSize/50)^2).*mean(Dft_prueba_2,3);
NPStotalPrueba3 = ((PixelSize/50)^2).*mean(Dft_prueba_3,3);
NPStotalPrueba4 = ((PixelSize/50)^2).*mean(Dft_prueba_4,3);
NPStotalPrueba5 = ((PixelSize/50)^2).*mean(Dft_prueba_5,3);

figure
subplot(2,3,1)
imagesc(NPStotalPrueba1)
axis image 
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
title('NPS de la ROI1 roja '); 
subplot(2,3,2)
imagesc(NPStotalPrueba2)
axis image 
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
title('NPS de la ROI2 verde');
subplot(2,3,3)
imagesc(NPStotalPrueba3)
axis image 
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
title('NPS de la ROI3 azul ');
subplot(2,3,4)
imagesc(NPStotalPrueba4)
axis image 
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
title('NPS de la ROI4 magenta '); 
subplot(2,3,5)
imagesc(NPStotalPrueba5)
axis image 
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
title('NPS de la ROI5 amarilla '); 


%% Cálculo de NPS

% NPS = (delta_x*delta_y/NxNyM) * SUM (abs(DFT(ROI - ROI_avg)) )
% en otras palabras:
% NPS = (delta_x*delta_y / NxNy) * AVG ((abs(DFT(ROI- ROI_avg))*0.4454^2)^2)

% Resolution:  2.2452 pixels per mm
%Voxel size: 0.4454x0.4454x3 mm^3


Dft_total_1 = zeros(size(Padtotal1,1),size(Padtotal1,2));
Dft_total_2 = zeros(size(Padtotal1,1),size(Padtotal1,2));
bigDFT_1 = zeros(size(bigPad1,1),size(bigPad1,2));
bigDFT_2 = zeros(size(bigPad1,1),size(bigPad1,2));


for i = 1:size(Padtotal1,3)
Dft_total_1(:,:,i) = (abs(fftshift(fft2(Padtotal1(:,:,i))))).^2;
Dft_total_2(:,:,i) = (abs(fftshift(fft2(Padtotal2(:,:,i))))).^2;
end

for i = 1:size(bigPad1,3)
bigDFT_1(:,:,i) = (abs(fftshift(fft2(bigPad1(:,:,i))))).^2;
bigDFT_2(:,:,i) = (abs(fftshift(fft2(bigPad2(:,:,i))))).^2;
end 

PixelSize = 0.4454;
 
NPStotal_K1 = ((0.4454/50)^2).*mean(Dft_total_1,3);
NPStotal_K2 = ((0.4454/50)^2).*mean(Dft_total_2,3);
bigNPS_K1 = ((0.4454/50)^2).*mean(bigDFT_1,3);
bigNPS_K2 = ((0.4454/50)^2).*mean(bigDFT_2,3);

figure(5)
subplot(2,2,1), imagesc(NPStotal_K1)
axis image
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
title('NPS del Kernel 1 '); 
subplot(2,2,2), imagesc(NPStotal_K2)
axis image
title('NPS total del kernel 2');
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')

subplot(2,2,3), imagesc(bigNPS_K1)
axis image
title('NPS total del kernel 1 con ROI grande');
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')

subplot(2,2,4), imagesc(bigNPS_K2)
axis image
title('NPS total del kernel 2 con ROI grande');
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
%% Cálculo de la desviación estándar de uno de los cortes 

%Desviación estándar del kernel 1 

sd1 = std2(K1_array_correc(231:280,231:280,1));
sd2 = std2(K2_array_correc(231:280,231:280,1));

% Mostramos los resultados
disp(['Sigma para Kernel 1: ', num2str(sd1)]);
disp(['Sigma para Kernel 2: ', num2str(sd2)]);

%% Conversión de ejes 

%Necesitamos el valor del límite de Nyquist
% w < = 1 / 2u_b

%sanity check: voxel pize * # pixeles = tamaño real catphan 
% tamaño real módulo catphan 150 mm 
% 0.4454 mm * 339 pixeles = 150.9 mm 


% delta_f = 1 / Delta_F 
% delta_f = 1 / 2* Delta_x = 1 / N*delta_x 
ImageSize = size(Padtotal1);

a = 1050 * PixelSize;
% Límite de Nyquist en ciclos por mm
delta_f = 1 / a;

% Mostrar el valor
disp(['El límite de Nyquist es: ', num2str(delta_f), ' 1/mm']);

[dimX, dimY] = size(NPStotal_K1);  % Dimensiones de la NPS 

%Genera vectores de frecuencias espaciales que van desde −w hasta w 
%en ambos ejes, con una longitud igual al número de píxeles en esos ejes.

freqX = linspace(-delta_f*(dimX/2), delta_f*(dimX/2), dimX);  % Frecuencias en el eje X
freqY = linspace(-delta_f*(dimY/2), delta_f*(dimY/2), dimY);  % Frecuencias en el eje Y

%Creamos el meshgrid
[FreqX, FreqY] = meshgrid(freqX, freqY); 

figure;
subplot(1,2,1)
imagesc(freqX, freqY, NPStotal_K1);
axis image;
xlabel('Frecuencia en X (1/mm)');
ylabel('Frecuencia en Y (1/mm)');
title('NPS del Kernel 1 en el Dominio de la Frecuencia');
colorbar;
subplot(1,2,2)
imagesc(freqX, freqY, NPStotal_K2); 
axis image;
xlabel('Frecuencia en X (1/mm)');
ylabel('Frecuencia en Y (1/mm)');
title('NPS del Kernel 2 en el Dominio de la Frecuencia');
colorbar;


%% Figuritas de las ROIS

imshow(avg1,[]) % Muestra la imagen
axis image % Ajusta las proporciones del eje
hold on % Mantén la imagen para graficar sobre ella

% Coordenadas de las esquinas del cuadrado
x1 = [231, 280, 280, 231, 231]; % Cierra el cuadrado al volver al primer punto
y1 = [231, 231, 280, 280, 231]; % Lo mismo aquí para la coordenada 'y'

x2 = [231, 280, 280, 231, 231];
y2 = [321, 321, 370, 370, 321]; 


x3 = [321, 370, 370, 321, 321];
y3 = [231, 231, 280, 280, 231]; % Lo mismo aquí para la coordenada 'y'

x4 = [231, 280, 280, 231, 231];
y4 = [141, 141, 190, 190, 141]; 

x5 = [141, 190, 190, 141, 141]; % Cierra el cuadrado al volver al primer punto
y5 = [231, 231, 280, 280, 231]; % Lo mismo aquí para la coordenada 'y'

x6 = [158, 357, 357, 158, 158];
y6 = [158, 158,357,357,158];

plot(x1, y1, 'r', 'LineWidth', 2)
plot(x2,y2, 'g','LineWidth',2);
plot(x3,y3,'b','LineWidth',2);
plot(x4, y4, 'm','LineWidth',2);
plot(x5,y5,'y','LineWidth',2);
plot(x6,y6,'k','LineWidth',3);
hold off 

%% Cálculo de factores de normalización 

% Sabemos que N* nps = NPS 
% N \Sigma \Sigma nps * (delta_f )^2 = varianza^2

% N = varianza^2 N_x^2 delta_f^2 / \Sigma\Sigma nps
%lo cual es aproximadamente igual a 
% N = var(ROI)^2 N^2 delta_f^2 / \Sigma\Sigma nps

sumatoria = sum(NPStotal_K1(:));

sd1 = std2(K1_array_correc(231:280,231:280,1));


Nx = size(NPStotal_K1,1);


Norma = (sd1*Nx*delta_f)^2/sumatoria;

% Mostrar el valor
disp(['El valor de la norma es: ', num2str(Norma), ]);

% También N = (delta_x*delta_y / NxNy)
Norma_prima = (delta_f*delta_f) / dimX*dimY; 

% Mostrar el valor
disp(['El valor de la norma prima es: ', num2str(Norma_prima), ]);