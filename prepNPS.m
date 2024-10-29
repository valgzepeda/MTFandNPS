%% Lectura de archivos para MTF y NPS 

clear;
clc;
close all;

base_path = 'C:\Users\valgz\OneDrive\Documentos\MAESTRIA\Tesis_maestria\RTPAbdom_1700\';
file_numbers1 = 26632:26641; % Números de archivo
num_files = length(file_numbers1);
file_numbers2 = 26642:26651;

% Inicializamos la celda para cada Kernel 
K1 = cell(1, num_files);
K2 = cell(1, num_files);

% Leer volumenes DICOM Kernel 1 y Kernel 2
for i = 1:num_files
    K1{1,i} = squeeze(dicomreadVolume([base_path 'Abdomen3_' num2str(file_numbers1(i)) '\']));
    K2{1,i} = squeeze(dicomreadVolume([base_path 'Abdomen3_' num2str(file_numbers2(i)) '\']));
end

%% Función Selección de cortes NPS 

% Rango de cortes a usar 
range = 17:28;

function [K_noise] = SeleccionarCortes(K, range)
    % Definimos el número de archivos que habrá en la celda 
    num_files = size(K,2);
    % Inicializamos el espacio para las matrices de cortes
    K_noise = cell(1, num_files);

    for i = 1:num_files
        K_noise{1,i} = K{1,i}(:,:,range);
    end

end

% Llamada a la función
[K1_noise] = SeleccionarCortes(K1, range); 
[K2_noise] = SeleccionarCortes(K2, range);

%% Función Covertir a Array, corrección de escala y promedio NPS

function [K_array, K_array_correc, avg] = ArrayEscalaAvg(K,K_noise)

range = 17:28;
num_files = size(K_noise,2); %10
num_slices = length(range); %12

RescaleSlope = 1; 
RescaleIntercept = -8192; 

K_array = zeros(size(K{1,1},1), size(K{1,1},2), num_slices*num_files);

    for i = 1:num_files

        start_slice = (i-1)*num_slices + 1;
        end_slice = i*num_slices;

        K_array(:,:,start_slice:end_slice) = K_noise{1,i};

    end 

K_array_correc = RescaleSlope.*K_array + RescaleIntercept; 
avg = mean(K_array_correc,3);

end

[K1_array, K1_array_correc, avg1] = ArrayEscalaAvg(K1,K1_noise); 
[K2_array, K2_array_correc, avg2] = ArrayEscalaAvg(K2,K2_noise); 

%% Figuritas 

figure(10) 
subplot(2,2,1), imshow(K1_array_correc(:,:,3),[]);
axis image
title('corte para el Kernel 1 '); 
colorbar
subplot(2,2,2), imshow(K2_array_correc(:,:,3),[]);
axis image
title('corte para el Kernel 2 '); 
colorbar

subplot(2,2,3), imshow(avg1,[]);
axis image
title('promedio de los slices para el kernel 1 (Sa36) corregido');
colorbar
subplot(2,2,4), imshow(avg2,[]);
axis image
title('promedio de los slices para el kernel 2 (Hn 44) corregido');
colorbar