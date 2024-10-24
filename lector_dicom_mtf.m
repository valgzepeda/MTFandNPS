%% Lectura de archivos DICOM

clear;
clc;
close all;

base_path = 'C:\Users\valgz\OneDrive\Documentos\MAESTRÍA\Tesis_maestría\RTPAbdom_1700\';
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

%% Función Seleccion de cortes para medir MTF 

range_mtf = 54:60;  %Rango que cortes que usaremos para medir la MTF

%Coordenadas de las esquinas de las ROIs
% 294 110 : coordenadas pixel sup izq
% 343 159 : coordenadas pixel inf der

y_range = 294:343; % Rango de Y para delimitar la ROI
x_range = 110:159; % Rango de X para delimitar la ROI 

function [K_mtf, K_rois] = CortesYROIs(K,num_files,range_mtf, x_range,y_range)

    %Inicializamos la celda donde guardaremos los cortes que nos importan
    K_mtf = cell(1, num_files);
        %Llenamos la celda con los cortes de range_mtf
        for i = 1:num_files
            K_mtf{1,i} = K{1,i}(:,:,range_mtf);
        end
    
    %Inicializamos la celda donde guardaremos las ROIs de los cortes 
    K_rois = cell(1, num_files);
        %Llenamos K_rois con las ROIs de los cortes en K_mtf
        for i = 1:num_files
        K_rois{1,i} = K_mtf{1,i}(x_range, y_range, :);
        end
end

    [K1_mtf,K1_rois] = CortesYROIs(K1,num_files,range_mtf,x_range,y_range);
    [K2_mtf,K2_rois] = CortesYROIs(K2,num_files,range_mtf,x_range,y_range);
    

%% Figuritas ilustrativas 

figure(1)
subplot(1,2,1), imshow(K1_rois{1,1}(:,:,1), []);
title('ROI de K1\_rois');
axis image
colorbar

subplot(1,2,2), imshow(K1{1,1}(:,:,range_mtf(1)), []);
title('Corte original de K1\_mtf');
axis image
colorbar

%% Función Concatenar ROIs y ajustar escala de HU 

function [K_array,K_array_correc, avg] = ConcatEscala(K_rois, num_files, range_mtf)
    
    % Inicializar una matriz 3D para almacenar las ROIs concatenadas
    K_array = zeros(size(K_rois{1,1},1),size(K_rois{1,1},2), length(range_mtf)*num_files);

        % Rellenar la matriz K_array con las ROIs concatenadas
        for i = 1:num_files
            start_slice = (i-1)*length(range_mtf) + 1;
            end_slice = i*length(range_mtf);
            K_array(:,:,start_slice:end_slice) = K_rois{1,i};
        end

    % Hay que hacer un ajuste en la escala de HU   

    %0028, 1052: RescaleIntercept = -8192
    %0028, 1053: RescaleSlope = 1

    RescaleSlope = 1; 
    RescaleIntercept = -8192; 
    %Corrección UH = RescaleSlope*(pixel value) + RescaleIntercept
    K_array_correc = RescaleSlope.*K_array+ RescaleIntercept;
    avg = mean(K_array_correc,3);

end

[K1_array,K1_array_correc,avg1] = ConcatEscala(K1_rois,num_files,range_mtf);
[K2_array,K2_array_correc,avg2] = ConcatEscala(K2_rois,num_files,range_mtf);

%% Figuritas ilustrativas parte2

figure(2) 
subplot(2,2,1), imshow(K1_rois{1,1}(:,:,3),[]);
axis image 
title('corte de una ROI de K1\_MTF\_rois'); 
colorbar
subplot(2,2,2), imshow(K2_rois{1,1}(:,:,3),[]);
axis image
title('corte de una ROI del K2\_MTF'); 
colorbar

subplot(2,2,3), imshow(avg1,[]);
axis image
title('promedio de las ROIs para el kernel 1 (Sa36) corregido');
colorbar
subplot(2,2,4), imshow(avg2,[]);
axis image
title('promedio de las ROIs para el kernel 2 (Hn 44) corregido');
colorbar

%% Obtención de centroides para MTF

% %Guardamos la imagen promedidada en una variable
% grayImage = avg1;
% %Binarizamos la imagen
% grayImage_bin = imbinarize(rescale(grayImage));
% %Los datos obtenidos de regionprops los guardamos
% measurements = regionprops(grayImage_bin, 'centroid');
% %Guardamos las coordenadas del centroide 
% centerOfMass = measurements.Centroid; %Coordenadas centroide
% 

%% Función Obtener centroides 

function centerOfMass = Centroides(avg)
    
    grayImage_bin = imbinarize(rescale(avg)); 
    measurements = regionprops(grayImage_bin, 'centroid');
    %Guardamos las coordenadas del centroide 
    centerOfMass = measurements.Centroid; %Coordenadas centroide

end 

centerOfMass_K1 = Centroides(avg1); 
centerOfMass_K2 = Centroides(avg2); 

% Figurita: Muestra el centroide sobre la imagen de la ROI
% figure(3)
% title('Ubicación del centroide en la imagen')
% imshow(grayImage_bin, []); hold on;
% axis image
% colorbar
% hold on
% for x = 1: numel(measurements)
%      plot(measurements(x).Centroid(1), measurements(x).Centroid(2),'ro');
% end
% hold off

%% Función Perfiles radiales en coord polares

function [radial_profiles,max_radius] = RadialProfiles(avg,centerOfMass)
    xCentroide = centerOfMass(1);
    yCentroide = centerOfMass(2);

    % Definimos los ángulos desde 0° a 359° con un paso de 1°
    angles = 0:1:359;  % Ángulos en grados

    % Convertimos los ángulos a radianes
    angles_rad = deg2rad(angles);
    
    % Definimos el vector r (distancia radial) desde el centroide
    r = 0:49;  % El radio máximo es 49, dado que la imagen es de 50x50

    % Inicializamos la matriz para almacenar los perfiles radiales
    radial_profiles = zeros(length(angles_rad), length(r));

    % Recorremos cada ángulo para obtener los perfiles
    for i = 1:length(angles_rad)
        theta = angles_rad(i);  % Ángulo actual en radianes

        % Convertimos las coordenadas polares a cartesianas
        [x, y] = pol2cart(theta, r);

        % hacemos el desplazamiento pertinente para que las coordenadas sean
        % respecto al centroide 
        x = x +xCentroide;
        y = y + yCentroide;

        % Usamos improfile para obtener el perfil a lo largo de la línea
        profile = improfile(avg, [xCentroide x(end)], [yCentroide y(end)], length(r),'bicubic');

        % Guardamos el perfil en la matriz radial_profiles
        radial_profiles(i, :) = profile;
    end
    
    % Definimos la variable max_radius como el tamaño del vector r
    max_radius = length(r);

end

[radial_profiles_K1,max_radius_K1] = RadialProfiles(avg1, centerOfMass_K1);
[radial_profiles_K2,max_radius_K2] = RadialProfiles(avg2, centerOfMass_K2);

% Graficamos los perfiles radiales
    figure(3)
    subplot(1,2,1)
    hold on
    % Recorremos cada perfil radial y lo graficamos
    for i = 1:size(radial_profiles_K1, 1)
        plot(0:max_radius_K1-1, radial_profiles_K1(i, :));
    end

    xlabel('Distancia radial (píxeles)');
    ylabel('Intensidad');
    title('Perfiles radiales desde 0° a 360° para el Kernel 1');
    xlim([0 max_radius_K1-1]);
    grid on;
    hold off
    subplot(1,2,2)
    hold on
    % Recorremos cada perfil radial y lo graficamos
    for i = 1:size(radial_profiles_K2, 1)
        plot(0:max_radius_K2-1, radial_profiles_K2(i, :));
    end

    xlabel('Distancia radial (píxeles)');
    ylabel('Intensidad');
    title('Perfiles radiales desde 0° a 360° para el Kernel 2');
    xlim([0 max_radius_K2-1]);
    grid on;
    hold off

%% Función Cálculo LSF 

function lsf = CalcularLSF(radial_profiles)
    %Inicializamos la matriz donde guardaremos la lsf 
    lsf = zeros(size(radial_profiles,1),size(radial_profiles,2)-1);

    for i = 1:size(radial_profiles,1)
        lsf(i,:) = abs(diff(radial_profiles(i,:)));
   
        nanIndices = isnan(lsf(i, :));
    
        % Reemplazar NaNs por 0
        lsf(i, nanIndices) = 0; 
   
        lsf(i,:) = lsf(i,:)/sum(lsf(i,:));  % Normalización
    end
end

lsf1 = CalcularLSF(radial_profiles_K1); 
lsf2 = CalcularLSF(radial_profiles_K2);

%Graficamos 
figure(4)
subplot(1,2,1)
hold on
% Recorremos cada LSF y lo graficamos
for i = 1:size(lsf1, 1)
    plot(lsf1(i, :));
end
xlabel('Distancia radial (píxeles)');
ylabel('LSF (Intensidad normalizada)');
title('LSF para cada ángulo de 0° a 360° Kernel 1');
xlim([0 size(lsf1, 2)-1]);
grid on;
hold off

subplot(1,2,2)
hold on
% Recorremos cada LSF y lo graficamos
for i = 1:size(lsf2, 1)
    plot(lsf2(i, :));
end
xlabel('Distancia radial (píxeles)');
ylabel('LSF (Intensidad normalizada)');
title('LSF para cada ángulo de 0° a 360° Kernel 2');
xlim([0 size(lsf2, 2)-1]);
grid on;
hold off

%% Función Cálculo de MTF 

function MTF = CalculoMTF(lsf)
    % Define el tamaño del padding
    padding_size = 1000;

    % Crea una nueva matriz con el tamaño extendido
    padded_lsf = zeros(size(lsf, 1), size(lsf, 2) + padding_size);
    MTF = zeros(size(lsf,1),size(lsf, 2) + padding_size);

    % Recorre cada fila de la matriz lsf y agrega el padding
    for i = 1:size(lsf, 1)
        % Coloca el perfil original en las primeras columnas
        padded_lsf(i, 1:size(lsf, 2)) = lsf(i, :);
        MTF(i, :) = abs(fftshift((fft(padded_lsf(i,:)))));
    end
end 

MTF1 = CalculoMTF(lsf1);
MTF2 = CalculoMTF(lsf2); 

% Graficamos la MTF para cada perfil y cada kernel 
figure
subplot(1,2,1)
hold on
for i = 1:size(MTF1, 1)
    plot(MTF1(i, :));
end
xlabel('Pixeles');
ylabel('MTF (Amplitud normalizada)');
title('MTF para cada ángulo de 0° a 360° Kernel 1');
xlim([0 size(MTF1, 2)-1]);
grid on
hold off

subplot(1,2,2)
hold on
for i = 1:size(MTF2, 1)
    plot(MTF2(i, :));
end
xlabel('Pixeles');
ylabel('MTF (Amplitud normalizada)');
title('MTF para cada ángulo de 0° a 360° Kernel 2');
xlim([0 size(MTF2, 2)-1]);
grid on
hold off

%% Conversión de ejes 

% 0028,0030  Pixel Spacing: 0.445390625\0.445390625 

% Resolution:  2.2452 pixels per mm
%Voxel size: 0.4454x0.4454x3 mm^3

PixelSize = 0.4454; %mm; 

%Necesitamos el valor del límite de Nyquist
% w < = 1 / 2u_b
% Es decir: 
% delta_f = 1 / 2*Delta_x
%En otras palabras: 
% delta_f = 1 / ((tamaño del pixel) * (número de pixeles) = 1 / 2u_b

ImageSize = size(MTF1,2);
a = PixelSize * ImageSize;

% Límite de Nyquist en ciclos por mm
w = 1 / a;

% El soporte compacto de la frecuencia espacial es W 
%Para definir el eje de frecuencias definimos 
frecuencia_espacial = linspace(-w*(1049/2),w*(1049/2),ImageSize); 
figure
subplot(1,2,1)
hold on 
for i =1:size(MTF1,1)
    plot(frecuencia_espacial,MTF1(i,:))
end
%xlim([0 w])
xlabel('Frecuencia espacial [1/mm]');
ylabel('MTF (Amplitud normalizada)');
title('MTF Kernel 1 ');
hold off

subplot(1,2,2)
hold on 
for i =1:size(MTF2,1)
    plot(frecuencia_espacial,MTF1(i,:))
end
%xlim([0 w])
xlabel('Frecuencia espacial [1/mm]');
ylabel('MTF (Amplitud normalizada)');
title('MTF kernel 2 ');
hold off

%% Función MTF frecuencias positivas 
function [MTF_positiva,frecuencia_espacial_positiva] = FrecsPositiv(MTF,frecuencia_espacial)
% Encuentra el índice correspondiente a la frecuencia 0
[~, idx_cero] = min(abs(frecuencia_espacial)); % Encuentra el índice más cercano a 0

% Extrae las frecuencias espaciales positivas
MTF_positiva = MTF(:, idx_cero:end);

%Extrae las frecuencias espaciales correspondientes a los valores positivos
frecuencia_espacial_positiva = frecuencia_espacial(idx_cero:end);
end

[MTF_positiva1,frecuencia_espacial_positiva] = FrecsPositiv(MTF1,frecuencia_espacial);
[MTF_positiva2,frecuencia_espacial_positiva] = FrecsPositiv(MTF1,frecuencia_espacial);

% Graficar las frecuencias espaciales positivas
figure
subplot(1,2,1)
hold on 
for i = 1:size(MTF_positiva1, 1)
    plot(frecuencia_espacial_positiva, MTF_positiva1(i, :))
end

xlabel('Frecuencia espacial [1/mm]');
ylabel('MTF (Amplitud normalizada)');
title('MTF para cada ángulo de 0° a 360° (frecuencias positivas) Kernel 1');
hold off
subplot(1,2,2)
hold on 
for i = 1:size(MTF_positiva2, 1)
    plot(frecuencia_espacial_positiva, MTF_positiva2(i, :))
end

xlabel('Frecuencia espacial [1/mm]');
ylabel('MTF (Amplitud normalizada)');
title('MTF para cada ángulo de 0° a 360° (frecuencias positivas) Kernel 2');
hold off

%% Distribución angular de la MTF 

[num_angles,profile_size] = size(MTF_positiva1);
size_disp = 150;

%Creamos el meshgrid auxiliar. Definimos un conjunto de puntos para
%x y para y 
x_prima = linspace(-profile_size,profile_size,size_disp);
y_prima = linspace(-profile_size,profile_size,size_disp);
[x_mesh,y_mesh] = meshgrid(x_prima,y_prima);

%Ahora hay que colocar la MTF en estos nuevos puntos. Para ello hay que 
%crear un conjunto de puntos del vector R y de angulos para mapear 

%Conjunto de ángulos (de 0 a 2pi con 360 pasos) 
valores_angulos = linspace(0,2*pi,num_angles);

%R va de 1 hasta el tamaño del perfil porque en la MTF positiva ya todos los
%perfiles están centrados en el mismo punto
R_auxiliar = linspace(1,profile_size,profile_size); 

%Definimos un vector radial R y un rango de angulos Theta en términos de
%los puntos del meshgrid (coord polares en este nuevo plano cartesiano)

R = sqrt(x_mesh.^2+y_mesh.^2);
Theta = atan2(x_mesh, y_mesh); 
    
%Interpolamos los valores de R y Theta para que nos den valores angulares
%entre los perfiles que no tenemos, para crear un set de índices radial y
%angular completos

%interp1(x,v,xq) devuelve valores interpolados de una función 1D en puntos 
% de consulta específicos utilizando la interpolación lineal. 
% El vector x contiene puntos de muestra y v contiene los valores 
% correspondientes, v (x). El vector xq contiene las coordenadas de los 
% puntos de consulta.

R_index = interp1(R_auxiliar,1:profile_size,R,'linear',NaN);
Theta_index = interp1(valores_angulos,1:num_angles,Theta,'linear',NaN);

%Delimitamos apropiadamente los valores de R_index y Theta_index

R_index = max(min(R_index,profile_size),1);
Theta_index = max(min(Theta_index,num_angles),1);

%Inicializamos la matriz donde guardaremos todo 

ang_MTF = NaN(size(x_mesh));

for i = 1:size_disp

    for j = 1:size_disp
        if ~isnan(R_index(i,j)) && ~isnan(Theta_index(i,j))
        ang_MTF(i,j) = interp2(MTF_positiva1,R_index(i,j),Theta_index(i,j),'linear');
        end
    end
end

 ang_MTF(isnan(ang_MTF)) = 0;

 
%% Gráfica distribución angular de la MTF en 2D
figure;
imagesc(ang_MTF);
axis image; 
colormap('jet'); 
colorbar;  
title('Distribución angular de la MTF en 2D');
xlabel('X');
ylabel('Y');

%% Cambio de ejes a frecuencia espacial 
 
PixelSize = 0.4454; %mm; 
%Necesitamos el valor del límite de Nyquist
% w < = 1 / 2u_b

ImageSize_positivo = size(MTF_positiva1,2);
a = ImageSize_positivo * PixelSize;

% Límite de Nyquist en ciclos por mm
w = 1 / a;
w = w*(1049/2);
frec_es_pos = linspace(-w,w,ImageSize_positivo); 

% Graficar la distribución angular de la MTF en 2D con ejes en frecuencia espacial
figure;
imagesc(frec_es_pos, frec_es_pos, ang_MTF);  % Mapeo de las frecuencias espaciales
axis image;
colormap('jet');
colorbar;
title('Distribución angular de la MTF en 2D');
xlabel('Frecuencia espacial [1/mm] (Eje X)');
ylabel('Frecuencia espacial [1/mm] (Eje Y)');
