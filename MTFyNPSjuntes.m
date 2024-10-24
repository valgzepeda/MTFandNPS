%% Lectura de archivos DICOM MTF

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

    % rescaleSlope = 1; 
    % rescaleIntercept = -8192
    % 
    % % Aplicar la corrección a cada entrada de K1 y K2
    % K1{1,i} = rescaleSlope .* K1{1,i} + rescaleIntercept;
    % K2{1,i} = rescaleSlope .* K2{1,i} + rescaleIntercept;
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
    K_rois = cell(1, num_files);

        for i = 1:num_files
            %Llenamos la celda con los cortes de range_mtf
            K_mtf{1,i} = K{1,i}(:,:,range_mtf);
            %Llenamos K_rois con las ROIs de los cortes en K_mtf
            K_rois{1,i} = K_mtf{1,i}(x_range, y_range, :);
        end
end

 [K1_mtf,K1_rois] = CortesYROIs(K1,num_files,range_mtf,x_range,y_range);
 [K2_mtf,K2_rois] = CortesYROIs(K2,num_files,range_mtf,x_range,y_range);
    

%% Figuritas ilustrativas 

figure(1)
subplot(2,2,1), imshow(K1_rois{1,1}(:,:,1), []);
title('ROI de K1\_rois');
axis image
colorbar

subplot(2,2,2), imshow(K1{1,1}(:,:,range_mtf(1)), []);
title('Corte original de K1\_mtf');
axis image
colorbar

subplot(2,2,3), imshow(K2_rois{1,1}(:,:,1), []);
title('ROI de K2\_rois');
axis image
colorbar

subplot(2,2,4), imshow(K2{1,1}(:,:,range_mtf(1)), []);
title('Corte original de K2\_mtf');
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
title('promedio de las ROIs para el kernel 1 (Sa36)');
colorbar
subplot(2,2,4), imshow(avg2,[]);
axis image
title('promedio de las ROIs para el kernel 2 (Hn 44)');
colorbar

%% Función Obtener centroides MTF

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

%% Función Perfiles radiales en coord polares MTF

function [radial_profiles,max_radius] = RadialProfiles(avg,centerOfMass)
    xCentroide = centerOfMass(1);
    yCentroide = centerOfMass(2);

    % Definimos los ángulos desde 0° a 359° con un paso de 1°
    angles = 0:1:359;  % Ángulos en grados

    % Convertimos los ángulos a radianes
    angles_rad = deg2rad(angles);
    
    % Definimos el vector r (distancia radial) desde el centroide
    r = 0:50;  % El radio máximo es 49, dado que la imagen es de 50x50

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
    grid on

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
%xlim([0 size(lsf1, 2)-1]);
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
%xlim([0 size(lsf2, 2)-1]);
grid on;
hold off

%% Función Cálculo de MTF para cada ángulo

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
figure(5)
subplot(1,2,1)
hold on
for i = 1:size(MTF1, 1)
    plot(MTF1(i, :));
end
xlabel('Pixeles');
ylabel('MTF (Amplitud normalizada)');
title('MTF para cada ángulo de 0° a 360° Kernel 1');
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
grid on
hold off


%% Conversión de ejes 

% % El soporte compacto de la frecuencia espacial es W 
% % El soporte compacto es el valor de w por el número de pixeles en la imagen
% %Para definir el eje de frecuencias definimos: 
% 
% frecuencia_espacial = linspace(-delta_f_MTF1*(PixelNumber1/2),delta_f_MTF1*(PixelNumber1/2),PixelNumber1); 
% figure(6)
% subplot(1,2,1)
% hold on 
% for i =1:size(MTF1,1)
%     plot(frecuencia_espacial,MTF1(i,:))
% end
% xlabel('Frecuencia espacial [1/mm]');
% ylabel('MTF (Amplitud normalizada)');
% title('MTF Kernel 1 ');
% grid on
% hold off
% 
% subplot(1,2,2)
% hold on 
% for i =1:size(MTF2,1)
%     plot(frecuencia_espacial,MTF1(i,:))
% end
% xlabel('Frecuencia espacial [1/mm]');
% ylabel('MTF (Amplitud normalizada)');
% title('MTF kernel 2 ');
% grid on
% hold off

%% Seleccion MTF positivas 

function MTFp = ExtraerFrecuenciasPositivas(MTF)
    % Inicializamos la matriz MTF1p con un tamaño provisional
    MTFp = zeros(size(MTF, 1), size(MTF, 2)/2);  % Aprox. la mitad del tamaño de MTF

    % Recorremos cada fila de la matriz MTF
    for i = 1:size(MTF, 1)
        % Buscamos el índice del valor más cercano a 1 en cada fila
        [~, index_cercano_a_uno] = min(abs(MTF(i,:) - 1));
        
        % Almacenamos los valores desde ese índice hasta el final
        valores_positivos = MTF(i, index_cercano_a_uno:end);
        
        % Ajustamos el tamaño de la fila en MTF1p para que coincida con los valores extraídos
        MTFp(i, 1:length(valores_positivos)) = valores_positivos;
    end
end

% Ahora aplicamos esta función para obtener MTF1p
%frecuencia_espacial = linspace(-delta_f_MTF1*(PixelNumber1/2),delta_f_MTF1*(PixelNumber1/2),PixelNumber1); 
MTF1p = ExtraerFrecuenciasPositivas(MTF1);
MTF2p = ExtraerFrecuenciasPositivas(MTF2);

figure(6)
subplot(1,2,1)
hold on 
for i = 1:size(MTF1p,1)
plot(MTF1p(i,:))
end
hold off
xlabel('pixeles')
title('MTF (frecuencias positivas) Kernel 1');
ylabel('MTF1p')
grid on 
subplot(1,2,2)
grid on 
hold on 
for i = 1:size(MTF1p,1)
plot(MTF2p(i,:))
end
hold off
xlabel('pixeles')
ylabel('MTF2p')
title('MTF (frecuencias positivas) Kernel 2');

%% Función MTF frecuencias positivas 
function [MTF_positiva,frecuencia_espacial_positiva] = FrecsPositiv(MTF)
PixelNumber1 = size(MTF,2);
PixelSize = 0.4454;
delta_f_MTF1 = 1 / (PixelNumber1*PixelSize);
frecuencia_espacial = linspace(-delta_f_MTF1*(PixelNumber1/2),delta_f_MTF1*(PixelNumber1/2),PixelNumber1);
% Encuentra el índice correspondiente a la frecuencia 0
[~,idx_cero] = min(abs(frecuencia_espacial)); % Encuentra el índice más cercano a 0

% Extrae las frecuencias espaciales positivas
MTF_positiva = MTF(:, idx_cero:end);

%Extrae las frecuencias espaciales correspondientes a los valores positivos
frecuencia_espacial_positiva = frecuencia_espacial(idx_cero:end);
end

[MTF_positiva1,frecuencia_espacial_positiva] = FrecsPositiv(MTF1);
[MTF_positiva2,frecuencia_espacial_positiva] = FrecsPositiv(MTF2);

% Graficar las frecuencias espaciales positivas
figure(7)
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
ylabel('Y')


%% Distribución 2D de la MTF 

function MTF2D = DistribucionAngular(MTF_positiva, size_disp)
[num_angles,profile_size] = size(MTF_positiva);

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

MTF2D = NaN(size(x_mesh));

for i = 1:size_disp

    for j = 1:size_disp
        if ~isnan(R_index(i,j)) && ~isnan(Theta_index(i,j))
        MTF2D(i,j) = interp2(MTF_positiva,R_index(i,j),Theta_index(i,j),'linear');
        end
    end
end

 MTF2D(isnan(MTF2D)) = 0;

end

size_disp1 = 571;
size_disp2 = 761;

MTF2D1_positiva = DistribucionAngular(MTF_positiva1,size_disp1);
MTF2D2_positiva = DistribucionAngular(MTF_positiva2,size_disp2);

% MTF2D1 = DistribucionAngular(MTF1p,size_disp1);
% MTF2D2 = DistribucionAngular(MTF2p,size_disp2);

%%

size_disp1 = 571;
size_disp2 = 761;
 MTF2D1 = polarToCartesian(MTF_positiva1,size_disp1);
 MTF2D2 = polarToCartesian(MTF_positiva2,size_disp2);

%% Gráfica distribución angular de la MTF en 2D

figure(8)
subplot(1,2,1)
imagesc(MTF2D1);
axis image; 
colormap('jet'); 
colorbar;  
title(' MTF en 2D para el Kernel 1 con polartoC ');
xlabel('X');
ylabel('Y');
subplot(1,2,2)
imagesc(MTF2D2);
axis image; 
colormap('jet'); 
colorbar;  
title('MTF en 2D para el Kernel 2 con polartoC');
xlabel('X');
ylabel('Y');


%%
figure(9)
subplot(1,2,1)
imagesc(MTF2D1_positiva);
axis image; 
colormap('jet'); 
colorbar;  
title(' MTF en 2D para el Kernel 1 con mi codigo');
xlabel('X');
ylabel('Y');
subplot(1,2,2)
imagesc(MTF2D2_positiva);
axis image; 
colormap('jet'); 
colorbar;  
title('MTF en 2D para el Kernel 2 con mi codigo ');
xlabel('X');
ylabel('Y');


%% Cálculo de frecuencias espaciales matrices MTF
% delta_f = 1 / 2*Delta_x
%En otras palabras: 
% delta_f = 1 / ((tamaño del pixel) * (número de pixeles) = 1 / 2u_b

PixelSize = 0.4454; % Este es el valor de delta_x de la CT 

PixelNumber1 = size(MTF1,2);

delta_f_MTF1 = 1 / (PixelNumber1*PixelSize);

Delta_F_MTF1 = delta_f_MTF1*PixelNumber1;

% Mostrar el valor
disp(['Frecuencia espacial para MTF1 y MTF2: ', num2str(delta_f_MTF1),' 1/mm' ]);
disp(['Soporte compacto para MTF1 y MTF2: ', num2str(Delta_F_MTF1),' mm' ]);

% Frecuencia espacial MTF positiva 

%Aquí la frecuencia espacial delta_f_MTF se conserva, lo que cambia es el
%soporte compacto de la matriz MTFp

PixelNumberMTFp = size(MTF1p,2);

Delta_F_MTF1_p = 2*PixelNumberMTFp*delta_f_MTF1; 

% Mostrar el valor
%disp(['Frecuencia espacial para MTF1p y MTF2p: ', num2str(delta_f_MTF1_p),' 1/mm' ]);
disp(['Soporte compacto para MTF1p y MTF2p: ', num2str(Delta_F_MTF1_p),' mm' ]);

% Frecuencia espacial MTF2D

PixelNumerMTF2D1 = size(MTF2D1,1); 
PixelNumerMTF2D2 = size(MTF2D2,1); 

% Delta_F_MTF1_2D = PixelNumerMTF2D1*PixelSize;
% Delta_F_MTF2_2D = PixelNumerMTF2D2*PixelSize;

% delta_f_MTF1_2D = Delta_F_MTF1_positiva/PixelNumerMTF2D1;
% delta_f_MTF2_2D = Delta_F_MTF1_positiva/PixelNumerMTF2D2;


% Mostrar el valor
% disp(['Frecuencia espacial para MTF1_2D: ', num2str(delta_f_MTF1_2D),' 1/mm' ]);
% disp(['Soporte compacto para MTF1_2D: ', num2str(Delta_F_MTF1_2D),' mm' ]);
% disp(['Frecuencia espacial para MTF2_2D: ', num2str(delta_f_MTF2_2D),' 1/mm' ]);
% disp(['Soporte compacto para MTF2_2D: ', num2str(Delta_F_MTF2_2D),' mm' ]);


%% --------------------CÓDIGO CÁLCULO DE LA NPS-------------------------------------------
%
%
%% Función Selección de cortes NPS 

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
end

[K1_noise] = SeleccionarCortes(K1,range); 
[K2_noise] = SeleccionarCortes(K2,range); 

%% Función Covertir a Array, corrección de escala y promedio NPS

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


%% Figuritas 

figure(10) 
subplot(2,2,1), imshow(K1_array(:,:,3),[]);
axis image
title('corte para el Kernel 1 '); 
colorbar
subplot(2,2,2), imshow(K2_array(:,:,3),[]);
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

%% Función Ensamble de ruido (Corte - avg) 

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

            %bigROI(:,:,i)=resta(148:348,148:348,i);
        end
    %Concatenamos todas las ROIs obtenidas 
    ROItotales = cat(3,ROI1,ROI2,ROI3,ROI4,ROI5);
    Padtotal = padarray(ROItotales,[500 500],0,"both");
end

  Padtotal1 = ROIsAndPadding(resta1);
  Padtotal2 = ROIsAndPadding(resta2); 

  % function bigPad = bigROIandPadding(resta)
  %   bigROI = zeros(200,200,size(resta,3));
  %   for i = 1:size(resta,3)
  %       bigROI(:,:,i)=resta(158:357,158:357,i);
  %   end 
  %   bigPad = padarray(bigROI,[500 500],'both');
  % end
     
    % bigPad1 = bigROIandPadding(resta1);
    % bigPad2 = bigROIandPadding(resta2);

%% Pruebas de ROIs
% 
% ROI1 = zeros(50,50,size(resta1,3));
% for i = 1:size(resta1,3)
%             ROI1(:,:,i) = resta1(231:280,231:280,i);
% end
% PadTROI1 = padarray(ROI1,[500 500],'both');
% 
% ROI2 = zeros(50,50,size(resta1,3));
% for i = 1:size(resta1,3)
%             ROI2(:,:,i) = resta1(231:280,331:380,i);
% end
% PadTROI2 = padarray(ROI2,[500 500],'both');
% 
% ROI3 = zeros(50,50,size(resta1,3));
% for i = 1:size(resta1,3)
%             ROI3(:,:,i) = resta1(331:380,231:280,i);
% end
% PadTROI3 = padarray(ROI3,[500 500],'both');
% 
% ROI4 = zeros(50,50,size(resta1,3));
% for i = 1:size(resta1,3)
%             ROI4(:,:,i) = resta1(231:280,95:144,i);
% end
% PadTROI4 = padarray(ROI4,[500 500],'both');
% 
% ROI5 = zeros(50,50,size(resta1,3));
% for i = 1:size(resta1,3)
%             ROI5(:,:,i) = resta1(95:144,231:280,i);
% end
% PadTROI5 = padarray(ROI5,[500 500],'both');
% 
% 
% %% Prueba NPS para cada ROI 
% 
% Dft_prueba_1 = zeros(size(PadTROI1,1),size(PadTROI1,2));
% Dft_prueba_2 = zeros(size(PadTROI2,1),size(PadTROI2,2));
% Dft_prueba_3 = zeros(size(PadTROI3,1),size(PadTROI3,2));
% Dft_prueba_4 = zeros(size(PadTROI4,1),size(PadTROI4,2));
% Dft_prueba_5 = zeros(size(PadTROI5,1),size(PadTROI5,2));
% 
% for i = 1:size(PadTROI1,3)
% Dft_prueba_1(:,:,i) = (abs(fftshift(fft2(PadTROI1(:,:,i))))).^2;
% Dft_prueba_2(:,:,i) = (abs(fftshift(fft2(PadTROI2(:,:,i))))).^2;
% Dft_prueba_3(:,:,i) = (abs(fftshift(fft2(PadTROI3(:,:,i))))).^2;
% Dft_prueba_4(:,:,i) = (abs(fftshift(fft2(PadTROI4(:,:,i))))).^2;
% Dft_prueba_5(:,:,i) = (abs(fftshift(fft2(PadTROI5(:,:,i))))).^2;
% 
% end 
% 
% PixelSize = 0.4454;
% 
% NPStotalPrueba1 = ((PixelSize/50)^2).*mean(Dft_prueba_1,3);
% NPStotalPrueba2 = ((PixelSize/50)^2).*mean(Dft_prueba_2,3);
% NPStotalPrueba3 = ((PixelSize/50)^2).*mean(Dft_prueba_3,3);
% NPStotalPrueba4 = ((PixelSize/50)^2).*mean(Dft_prueba_4,3);
% NPStotalPrueba5 = ((PixelSize/50)^2).*mean(Dft_prueba_5,3);
% 
% figure
% subplot(2,3,1)
% imagesc(NPStotalPrueba1)
% axis image 
% xlabel('Frecuencia en X ')
% ylabel('Frecuencia en Y')
% title('NPS de la ROI1 roja '); 
% subplot(2,3,2)
% imagesc(NPStotalPrueba2)
% axis image 
% xlabel('Frecuencia en X ')
% ylabel('Frecuencia en Y')
% title('NPS de la ROI2 verde');
% subplot(2,3,3)
% imagesc(NPStotalPrueba3)
% axis image 
% xlabel('Frecuencia en X ')
% ylabel('Frecuencia en Y')
% title('NPS de la ROI3 azul ');
% subplot(2,3,4)
% imagesc(NPStotalPrueba4)
% axis image 
% xlabel('Frecuencia en X ')
% ylabel('Frecuencia en Y')
% title('NPS de la ROI4 magenta '); 
% subplot(2,3,5)
% imagesc(NPStotalPrueba5)
% axis image 
% xlabel('Frecuencia en X ')
% ylabel('Frecuencia en Y')
% title('NPS de la ROI5 amarilla '); 
% 

%% Cálculo de NPS

% NPS = (delta_x*delta_y/NxNyM) * SUM (abs(DFT(ROI - ROI_avg)) )
% en otras palabras:
% NPS = (delta_x*delta_y / NxNy) * AVG ((abs(DFT(ROI- ROI_avg))*0.4454^2)^2)

% Resolution:  2.2452 pixels per mm
%Voxel size: 0.4454x0.4454x3 mm^3

Dft_total_1 = zeros(size(Padtotal1,1),size(Padtotal1,2));
Dft_total_2 = zeros(size(Padtotal1,1),size(Padtotal1,2));
% bigDFT_1 = zeros(size(bigPad1,1),size(bigPad1,2));
% bigDFT_2 = zeros(size(bigPad1,1),size(bigPad1,2));


for i = 1:size(Padtotal1,3)
Dft_total_1(:,:,i) = (abs(fftshift(fft2(Padtotal1(:,:,i))))).^2;
Dft_total_2(:,:,i) = (abs(fftshift(fft2(Padtotal2(:,:,i))))).^2;
end

% for i = 1:size(bigPad1,3)
% bigDFT_1(:,:,i) = (abs(fftshift(fft2(bigPad1(:,:,i))))).^2;
% bigDFT_2(:,:,i) = (abs(fftshift(fft2(bigPad2(:,:,i))))).^2;
% end 
% bignps_K1 = ((0.4454/50)^2).*mean(bigDFT_1,3);
% bignps_K2 = ((0.4454/50)^2).*mean(bigDFT_2,3);
 
npsTotal_K1 = ((PixelSize/50)^2).*mean(Dft_total_1,3);
npsTotal_K2 = ((PixelSize/50)^2).*mean(Dft_total_2,3);

%%
figure(11)
subplot(1,2,1), imagesc(npsTotal_K1)
axis image
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')
title('NPS del Kernel 1 '); 
subplot(1,2,2), imagesc(npsTotal_K2)
axis image
title('NPS total del kernel 2');
xlabel('Frecuencia en X ')
ylabel('Frecuencia en Y')

% subplot(2,2,3), imagesc(bignps_K1)
% axis image
% title('NPS total del kernel 1 con ROI grande');
% xlabel('Frecuencia en X ')
% ylabel('Frecuencia en Y')
% 
% subplot(2,2,4), imagesc(bignps_K2)
% axis image
% title('NPS total del kernel 2 con ROI grande');
% xlabel('Frecuencia en X ')
% ylabel('Frecuencia en Y')

%% Conversión de ejes 

%Necesitamos el valor del límite de Nyquist
% w < = 1 / 2u_b

% delta_f = 1 / Delta_F 
% delta_f = 1 / 2* Delta_x = 1 / N*delta_x 
ImageSize = size(Padtotal1);

a = 1050 * PixelSize;
% Límite de Nyquist en ciclos por mm
delta_f = 1 / a;

% Mostrar el valor
disp(['El límite de Nyquist es: ', num2str(delta_f), ' 1/mm']);

[dimX, dimY] = size(npsTotal_K1);  % Dimensiones de la NPS 

%Genera vectores de frecuencias espaciales que van desde −w hasta w 
%en ambos ejes, con una longitud igual al número de píxeles en esos ejes.

freqX = linspace(-delta_f*(dimX/2), delta_f*(dimX/2), dimX);  % Frecuencias en el eje X
freqY = linspace(-delta_f*(dimY/2), delta_f*(dimY/2), dimY);  % Frecuencias en el eje Y

%Creamos el meshgrid
[FreqX, FreqY] = meshgrid(freqX, freqY); 

figure(12)
subplot(1,2,1)
imagesc(freqX, freqY, npsTotal_K1);
axis image;
xlabel('Frecuencia en X (1/mm)');
ylabel('Frecuencia en Y (1/mm)');
title('NPS del Kernel 1 en el Dominio de la Frecuencia');
colorbar;
subplot(1,2,2)
imagesc(freqX, freqY, npsTotal_K2); 
axis image;
xlabel('Frecuencia en X (1/mm)');
ylabel('Frecuencia en Y (1/mm)');
title('NPS del Kernel 2 en el Dominio de la Frecuencia');
colorbar;


%% Figuritas de las ROIS

% imshow(avg1,[]) % Muestra la imagen
% axis image % Ajusta las proporciones del eje
% hold on % Mantén la imagen para graficar sobre ella
% 
% % Coordenadas de las esquinas del cuadrado
% x1 = [231, 280, 280, 231, 231]; % Cierra el cuadrado al volver al primer punto
% y1 = [231, 231, 280, 280, 231]; % Lo mismo aquí para la coordenada 'y'
% 
% x2 = [231, 280, 280, 231, 231];
% y2 = [321, 321, 370, 370, 321]; 
% 
% x3 = [321, 370, 370, 321, 321];
% y3 = [231, 231, 280, 280, 231]; 
% 
% x4 = [231, 280, 280, 231, 231];
% y4 = [141, 141, 190, 190, 141]; 
% 
% x5 = [141, 190, 190, 141, 141];
% y5 = [231, 231, 280, 280, 231];
% 
% x6 = [158, 357, 357, 158, 158];
% y6 = [158, 158,357,357,158];
% 
% plot(x1, y1, 'r', 'LineWidth', 2)
% plot(x2,y2, 'g','LineWidth',2);
% plot(x3,y3,'b','LineWidth',2);
% plot(x4, y4, 'm','LineWidth',2);
% plot(x5,y5,'y','LineWidth',2);
% %plot(x6,y6,'k','LineWidth',3);
% hold off 

%% Cálculo de factores de normalización 

% Sabemos que N* nps = NPS 
% N \Sigma \Sigma nps * (delta_f )^2 = varianza^2

% N = varianza^2 N_x^2 delta_f^2 / \Sigma\Sigma nps
%lo cual es aproximadamente igual a 
% N = var(ROI)^2 N^2 delta_f^2 / \Sigma\Sigma nps

sumatoria1 = sum(npsTotal_K1(:));

sd1 = std2(resta1(231:280,231:280,1));


Nx1 = size(npsTotal_K1,1);


Norma1 = (sd1*(Nx1*PixelSize)^2)/sumatoria1;

sumatoria2 = sum(npsTotal_K2(:));

sd2 = std2(resta1(231:280,231:280,1));


Nx2 = size(npsTotal_K2,1);


Norma2 = (sd2*(Nx2*PixelSize)^2)/sumatoria2;

% Mostrar el valor
disp(['El valor de la norma es: ', num2str(Norma1), ]);

% También N = (delta_x*delta_y / NxNy)
Norma_prima = (PixelSize*PixelSize) / dimX*dimY; 

% Mostrar el valor
disp(['El valor de la norma prima es: ', num2str(Norma_prima), ]);


%% DQE 

width_nps1 = 285;
width_nps2 = 380;

function NPSforDQE = SizeReduction(NPStotal_K,Norma,width_nps)

NPSforDQE = NPStotal_K;

center = round(size(NPSforDQE,1)/2)+1; 
NPSforDQE = Norma.*NPSforDQE((center-width_nps):(center+width_nps),(center-width_nps):(center+width_nps));

end

NPSforDQE1 = SizeReduction(npsTotal_K1,Norma1,width_nps1);
NPSforDQE2 = SizeReduction(npsTotal_K2,Norma2,width_nps2);

function DQE = CalculoDQE(NPSforDQE,MTF2D)

% Calcula el valor máximo de la matriz NPS
NPS_max = max(NPSforDQE(:));

% Define el umbral como el 1% del valor máximo de NPS
threshold = 0.01 * NPS_max;

% Inicializa la matriz DQE
DQE = zeros(size(MTF2D));

% Realiza el cálculo de la DQE solo donde NPS sea mayor o igual al umbral
DQE(NPSforDQE >= threshold) = ((abs(MTF2D(NPSforDQE >= threshold)))).^2 ./ NPSforDQE(NPSforDQE >= threshold);

% Donde el NPS es menor que el umbral, DQE ya será cero por la inicialización

end

DQE1 = CalculoDQE(NPSforDQE1,MTF2D1);
DQE2 = CalculoDQE(NPSforDQE2,MTF2D2);

%% Verificamos la frecuencia espacial de MTF y NPS
% 
% Hacemos el cálculo de deltaf para la nueva matriz de NPS ya que fue
% recortada 
% 
% delta_f = 1 / DeltaF = 1 / N'*delta_x
% 
% NewN = size(NPSforDQE1,1); 
% 
% delta_x = PixelSize; 
% 
% delta_f_NPSforDQE = 1/(NewN*delta_x); 
% 
% Mostrar el valor
% disp(['El valor de la frecuencia espacial para la NPSforDQE es: ', num2str(delta_f_NPSforDQE), ' 1/mm']);
% 
% 
% N_MTF_2D = size(MTF2D1,1); 
% 
% delta_f_MTF_2D = 1/(N_MTF_2D*delta_x); 
% 
% Mostrar el valor
% disp(['El valor de la frecuencia espacial para MTF2D es: ', num2str(delta_f_MTF_2D),' 1/mm' ]);


%% Gráficas DQE en frecuencia espacial 
function frec_esp_DQE = FrecEspDQE(DQE)

PixelSize = 0.4454; %mm
DeltaF_DQE = size(DQE,1)*PixelSize;
delta_f_DQE = 1 / DeltaF_DQE;
W = delta_f_DQE*(size(DQE,1)/2);

% El soporte compacto de la frecuencia espacial es DeltaF
% El soporte compacto es el valor de delta_f por el número de pixeles en la imagen
%Para definir el eje de frecuencias definimos 
frec_esp_DQE = linspace(-W,W,size(DQE,1)); 

end 
frec_esp_DQE1 = FrecEspDQE(DQE1);
frec_esp_DQE2 = FrecEspDQE(DQE2); 

%% Visualizar la DQE 
figure(13)
subplot(1,2,1)
imagesc(frec_esp_DQE1,frec_esp_DQE1,DQE1)
axis image
xlabel('Frecuencia espacial X [1/mm]')
ylabel('Frecuencia espacial Y [1/mm]')
title('DQE1')
subplot(1,2,2)
imagesc(frec_esp_DQE2,frec_esp_DQE2,DQE2)
axis image 
xlabel('Frecuencia espacial X [1/mm]')
ylabel('Frecuencia espacial Y [1/mm]')
title('DQE 2')


%%
figure(14)
subplot(2,2,1)
imagesc(MTF2D1)
axis image 
title('MTF 2D 1')
subplot(2,2,2)
imagesc(NPSforDQE1)
axis image 
title('NPS1')
subplot(2,2,3)
imagesc(MTF2D2)
axis image 
title('MTF 2D 2')
subplot(2,2,4)
imagesc(NPSforDQE2)
axis image 
title('NPS2')


%% Definimos bolita

% ELIP= [
%     0,   5,  5,   0,   0,  0; %fondo
%    1,  0.2,  0.2, 0, 0, 0;];  %circulo central
% 
% phantomtest = phantom(ELIP,50); 
% imshow(phantomtest,[]);
% 
% %% convolución psf y bolita
% 
% % Extraemos la LSF de la primera fila
% lsfprueba = lsf1(1, :);
% 
% % Realizamos la convolución 2D entre la LSF y el círculo
% convResult = conv2(phantomtest, lsfprueba, 'same');
% 
% % Mostrar el resultado de la convolución
% imshow(convResult,[])
% axis equal
% title('Convolución entre PSF y bolita')
% colorbar

%%

% Coordenadas del centro de la matriz
centroDQEx = round(size(DQE1,1)/2)-1;
centroDQEy = round(size(DQE1,1)/2)-1;

dqe_0grados = DQE1(centroDQEx:end,centroDQEy);
dqe_90grados = DQE1(centroDQEx,centroDQEy:end); 

figure(15)
plot(dqe_0grados,'r','DisplayName','0 grados')
hold on 
plot(dqe_90grados,'b','DisplayName','90 grados')
legend
xlabel('Pixeles')
ylabel('Intensidad')
title('DQE a distintos ángulos')

%%

% Coordenadas del centro de la matriz
center_x = round(size(MTF2D1,1)/2)-1;
center_y = round(size(MTF2D1,1)/2)-1;

% Extraer el perfil a lo largo de 0 grados (eje x positivo)
perfil_0grados = MTF2D1(center_x:end, center_y);
perfil_90grados = MTF2D1(center_x,center_y:end);

% Puedes ahora graficar para ver si ambos perfiles coinciden
figure(16)
subplot(1,2,1)
plot(perfil_0grados, 'r', 'DisplayName', 'MTF2D1 0 grados');
hold on;
plot(MTF_positiva1(1,:), 'b--', 'DisplayName', 'MTFpositiva 1');
legend;
xlabel('Pixeles');
ylabel('Valor MTF');
title('Comparación del perfil MTF en 0 grados');

subplot(1,2,2)
plot(perfil_90grados, 'r', 'DisplayName', 'MTF2D1 90 grados');
hold on;
plot(MTF_positiva1(90,:), 'b--', 'DisplayName', 'MTFpositiva 1');
legend;
xlabel('Pixeles');
ylabel('Valor MTF');
title('Comparación del perfil MTF en 90 grados');
