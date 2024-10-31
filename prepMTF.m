dataprep; %Ejecuta el archivo dataprep que extrae los dicom a una celda por kernel

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

%% Función Concatenar ROIs y ajustar escala de HU 

function avg = ConcatEscala(K_rois, num_files, range_mtf,rescaleSlope,rescaleIntercept)
    
    % Inicializar una matriz 3D para almacenar las ROIs concatenadas
    K_array = zeros(size(K_rois{1,1},1),size(K_rois{1,1},2), length(range_mtf)*num_files);

        % Rellenar la matriz K_array con las ROIs concatenadas
        for i = 1:num_files
            start_slice = (i-1)*length(range_mtf) + 1;
            end_slice = i*length(range_mtf);
            K_array(:,:,start_slice:end_slice) = K_rois{1,i};
        end

    %Corrección UH = RescaleSlope*(pixel value) + RescaleIntercept
    K_array_correc = rescaleSlope.*K_array+ rescaleIntercept;
    avg = mean(K_array_correc,3);

end

avg1 = ConcatEscala(K1_rois,num_files,range_mtf,rescaleSlope1,rescaleIntercept1);
avg2 = ConcatEscala(K2_rois,num_files,range_mtf,rescaleSlope2,rescaleIntercept2);

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
%%
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
