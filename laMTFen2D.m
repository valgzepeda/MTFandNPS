
prepMTF; % Ejecutamos el archivo que prepara todos los datos bonitos

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

