prepNPS; 

%% Función Ensamble de ruido (Corte - avg) 

function resta = Ensamble(K_array_correc,avg)
    resta = zeros(size(K_array_correc,1),size(K_array_correc,2));

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
        end
    %Concatenamos todas las ROIs obtenidas 
    ROItotales = cat(3,ROI1,ROI2,ROI3,ROI4,ROI5);
    Padtotal = padarray(ROItotales,[500 500],0,"both");
end

  Padtotal1 = ROIsAndPadding(resta1);
  Padtotal2 = ROIsAndPadding(resta2); 

%% Cálculo nps (sin normalizar)

% NPS = (delta_x*delta_y/NxNyM) * SUM (abs(DFT(ROI - ROI_avg)) )
% en otras palabras:
% NPS = (delta_x*delta_y / NxNy) * AVG ((abs(DFT(ROI- ROI_avg))*0.4454^2)^2)

function npsTotal = Calcularnps(Padtotal)

Dft_total = zeros(size(Padtotal,1),size(Padtotal,2));

for i = 1:size(Padtotal,3)
Dft_total(:,:,i) = (abs(fftshift(fft2(Padtotal(:,:,i))))).^2;
end

npsTotal = mean(Dft_total,3); %Así viene en el paper  

end

npsTotal_K1 = Calcularnps(Padtotal1); 
npsTotal_K2 = Calcularnps(Padtotal2); 

%% nps de ambos kernels 
figure(11)
subplot(1,2,1), imagesc(npsTotal_K1)
axis image
xlabel('Pixel')
ylabel('Pixel')
title('NPS del Kernel 1 '); 
subplot(1,2,2), imagesc(npsTotal_K2)
axis image
title('NPS total del kernel 2');
xlabel('Pixel')
ylabel('Pixel')

%% Cálculo de factores de normalización 

% Sabemos que N* nps = NPS 
% N \Sigma \Sigma nps * (delta_f )^2 = varianza^2

% N = varianza^2 N_x^2 delta_f^2 / \Sigma\Sigma nps
%lo cual es aproximadamente igual a 
% N = var(ROI) N^2 delta_f^2 / \Sigma\Sigma nps

PixelSize = 0.4454;

function Norma = Normalizarnps(K_array_correc,npsTotal)

PixelSize = 0.4454;
[dimX,~] = size(npsTotal); 

varianza = var(K_array_correc(231:280,231:280,10), 0, 'all'); 

sumatoria = sum(npsTotal(:));

Norma = (varianza*(dimX*PixelSize)^2)/sumatoria;
end

Norma1 = Normalizarnps(K1_array_correc,npsTotal_K1);
Norma2 = Normalizarnps(K2_array_correc,npsTotal_K2);

% Mostrar el valor
disp(['El valor de la norma correspondiente al kernel 1 es: ', num2str(Norma1), ]);

[dimX,dimY] = size(npsTotal_K1);
% También N = delta_x*delta_y / NxNy de la roi 
Norma_prima = (PixelSize*PixelSize) / (50*50); 

% Mostrar el valor
disp(['El valor de la norma prima es: ', num2str(Norma_prima), ]);

% Mostrar el valor
disp(['El valor de la norma correspondiente al kernel 2 es: ', num2str(Norma2), ]);



%% Cálculo de frecuencia espacial 

%Necesitamos el valor del límite de Nyquist
% w < = 1 / 2u_b

% delta_f = 1 / Delta_F = 1 / 2* Delta_x = 1 / N*delta_x 

PixelSize = 0.4454;
ImageSize = size(npsTotal_K1,1);

% delta_f ciclos por mm
delta_f = 1 /( ImageSize*PixelSize);

Delta_F = delta_f*ImageSize;


% Mostrar el valor
disp(['La frecuencia espacial de la nps es: ', num2str(delta_f), ' 1/mm']);

% Mostrar el valor
disp(['El soporte compacto de la nps es: ', num2str(Delta_F), ' mm']);