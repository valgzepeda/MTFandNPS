%% Lectura de archivos DICOM 
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

% Leer volúmenes DICOM Kernel 1 y Kernel 2 y aplicar la corrección
for i = 1:num_files
    % Definir la ruta de la carpeta
    folder1 = [base_path 'Abdomen3_' num2str(file_numbers1(i)) '\'];
    folder2 = [base_path 'Abdomen3_' num2str(file_numbers2(i)) '\'];
    
    % Obtener el archivo DICOM dentro de la carpeta
    dicom_files1 = dir([folder1 '*.dcm']);
    dicom_files2 = dir([folder2 '*.dcm']);
    
    % Verificar que se haya encontrado al menos un archivo DICOM en cada carpeta
    if isempty(dicom_files1)
        error(['No se encontraron archivos DICOM en la carpeta: ' folder1]);
    end
    if isempty(dicom_files2)
        error(['No se encontraron archivos DICOM en la carpeta: ' folder2]);
    end
    
    % Usar el primer archivo encontrado para obtener los metadatos
    dicom_file1 = [folder1 dicom_files1(1).name];
    dicom_file2 = [folder2 dicom_files2(1).name];
    
    % Leer el volumen DICOM
    K1{1,i} = squeeze(dicomreadVolume(folder1));
    K2{1,i} = squeeze(dicomreadVolume(folder2));
    
    % Leer los metadatos para obtener RescaleSlope y RescaleIntercept
    info1 = dicominfo(dicom_file1);
    info2 = dicominfo(dicom_file2);
    
    % Extraer los valores de RescaleSlope y RescaleIntercept
    rescaleSlope1 = info1.RescaleSlope;
    rescaleIntercept1 = info1.RescaleIntercept;
    rescaleSlope2 = info2.RescaleSlope;
    rescaleIntercept2 = info2.RescaleIntercept;

end
