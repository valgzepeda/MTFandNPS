1. dataprep manipula los archivos dicom colocándolos en celdas en función del número de adquisiciones para cada kernel.
2. prepNPS y prepMTF requieren de dataprep para funcionar.
3. prep NPS manipula las celdas provenientes de dataprep para seleccionar los cortes de interés, corregir por UH y obtener el valor promedio en el eje axial de los cortes seleccionados.
4. theNPS utiliza los datos de prep NPS para obtener el ensamble de ruido y calcular la NPS para las ROIS de interés. También calcula las constantes de normalización, la frecuencia espacial y el soporte compacto.
