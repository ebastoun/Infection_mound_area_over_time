# Infection_mound_area_over_time

ReadMe

This script is meant to analyse the area of infected cells in cell monolayers over time from live cell imaging data of the bacterial fluorescence. The output is a binary mask of the infected cell area over time and the area over time plotted and as .xlsx file. 

Involved steps are: 
1. Binarize bacterial fluorescence images
2. Get border of infected cell area with alpha shapes

Created: 15/11/2023 by Lara Hundsdorfer, with help from Julio César Sánchez Rendón and Marie Münkel.

INPUT DATA: Imaging data of the bacterial fluorescence as multi-tif file. 
Example: we provide an example, see "my_sample.tif" file
