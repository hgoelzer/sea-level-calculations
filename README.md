# sea-level-calculations
## Compare sea-level calculation methods

Heiko Goelzer (heig@norceresearch.no), Feb 2023

## Sea-level calculation based on 
Vaf: Volume above floatation <br> 
H1991: https://doi.org/10.3189/S0260305500008387 <br>
A2020: https://tc.copernicus.org/articles/14/2819/2020/; Based on provided example https://doi.org/10.7910/DVN/9LUJTD by Surendra Adhikari, Jet Propulsion Lab, Caltech <br>
G2020: https://tc.copernicus.org/articles/14/833/2020/ <br>

## Test Vaf method for a number of simplified configurations
slc_driver_vaf.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  vaf_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  plot_columns.m <br>

## Test H1990 for a number of simplified configurations
slc_driver_h1990.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  h1990_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  plot_columns.m <br>

## Test A2020 for a number of simplified configurations
slc_driver_a2020.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  a2020_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  plot_columns.m <br>

## Test G2020 for a number of simplified configurations
slc_driver_g2020.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  g2020_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  plot_columns.m <br>


## Compare Vaf, H1990, A2020 and G2020 for a number of simplified configurations
slc_driver_comp.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  vaf_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  h1990_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  a2020_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  g2020_func.m <br>
&nbsp;&nbsp;&nbsp;&nbsp;  plot_columns.m <br>

