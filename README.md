# sea-level-calculations
Compare sea-level calculation methods

Heiko Goelzer (heig@norceresearch.no), Feb 2023

% Sea-level calculation based on 
% G2020: https://tc.copernicus.org/articles/14/833/2020/
% A2020: https://tc.copernicus.org/articles/14/2819/2020/
% Based on provided example https://doi.org/10.7910/DVN/9LUJTD
% Surendra Adhikari, Jet Propulsion Lab, Caltech, adhikari@jpl.nasa.gov

% Test A2020 for a number of simplified configurations
sle_driver_a2020.m
  a2020_func.m
  plot_columns.m


% Test G2020 for a number of simplified configurations
sle_driver_g2020.m
  g2020_func.m
  plot_columns.m


% Compare A2020 and G2020 for a number of simplified configurations
sle_driver_comp.m
  a2020_func.m
  g2020_func.m
  plot_columns.m

