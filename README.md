## Eco-Data Manage Toolkit  [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18914.svg)](https://zenodo.org/record/3941213#.XwtdeyhKhPY)
Eco-Data Manage Toolkit is a Python toolkit is a set of tools to facilitate data management for hydrology/limnology applications.

 In this first version, it is only capable of reading data in NetCDF format. Its features include:
- Reading data
- Spatial and temporal filters
- Spatial and temporal interpolation
- Export to .tif and .csv formats

## Installation
To install lastest release version, use `pip install ecodatatk`.

### Microsoft Windows installation
If the GeoPandas or/and Rasterio package is not installed, it is suggested that it be installed manually using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html).<br><br>
Geopandas package install:<br>
`conda install -c conda-forge geopandas` <br><br>
Rasterio package install:<br>
`conda install rasterio`<br><br>

Then, the package can be installed without any error:<br>
`pip install ecodatatk`

## Application Examples
Applications can be found in Jupyter Notebooks, in the `./Examples` directory in this repository.

## Help
You can use the `help()` to get more information about the methods uses and parameters.

## License
[BSD 3-Clause License](LICENSE)

## How to cite
Cayo Lopes. (2020, July 12). Eco-Data Manage Toolkit (v. 0.0.1) (Version 0.0.1). Zenodo. http://doi.org/10.5281/zenodo.3941213
