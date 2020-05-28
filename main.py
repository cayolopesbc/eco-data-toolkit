# -*- coding: utf-8 -*-
"""
Created on Sun May 10 23:22:28 2020

@author: cayoh
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:15:14 2020

@author: cayoh
"""
import xarray as xr
import rioxarray as rxy
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import datetime

from shapely.geometry import Point, MultiPoint

import shapely
################################################################################


class SpatialDataset:
    
    def __init__(self, datapath):    
        #Data Read:
        ds   = xr.open_dataset(datapath)
        self._data = ds.to_dataframe().reset_index() 
        if self._data.select_dtypes(include=[np.datetime64]).keys()[0]:
            self._data = self._data.set_index(self._data.select_dtypes(include=[np.datetime64]).keys()[0])
        self._data['geometry'] =  self._getgeom()
        self._variables = self._data.columns
        self._crs  = {'init': 'epsg:4326'}
        
    @property
    def data(self):
        return self._data
    
    @property
    def variables(self):
        return self._variables
    
    def _df2geodf(self, geom, crs = None):
        if crs != self._crs:
            self._crs = crs
        return gpd.GeoDataFrame(self._data, crs = self._crs, geometry = geom)
    
    def _getgeom(self):        
        try:
            geom = [Point(x,y) for x, y in zip(self._data['longitude'], self._data['latitude'])]
        except KeyError:
            geom = [Point(x,y) for x, y in zip(self._data['lon'], self._data['lat'])]
            self._data['longitude'] = self._data['lon']
            self._data['latitude']  = self._data['lat'] 
        return geom
        
    def var_choice(self, variables = []):

        #Variables list and choice:
        if not variables:            
            selected = []
            variables = np.array(self._data.columns)
            i = 1
            print("The file have following variables:\n")

            for var in variables:
                if var in ['latitude', 'longitude', 'lon', 'lat']:
                    selected.append(var)
                    variables = variables[variables!=var]
                else:    
                    print('{} - {}'.format(i,var))
                    i+=1          
            print('{} - All variables'.format(i))
            
            print("Choose variables by number, when enough type 'q' and press enter:\n")
            q = " "
            while q.upper() != "Q":
                q = input("Variable code: ")       
                if q == str(i):
                    selected  = list(selected) + list(variables)
                    q = 'Q'
                elif q.upper() != "Q":
                    selected.append(variables[int(q)-1])
        else:
            selected = list(variables)
        
        self._data = self._data[selected]

        try:
            self._data = self._data.drop('index')
        except:
            pass
        
        self._variables = self._data.columns
        return
    
    def centroidFilter(self, mask_path, radius = 1):
        '''
        It extracts the data closest to the study area to decrease 
        the amount of data and speed up processing:
        '''
        try:
            mask_shape = gpd.GeoDataFrame.from_file(mask_path)
        except:
            print("Invalid Shape or non-existent path.")
            return
        
        lat_filter_up   = mask_shape.centroid.y - radius
        lat_filter_down = mask_shape.centroid.y + radius
        lon_filter_up   = mask_shape.centroid.x - radius
        lon_filter_down = mask_shape.centroid.x + radius
        
        try:
            self._data = self._data[(lon_filter_down[0] >= self._data['longitude']) & (self._data['longitude'] >= lon_filter_up[0]) & (lat_filter_down[0] >= self._data['latitude']) & (self._data['latitude'] >= lat_filter_up[0])]
        except KeyError:
            self._data = self._data[(lon_filter_down[0] >= self._data['lon']) & (self._data['lon'] >= lon_filter_up[0]) & (lat_filter_down[0] >= self._data['lat']) & (self._data['lat'] >= lat_filter_up[0])]
            self._data['longitude'] = self._data['lon']
            self._data['latitude'] = self._data['lat']
            
        geom =  self._getgeom()
        crs = mask_shape.crs
        self._data = self._df2geodf(geom, crs)
        return
        
    def pointFilter(self, Coord):
        '''
        Find the nearest grid point to specified coordinates 
         and extract data.
        '''
        if isinstance(Coord,list):
            Point = shapely.geometry.Point(Coord[0], Coord[1])
        elif isinstance(Coord, shapely.geometry.Point):
            Point = Point
        else:
            print("Invalid coordinates format.")
            return
            
        if isinstance(self._data,gpd.GeoDataFrame):
            geom = list(self._data.geometry)           
        else:
            geom = self._getgeom()
            self._data = self._df2geodf(geom)
                
        geomeInt = shapely.ops.nearest_points(Point.centroid, shapely.geometry.MultiPoint(geom))
        self._data = self._data[self._data.geometry.intersects(geomeInt[1])]
        return
     
    def temporalFilter(self, StartDate, EndDate):
        '''
        Makes a temporal filter in dataset.
        '''
        if isinstance(self._data.index,pd.DatetimeIndex):
            self._data = self._data[StartDate:EndDate]
        else:
            print('No has DateTime index in the dataset.')
        return 

    def resample(self, freq = '1H', method = 'ffill', order = 1):

        #Converte GeoDataFrame em DataFrame para utilizar o metodo .resample():
        if isinstance(self._data, gpd.GeoDataFrame):    
            self._data = pd.DataFrame(self._data)
          
        #Filtra as geometrias para iniciar o resample:
        geometry = self._data['geometry'].drop_duplicates()
        
        #De acordo com a geometria (coordenadas) e feito o resample de sua respectiva
        #serie temporal e posteriormente todas as series novas sao unidas em um novo
        #dataframe
        resamples = []
        for geoPoint in geometry:
            
            resample_data = self._data[(self._data['geometry'] == geoPoint)].resample(freq).mean()
            
            if method == 'polynomial' or method == 'spline':
                resample_data = resample_data.interpolate(method = method, order = order)   
            elif method=='backfill' or method=='bfill':
                resample_data = resample_data.bbfill()            
            else:
                resample_data = resample_data.ffill()
            resamples.append(resample_data)  

        self._data = pd.concat(resamples)       
        
        geom = self._getgeom()
        self._data = self._df2geodf(geom)
        return



# Preenche as falhas das séries, sejam essas falhas oriundas dos dados de origem ou devido ao aumento da frequência utilizando a função groupData
# Métodos de preenchimento dos dados faltantes:
# Apenas substituicao:
# pad / ffill: propagate last valid observation forward to next valid 
# backfill / bfill: use next valid observation to fill gap.
        
# 'ffill'       - Propaga um valor conhecido "para frente" até um encontrar um novo valor conhecido.
# 'backfill'    - Propaga um valor conhecido "para trás" até um encontrar um novo valor conhecido.
# Interpolacoes:
# 'polynomial'  - Interpolação polinomial, deve se informar uma ordem (order = n): a) ordem máxima para os dados: 4; b) Ordem = 1 retorna uma interpolacao linear.
# 'spline'      - Interpolação usando splines, deve se informar uma ordem (order = n): a) ordem máxima para esses dados: 5; sugiro ordem 1 ou 3, melhor ler mais sobre o método.



def cropRst(raster, mask_shp, out_tif = None, remove = False, buffer_mask = 0):
    if isinstance(mask_shp,gpd.GeoDataFrame):
      pass  
    else:
        mask_shp = gpd.GeoDataFrame.from_file(mask_shp)
        mask_shp.geometry[0] = mask_shp.geometry[0].buffer(buffer_mask)
   
    if out_tif == None:
        out_tif = raster.replace('.tif','_crop.tif')
    
    
    ds = rxy.open_rasterio(raster, masked = True, chunks = True)
    clipped = ds.rio.clip(mask_shp.geometry.apply(mapping), mask_shp.crs, drop=False, invert=False)
    clipped.rio.to_raster(out_tif)
    ds.close() 
    
    if remove:
        os.remove(raster) 
        print("Aqui02")
          
    return out_tif


def pts2Raster(df, var_list, base_name, out_dir, outputBounds = [-52.40,-33.80,-53.80,-32.00], outCRS = 'WGS84', buffer_mask = 0, shp_mask = None):
    outs = []
    for varName in var_list:
        vrt_fn = os.path.join(out_dir,varName+'Vrt.vrt')
        lyr_name = varName
        out_tif = os.path.join(out_dir, base_name + '_'+ varName +'.tif')
        tempPath = os.path.join(out_dir, varName +'.csv')
        df[[varName,'latitude','longitude']].to_csv(tempPath,header = True, index = False)
        
        with open(vrt_fn, 'w') as fn_vrt:
            fn_vrt.write('<OGRVRTDataSource>\n')
            fn_vrt.write('\t<OGRVRTLayer name="%s">\n' % lyr_name)
            fn_vrt.write('\t\t<SrcDataSource>%s</SrcDataSource>\n' % tempPath)
            fn_vrt.write('\t\t<SrcLayer>%s</SrcLayer>\n' % lyr_name)
            fn_vrt.write('\t\t<GeometryType>wkbPoint</GeometryType>\n')
            fn_vrt.write('\t\t<GeometryField encoding="PointFromColumns" x="longitude" y="latitude" z="%s"/>\n'  %varName)
            fn_vrt.write('\t</OGRVRTLayer>\n')
            fn_vrt.write('</OGRVRTDataSource>\n')
        
#        https://gdal.org/programs/gdal_grid.html
        gridOp = gdal.GridOptions(format = 'Gtiff', outputBounds = outputBounds, algorithm = 'linear:radius=0.0:nodata = -9999', outputSRS = outCRS)
        
        if isinstance(shp_mask,gpd.GeoDataFrame):
            temp_tif = os.path.join(out_dir,"Temp_" + base_name + '_'+ varName +'.tif')
            gdal.Grid(temp_tif, vrt_fn, options = gridOp)
            cropRst(temp_tif, shp_mask, out_tif, remove = True, buffer_mask = buffer_mask)
        else:
            gdal.Grid(out_tif, vrt_fn, options = gridOp)
        
        os.remove(tempPath)
        outs.append(out_tif)

    return outs


################################
    
#Data Path/Name:
ei_path = 'D:\\OneDrive\\Arquivos OLD PC\\Cafofo\\Balanco de Calor\\ERA-Interim'
ei_name = 'Mangueira2000_01_01to2018_01_01.nc'
ei_data = os.path.join(ei_path, ei_name) 


#Data Read:
ds = xr.open_dataset(ei_data)
ei = ds.to_dataframe()
ei = ei.reset_index()

#Variables list/choose:
ei = variablesChoose(ei)
ei = ei.set_index('time')

#Convert data to GeoDataframe, set coordinates reference system code (crs) 
#(default is WGS84 - epgs:4326 ):
crs = {'init': 'epsg:4326'}
geomei = [Point(x,y) for x, y in zip(ei['longitude'], ei['latitude'])]
ei = gpd.GeoDataFrame(ei, crs = crs, geometry = geomei)

#Temporal Filter:
#Insert dates on format Month/Day/Year
StartDate = "01/21/2000"
EndDate   =  "01/21/2000"

ei = temporalFilter(df, StartDate, EndDate)

#Spatial Filter:

# 1º Mode - By shapes centroid 
# The spatial resolution in Data Grid may be roughly and an filter inside only
# study shape area can be return no enough points to adequate analysis. 
# For instance, radial filter can provide more points order to make a spatial interpolation more 'robust':

maskShapePath = "D:\\OneDrive\\Arquivos OLD PC\\Cafofo\\Balanco de Calor\\MangueiraShp\\polimirim.shp"
maskShape = gpd.GeoDataFrame.from_file(maskShapePath)

# Set radius value in radialFilter function. The value depend on crs (in this case radius = 1 is equivalent 1 degree)
ei = radialFilter(ei, maskShape, radius = 1)

# 2º Mode -  By Specific Point:
# Search the nearest point insert in Data Grid in relation the coordinates set.
# Set point coordinate in variable 'coord'(Lon/Lat):
coords = [-50.75, -33.55]
ei = pointFilter(ei, coords)

# Resample data:
# Resample temporal step. Set frequency value in "freq" variable.
# Frequencies: H = hourly, D = dayly, M = Monthly
# For example, freq = '0.5H' the data is resample in 30 minutes data (0.5 hour).

# Set method for resample:
# mean  - The data is set as mean of period. 
# ffill - 'FowardFill': The data is set as first data no null between two timestep. 
# bfill - 'BackFill': The data is set as second data no null between two timestep.
# sum   - Make values acumulation between timesteps. 

ei = resampleMI(ei, freq ='0.5H',  method = 'ffill')


# Rasterization:
# The data is interpolate on rectangular raster with bounds definite by user through outputBounds variables.
# outputBounds = [Longitude-Left Upper Point, Latitude-Left Upper Point, Longitude-Right Lower Point, Latitude-Right Lower Point]

