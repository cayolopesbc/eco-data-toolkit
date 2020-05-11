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
        self.data = ds.to_dataframe().reset_index()
        self.variables = self.data.columns
        
#        try:
#            self.geometry  = [shapely.geometry.Point(x,y) for x, y in zip(self.data['longitude'], self.data['latitude'])]
#        except KeyError:
#            self.geometry  = [shapely.geometry.Point(x,y) for x, y in zip(self.data['lon'], self.data['lat'])]
            

    def var_choice(self, variables = []):
        
        # Checks whether a datetime column exists. If it exists, 
        # it will be considered as a timeIndex.
        
        timestamp = [self.data.select_dtypes(include=[np.datetime64]).keys()[0]]
        if len(timestamp)>0:
            self.data = self.data.set_index(self.data.select_dtypes(include=[np.datetime64]).keys()[0])
        
        #Variables list and choice:
        
        if not variables:            
            selected = []
            variables = np.array(self.data.columns)
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
        
        self.data = self.data[selected]
        self.data = self.data.reset_index()
        try:
            self.data = self.data.drop('index')
        except:
            pass
        
        
        return
    
    def centroidFilter(self, mask_path, radius = 1):
        
        #Extraindo os dados  mais proximos da regiao de estudo para diminuir a quantidade
        # de dados e agilizar o processamento:
        
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
            self.data = self.data[(lon_filter_down[0] >= self.data['longitude']) & (self.data['longitude'] >= lon_filter_up[0]) & (lat_filter_down[0] >= self.data['latitude']) & (self.data['latitude'] >= lat_filter_up[0])]
        except KeyError:
            self.data = self.data[(lon_filter_down[0] >= self.data['lon']) & (self.data['lon'] >= lon_filter_up[0]) & (lat_filter_down[0] >= self.data['lat']) & (self.data['lat'] >= lat_filter_up[0])]
            self.data['longitude'] = self.data['lon']
            self.data['latitude'] = self.data['lat']
            
        geom = [shapely.geometry.Point(x,y) for x, y in zip(self.data['longitude'], self.data['latitude'])]
        self.data = gpd.GeoDataFrame(self.data, geometry = geom)
        self.data.crs = mask_shape.crs
        return
        
    def pointFilter(self, Coord):
        
        #Find the nearest grid point to specified coordinates 
        # and extract data.
        
        if isinstance(Coord,list):
            Point = shapely.geometry.Point(Coord[0], Coord[1])
        elif isinstance(Coord, shapely.geometry.Point):
            Point = Point
        else:
            print("Invalid coordinates format.")
            return
            
        try:
            geom = [shapely.geometry.Point(x,y) for x, y in zip(self.data['longitude'], self.data['latitude'])]
        except KeyError:
            geom = [shapely.geometry.Point(x,y) for x, y in zip(self.data['lon'], self.data['lat'])]   
            
        geom = list(set(geom))

        geomei = shapely.ops.nearest_points(Point.centroid, shapely.geometry.MultiPoint(geom))
        self.data = self.data.geometry.intersects(geomei[1])
        return self.data
     
    def temporalFilter(self, StartDate, EndDate):
        
        # Filter time series between specific dates.
        
        self.data = self.data[StartDate:EndDate]
        return 




  
def resampleMI(data, freq = '0.5H', method = 'ffill'):
# .mean()
# .ffill()
# .bbfill()
# .sum()    
    try:
        data_crs = data.crs
    except AttributeError:
        data_crs = {'init': 'epsg:4326'}
    
    #Converte GeoDataFrame em DataFrame para utilizar o metodo .resample():
    if isinstance(data, pd.DataFrame):    
        data = pd.DataFrame(data)
    
    #Filtra as geometrias para iniciar o resample:
    geometry = data[['geometry']].drop_duplicates(subset = ['geometry'])
    
    #De acordo com a geometria (coordenadas) e feito o resample de sua respectiva
    #serie temporal e posteriormente todas as series novas sao unidas em um novo
    #dataframe
    
    resamples = []
    for geom in geometry['geometry']:
        filt = data[(data['geometry'] == geom)]
        resample_data = filt.resample(freq).ffill()
        resamples.append(resample_data)
    df = pd.concat(resamples)

    #Convertendo o novo DF em GeoDataFrame:    
    df = df.drop(['geometry'], axis = 1)

    try:
        geom = [Point(x,y) for x, y in zip(df['longitude'], df['latitude'])]
    except KeyError:
        geom = [Point(x,y) for x, y in zip(df['lon'], df['lat'])]
        
    df = gpd.GeoDataFrame(df, crs = data_crs, geometry = geom)
        
    return df     

def pointFilter(df, Coord):
    geom = df.geometry.drop_duplicates()
    Point = Point(Coord[0], Coord[1])
    geomei = shapely.ops.nearest_points(Point.centroid, MultiPoint(geom.values))
    return df[df.geometry.intersects(geomei[1])]

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
EndDate   =  "02/21/2000"

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

