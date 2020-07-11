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
import shapely
from shapely.geometry import Point, MultiPoint, mapping
import gdal
import rasterio
################################################################################


class NetCDFSpatial:
    
    def __init__(self, datapath, geom_keys = None,crs ={'init': 'epsg:4326'} ):    
        '''
        class constructor
        '''
        ds   = xr.open_dataset(datapath)
        self._data = ds.to_dataframe().reset_index() 
        if self._data.select_dtypes(include=[np.datetime64]).keys()[0]:
            self._data = self._data.set_index(self._data.select_dtypes(include=[np.datetime64]).keys()[0])
        
        self._data['geometry'] =  self._getgeom()
        self._crs  = crs
        
        self._data = self._df2geodf(self._data['geometry'], crs = self._crs)
        self._variables = list(self._data.columns.drop(['geometry']))
        
        
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
    
    def _getgeom(self, geom_keys = None):
        
        if geom_keys:
            geom = [Point(x,y) for x, y in zip(geom_keys[0], geom_keys[1])]
            self._data['longitude'] = self._data[geom_keys[0]]
            self._data['latitude']  = self._data[geom_keys[1]]
        else:
            try:
                geom = [Point(x,y) for x, y in zip(self._data['longitude'], self._data['latitude'])]
            except KeyError:
                try:
                    geom = [Point(x,y) for x, y in zip(self._data['lon'], self._data['lat'])]
                    self._data['longitude'] = self._data['lon']
                    self._data['latitude']  = self._data['lat']
                except KeyError:
                    geom = [Point(x,y) for x, y in zip(self._data['x'], self._data['y'])]
                    self._data['longitude'] = self._data['x']
                    self._data['latitude']  = self._data['y']
             
        return geom
        
    def var_choice(self, variables = []):
        
        '''
        This method selects the data to be handled.
         
        input:
            :param variables = shapefile

        output:
            return: variables dataset (pandas.GeoDataFrame instance)
         
        '''

        if not variables:            
            selected = []
            variables = np.array(self._data.columns)
            i = 1
            print("Choose the variables from the numbers.\n When finished the choice type 'q'.\n")
            print("The file has the following variables:\n")

            for var in variables:
                if var in ['latitude', 'longitude', 'lon', 'lat']:
                    selected.append(var)
                    variables = variables[variables!=var]
                else:    
                    print('{} - {}'.format(i,var))
                    i+=1          
            print('{} - All variables'.format(i))           
    
            q = " "
            while q.upper() != "Q":
                q = input("Variable code: ")       
                if q == str(i):
                    selected  = list(selected) + list(variables)
                    q = 'Q'
                elif q.upper() != "Q":
                    selected.append(variables[int(q)-1])
            selected = list(set(selected))
        else:
            selected = list(variables) + ['longitude','latitude','geometry']
        
        self._data = self._data[selected]

        try:
            self._data = self._data.drop('index')
        except:
            pass
        
        self._variables = self._data.columns
        return
    
    def centroidFilter(self, mask_path, radius = 1):
        
        '''
        This method extracts the data in a search radius defined  from the 'mask shapefile' centroid.
         
        input:
            :param mask_path = shapefile
            :param radius = limit radius to buffer in the centroid shapefile mask, 
                            radius default value is equal = 1 degree
        output:
            return: filtered spatial dataset (pandas.GeoDataFrame instance)
            
        info:
             - in some case the study area has a much smaller scale than the data 
             used, so filtering it using only the points entered in a region can 
             return null values.Therefore, this is an alternative method of spatial 
             filtering.
         
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

        self._data = self._data[(lon_filter_down[0] >= self._data['longitude']) & (self._data['longitude'] >= lon_filter_up[0]) & (lat_filter_down[0] >= self._data['latitude']) & (self._data['latitude'] >= lat_filter_up[0])]
        
        geom =  self._getgeom()
        crs = mask_shape.crs
        self._data = self._df2geodf(geom, crs)
        return
        
    def pointFilter(self, Coord):
        
        '''
        This method find the nearest grid point to specified coordinates 
         and extract data.
         
        input:
            :param Coord =  coordinates latitude/Longitude (can be set as a 
                             list instance in [lat,Long] format or a 
                             shapely.geometry.Point instance)
        output:
            return: filtered spatial dataset (pandas.GeoDataFrame instance)
         
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
       This method applies a temporal filter in dataset.
        
        input:
            :param StarDate =  initial timestamp 
            :param EndDate = end timestamp
        output:
            return: filtered temporal dataset (pandas.GeoDataFrame instance)
        
        '''
        if isinstance(self._data.index,pd.DatetimeIndex):
            self._data = self._data[StartDate:EndDate]
        else:
            print('No has DateTime index in the dataset.')
#        self._data = self._df2geodf(self._data['geometry'])
        return 

    def resampleDS(self, freq = '1H', method = 'ffill', order = 1):
        '''
        This method applies a temporal resample.
        
        input:
            :param freq = resample frequencie, default frequencie is set as freq = '1H'
                          frequencies: H = hourly, D = dayly, M = Monthly
            
            :param method = resample method  (forwardfill, backfill, polynomial and spline methods is able),
                            default method is set as 'fowardfill'
                            
            :param order = remsample method order (is used only in 'polynomial' and 'spline' methods)
                            default order is set as order = 1
        output:
            return: filtered temporal dataset (pandas.GeoDataFrame instance)
            
        info:
            resample case:
            - in downscale resample,the outputsare defined by the period average
            - in upscale resample, the outputs are defined using the chosen method
        
            resample methods:
                - fowardfill/ ffill: propagate last valid observation forward to next valid.
                - backfill / bfill: use next valid observation to fill gap.
                - polynomial: fill data through polynomial interpolation
        
        '''

        if isinstance(self._data, gpd.GeoDataFrame):    
            self._data = pd.DataFrame(self._data)

        geometry = self._data['geometry'].drop_duplicates()
 
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
    
    def _maskCheck(self, mask_shape, buffer_mask):
        if not isinstance(mask_shape,gpd.GeoDataFrame):
            mask_shape = gpd.GeoDataFrame.from_file(mask_shape)
            mask_shape.geometry[0] = mask_shape.geometry[0].buffer(buffer_mask)
        return mask_shape
        
    def pts2csv(self, out_name):
        '''
        This method exports the points dataset to comma separete file (*.csv).
        
        input:
            param: out_name = output file name.
        output:
            return: .csv file.
        
        info:
            This method can apply in cases where its wish to work with data at one point.
            
        '''
        
        aux = self._data
#        aux.drop('geometry', axis = 1, inplace = True)
        aux.to_csv(out_name, sep =';')
        return

    def _cropRst(self,raster, mask_shp, out_tif = None, remove = False, buffer_mask = 0):
        
        mask_shp = self._maskCheck(mask_shp, buffer_mask)
       
        if out_tif == None:
            out_tif = raster.replace('.tif','_crop.tif')
        
        with rxy.open_rasterio(raster, masked = True, chunks = True) as ds:
            clipped = ds.rio.clip(mask_shp.geometry.apply(mapping), mask_shp.crs, drop=False, invert=False)
            clipped.rio.to_raster(out_tif)
        
        if remove:
            os.remove(raster)             
        return out_tif


    def ptsTime2Raster(self, out_name, var_list = None, outputBounds = None, outCRS = 'WGS84', mask_shp = None, buffer_mask = 0):
        
        '''
        This method convert points time series in rasters data series by linear interpolation.
        
        input:
            :param out_name  = output name base name
                If out_name = 'D:/output/dir/path/raster_name.tif',this implies that the output 
                raster file names will have the following format: 
                'D:/output/dir/path/raster_name_time_stamp_variable_name.tif')
                
            :param var_list (optional) = variables list to convert in rasters files.
                If not informed, all available variables in dataset will be converted.
                
           
            :param outputBounds (optional) = is set as a list instance with the following format: 
                
                [upperLeft Longitude, upperLeft Latitude, lowerRight Longitude, lowerRight Latitude]
                
                This defines the interpolation area in a rectangle delimited by its upper left point 
                coordinates and lower right point coordinates.
                
                default value: the rectangle area is defined by upper left and lower right dataset 
                coordinates increased by 1%.
                
                
            :param outCRS (optional) = output Coordinate Reference System.
                                        default value: WGS84.
            
            :param mask_shp (optional) = If informed, the interpolation raster is croped to shapefile boundaries.
            
            :param buffer_mask (optional) = buffer in shapefile mask area ( in %percentage).
                                            It is only used if the mask shapefile is inserted.
                                            default value: 0%.
        
        output:
            return: multiples raster files ( format .tif)
        
        '''
        
        if not var_list:
            var_list = self._data.columns.drop(['latitude', 'longitude', 'geometry'], errors = 'ignore')
                
        if not outputBounds:
            geom = self._data.geometry.drop_duplicates()
            x1 = geom.x.max()
            x2 = geom.x.min()
            y1 = geom.y.max()
            y2 = geom.y.min()
            
            if x1*x2>=0:
                if abs(x1)<abs(x2):
                    aux = x1
                    x1 = x2
                    x2 = aux
            if y1*y2>=0:      
                if abs(y1)<abs(y2):
                    aux = y1
                    y1 = y2
                    y2 = aux
 
#            outputBounds = [x*1.02 for x in [geom.x.max(), geom.y.min(), geom.x.min(), geom.y.max()]]
            outputBounds = [x1, y2, x2, y1]

            
        if '.tif' in out_name:
            out_name = out_name.replace('.tif','')
        
        if mask_shp:
            mask_shp = gpd.GeoDataFrame.from_file(mask_shp)
            
        out_dir = '\\'.join(out_name.split('\\')[:-1])
        timeSteps = self._data.sort_index().index.drop_duplicates()
        
        for time in timeSteps:
            timeName = str(time).replace(' ','_').replace('-','_').replace(':','_')
        
            for varName in var_list:
                vrt_fn = os.path.join(out_dir,varName+'Vrt.vrt')
                lyr_name = varName
                out_tif = '_'.join([out_name,varName,timeName,'.tif'])
                tempPath = os.path.join(out_dir, varName +'.csv')
                self._data[[varName,'latitude','longitude']].loc[time].to_csv(tempPath,header = True, index = False)
                
                with open(vrt_fn, 'w') as fn_vrt:
                    fn_vrt.write('<OGRVRTDataSource>\n')
                    fn_vrt.write('\t<OGRVRTLayer name="%s">\n' % lyr_name)
                    fn_vrt.write('\t\t<SrcDataSource>%s</SrcDataSource>\n' % tempPath)
                    fn_vrt.write('\t\t<SrcLayer>%s</SrcLayer>\n' % lyr_name)
                    fn_vrt.write('\t\t<GeometryType>wkbPoint</GeometryType>\n')
                    fn_vrt.write('\t\t<GeometryField encoding="PointFromColumns" x="longitude" y="latitude" z="%s"/>\n'  %varName)
                    fn_vrt.write('\t</OGRVRTLayer>\n')
                    fn_vrt.write('</OGRVRTDataSource>\n')
                
                gridOp = gdal.GridOptions(format = 'Gtiff', outputBounds = outputBounds, algorithm = 'linear:radius=0.0:nodata = -9999', outputSRS = outCRS)
                
                if isinstance(mask_shp,gpd.GeoDataFrame):
                    temp_tif = out_name + '_' + varName +'.tif'
                    gdal.Grid(temp_tif, vrt_fn, options = gridOp)
                    self._cropRst(temp_tif, mask_shp, out_tif, remove = True, buffer_mask = buffer_mask)
                else:
                    gdal.Grid(out_tif, vrt_fn, options = gridOp)
                
                os.remove(tempPath)
                os.remove(vrt_fn)
        return 


################################
if __name__ == "__main__":
    
    ### EXAMPLE 1: MANAGE NetCDF SPATIAL DATASET:
    
    # Initialize a NetCDFSpatial instance with the file path:
    file = r'D:\\OneDrive\\Arquivos OLD PC\\Cafofo\\Balanco de Calor\\ERA-Interim\\Mangueira2000_01_01to2018_01_01.nc'
    file = r'D:\OneDrive\Arquivos OLD PC\Cafofo\interim_1979-01-01to2019-09-31_lakes.nc'
    dataset = NetCDFSpatial(file)
    
    # We can visualize the reader data and variables in file:
    dataset.data
    dataset.variables
    
    # We can choice work with only interest variables by:
    dataset.var_choice()
    dataset.variables
    # If you already know the variables names, can insert through a list by:
    dataset.var_choice(variables = ['sp', 'u10','v10','t2m'])
    dataset.variables
    # This method filter our dataset for only the specific variables.
    
    # We can apply a temporal filter in dataset through by:
    #Insert dates on format Month/Day/Year
    StartDate = "01/21/2000"
    EndDate   =  "01/21/2000"
    dataset.temporalFilter(StartDate,EndDate)
    
    #It is also possible to perform a spatial filter by two methods.
    #The first method is to inform the coordinates of a point of interest, where a search will be carried out in the dataset looking for the data closest to that point and extracting its data series:      
    interest_point = [-50.75, -33.55]
    dataset.pointFilter(interest_point)
    dataset.data
    # The dataset returned by this method is just more than the time series of the variables chosen for the point found.
    
    # The second method performs the spatial filter using a search radius starting from the center of an informed shapefile: 
    # Read dataset again and apply a temporal filter:
    dataset = NetCDFSpatial(file)
    dataset.temporalFilter(StartDate,EndDate)
    
    # Shapefile path:
    mask_path = r"Example/Example Files/shapefile/polimirim.shp"
    
    # Apply the filter:
    dataset.centroidFilter(mask_path, radius=1) # Note that radius unit is dependent from crs, in this case radius = 1 is equivalent 1 degree
    
    # Some data sets have a coarse grid, this implies that in smaller scale study areas, a filter using only the shapefile area could return a null result.
    
    # As we are working with time series, we may wish to perform a time resample. We can do this by:
    dataset.resampleDS(freq='1D', method='ffill', order=1)
    # We can inform resample frequency by 'freq' and resample method by 'method' arguments.
    # The frequency argument can be a hourly, daily or monthly fraction to a specific temporal resample (e.g. freq = '0.5H' makes a resample to each 30 minutes)
    # More information on the resample methods available, see help(NetCDFSpatial.resample).
    
    # We can export dataset to a .csv format:
    interest_point = [-50.75, -33.55]
    a = dataset
    a.pointFilter(interest_point)
    dataset.pts2csv('Examples/Example Files/dataset-to-csv.csv')
     
    # Besides, we can export the data to .tiff format:
    # The data is interpolate on rectangular raster with bounds definite by user through outputBounds argument. The default value
    # catch the limit points in dataset and create a interpolation area.
    # outputBounds = [Longitude-Left Upper Point, Latitude-Left Upper Point, Longitude-Right Lower Point, Latitude-Right Lower Point]
    
    # We can inform a shapefile mask to crop the resulting .tif file through 'mask_shp' argument. On this mask, we can also add a buffer in your area through buffer_mask (percentage) argument.
    
    dataset.ptsTime2Raster('Examples/Example Files/raster-output/raster_base_name', var_list = ['sp','t2m'], outputBounds = None, outCRS = 'WGS84', mask_shp = mask_path, buffer_mask = 0)
    

##########################################################


class RasterTemporalSeries():
       
    def __init__(self, datadir, temporal_standard = None, crs ={'init': 'epsg:4326'} ):    
        '''
        class constructor
        '''
        self._data = self._readFiles(self, datadir)
        
        ds   = xr.open_dataset(datapath)
        self._data = ds.to_dataframe().reset_index() 
        if self._data.select_dtypes(include=[np.datetime64]).keys()[0]:
            self._data = self._data.set_index(self._data.select_dtypes(include=[np.datetime64]).keys()[0])
        
        self._data['geometry'] =  self._getgeom()
        self._crs  = crs
        
        self._data = self._df2geodf(self._data['geometry'], crs = self._crs)
        self._variables = list(self._data.columns.drop(['geometry']))
        
        
    @property
    def data(self):
        return self._data
    
    
# Read files: extract time, mean, median, std, max, min, path
# Crop all rasters
# Mean Raster from periods (resample - temporal upscale) by sum, min or mean
# Silhoutte Analysis
# Kmeans classifier

def _getstatistc(self, file):
    with rasterio.open(file) as src:
        data = src.read(1)
#        mean = data.mean()
#        median = data.median()
#        maxv = data.max()
#        minv = data.min()
#        stdv = data.std()
        return [data.mean(), data.median(), data.max(), data.min(), data.std(), file]

def _readFiles(self, datadir):
    statistics = []
    for file in os.listdir(datadir):
        file = os.path.join(datadir,file)
        statistics.append(_getstatistc(self, file)) 
    return pd.DataFrame(statistics, columns = ['Mean', 'Median', 'STD', 'Max', 'Min', 'FilePath'])
        

#def _readFiles(self, datadir):
#    means = maxv = minv = stdv = filep = []
#    for file in os.listdir(datadir):
#        file = os.path.join(datadir,file)
#        with rasterio.open(file) as src:
#            data = src.read(1)
#            means.append(data.mean())
#            median.append(data.median())
#            maxv.append(data.max())
#            minv.appned(data.min())
#            stdv.append(data.std())
#            filep.append(file)    
#    return pd.DataFrame(zip(means,maxv,minv,stdv,filep), columns = ['Mean', 'Median', 'STD', 'Max', 'Min', 'FilePath'])
#

## Base spatial _cropRst, _maskCheck

def cropAllRst(self, mask_shp, out_dir = None, remove = False, buffer_mask = 0):
   statistics = []
   for raster in tqdm(self._data['FilePath']):
       file = self._cropRst(temp_tif, mask_shp, out_tif, remove = remove, buffer_mask = buffer_mask)
       statistics.append(_getstatistc(self, file))   
    self._data = pd.DataFrame(statistics, columns = ['Mean', 'Median', 'STD', 'Max', 'Min', 'FilePath'])
    return




     
def MeanArray(files):
    array_list = np.array([read_file(x) for x in files])
    return np.nanmean(array_list, axis = 0)

        # Merge rasters:
        out_name_merge = []
        for r in range(0, int(len(rasters)/2)):
            
            out_name = os.path.join(out_dir,directories[r], base_name.split('_')[0]+'_'+ '_'.join([ s.zfill(2) for s in [str(dateGreg.year), str(dateGreg.day), str(dateGreg.month), str(hh), str(mm), directories[r]]]) + '.tif')
            out_name_merge.append(out_name)
            
            mergeRst = [rasterio.open(rasters[r]).read(1), rasterio.open(rasters[r+4]).read(1)]
            
            rst1 = mergeRst[0]
            Pos1 = np.argwhere(~np.isnan(rst1))
            rows, cols = zip(*Pos1)
            
            outRst = mergeRst[1]
            outRst[rows, cols] = rst1[rows, cols]
            
            with rasterio.open( out_name , "w", **out_meta) as dest:
                dest.write(outRst,1)
        
        dic['mergeFiles'] = out_name_merge
        statistics.append(dic)

def arraysMean(files):
    # Calcula valor mediorasters
    # Baseado em https://gis.stackexchange.com/questions/244376/computing-mean-of-all-rasters-in-a-directory-using-python

    if isinstance(files, np.ndarray):
        array_list = [read_file(x) for x in files]
    else:
        files = np.array(files)
        array_list = [read_file(x) for x in files]
        
    meanValues = []
    for array in tqdm(array_list):
        array = array[~np.isnan(array)]
        meanValues.append(np.mean(array)) 
    return np.array(meanValues)

def MeanArray(files):
    array_list = np.array([read_file(x) for x in files])
    return np.nanmean(array_list, axis = 0)



class BaseSpatial:
    
    def _maskCheck(self, mask_shape, buffer_mask):
        if not isinstance(mask_shape,gpd.GeoDataFrame):
            mask_shape = gpd.GeoDataFrame.from_file(mask_shape)
            mask_shape.geometry[0] = mask_shape.geometry[0].buffer(buffer_mask)
        return mask_shape

    def _cropRst(self,raster, mask_shp, out_tif = None, remove = False, buffer_mask = 0):
        
        mask_shp = self._maskCheck(mask_shp, buffer_mask)
       
        if out_tif == None:
            out_tif = raster.replace('.tif','_crop.tif')
        
        with rxy.open_rasterio(raster, masked = True, chunks = True) as ds:
            clipped = ds.rio.clip(mask_shp.geometry.apply(mapping), mask_shp.crs, drop=False, invert=False)
            clipped.rio.to_raster(out_tif)
        
        if remove:
            os.remove(raster)             
        return out_tif


#########Media Rasters ############

chd = os.chdir("D:\OneDrive\Mestrado\Plano Dissertacao\Temperatura")
df = pd.read_pickle('VariblesTemp02.pkl')

out_dir = "D:\OneDrive\Mestrado\Plano Dissertacao\Temperatura\Modis\Lagoas"

df['DateTime'] = pd.to_datetime(df['Date'].astype(str) + ' '+df['TimeObs_Mirim'].astype(str))
df = df.set_index('DateTime')
