import os 
os.environ['USE_PYGEOS'] = '0'
import sys 
import h5py
import copy
import time
import pandas as pd 
import numpy as np 
import geopandas as gpd 
from shapely.geometry import Polygon
import warnings
warnings.filterwarnings("ignore")


# --------------------DEFINE PRESET BAND/LAYER SUBSETS ------------------------------------------ #
# Default layers to be subset and exported, see README for information on how to add additional layers


def convert_vector_fmt(vi_path, vo_path, vo_fmt):
    g = gpd.read_file(vi_path)
    if not os.path.exists(vo_path):
        g.to_file(vo_path, driver=vo_fmt)

def wget_download(url, odir, uname, pword):
    time.sleep(0.5)
    cmd  = f"wget --user={uname} --password={pword} --directory-prefix={odir} {url}"
    os.system(cmd)

def gedi_hdf5_to_vector(odir,idir,roi_file, extra_layers,g, l1bSubset,l2aSubset,l2bSubset,out_fmt = 'gpkg'):

    ROI,finalClip = get_roi_gdf(roi_file)
    inDir = make_dirpath_system_agnostic(idir)
    beamSubset = subset_beams() 
    create_directory(odir)
    
    gedi, gediName, sdsSubset, gediSDS,beams = gedi_load_metadata(g, l1bSubset,l2aSubset,l2bSubset,extra_layers,beamSubset)
    #gedi_load_metadata(g, l1bSubset,l2aSubset,l2bSubset)
    gedi_bdf = add_beams_to_gedi_base_df(beams,gediSDS, sdsSubset, gedi, ROI)
    ##@ open the short following the tutorial 2
    if gedi_bdf.shape[0] == 0:
        print(f"No intersecting shots were found between {g} and the region of interest submitted.")
        #continue
    gedi_ldf = gedi_to_layersdf(gedi,beams,gediSDS, sdsSubset,gedi_bdf)
    outDF = pd.merge(gedi_bdf, gedi_ldf, left_on='shot_number', 
                 right_on=[sn for sn in gedi_ldf.columns if sn.endswith('shot_number')][0])
    
    print(outDF.shape, gedi_bdf.shape, gedi_ldf.shape)
    del gedi_bdf, gedi_ldf   
    outDF = gpd.overlay(outDF, finalClip)
    print(outDF.shape)
    outfile = f"{odir}{g}"
    outfile = outfile.replace('.h5',f'.{out_fmt}')

    cols2rename = {'geolocation_digital_elevation_model':'tdemx',
               'geolocation_elev_highestreturn':'dsm',
               'geolocation_elev_lowestmode':'dtm',
               'Latitude':'lat','Longitude':'lon',
               'BEAM':'beam',
               'l2a_quality_flag':'qlty_l2a',
               'l2b_quality_flag':'qlty_l2b',

               }

    cols2drop = ['geolocation_shot_number', 'beam','rhog']
    outDF= outDF.drop(cols2drop, axis=1)
    outDF = outDF.rename(columns=cols2rename)
    outDF['chm'] = outDF['dsm'] - outDF['dtm'] # should it be absolute diff?
    colsinorder = [ 'tdemx', 'dsm', 'dtm', 'chm','qlty_l2a', 'qlty_l2b','shot_number','index', 'beam', 'lat', 'lon','geometry']
    outDF = outDF[colsinorder]

    print(f"{g} {os.path.basename(outfile)} saved at: {odir}")
    if out_fmt == 'gpkg':
        if not os.path.exists(outfile):
            outDF.to_file(outfile, driver='GPKG')
    elif out_fmt == 'geojson':
         if not os.path.exists(outfile):
            outDF.to_file(outfile, driver='GeoJSON')

    vo_path = outfile.replace('.geojson', '.gpkg')
    try:
        convert_vector_fmt(outfile, vo_path, 'GPKG')
    except:
        pass
    return outDF 

    print('################ gedi_hdf5_to_vector ######################')

def gedi_to_layersdf(gedi,beams,gediSDS, sdsSubset,gediDF):
    beamsDF = pd.DataFrame()
    j = 0
    for b in beams:
        ### f1
        beamDF = pd.DataFrame()
        beamSDS = [s for s in gediSDS if b in s and not any(s.endswith(d) for d in sdsSubset[0:3])]
        shot = f'{b}/shot_number'
        try:
            # set up indexes in order to retrieve SDS data only within the clipped subset from above
            mindex = min(gediDF[gediDF['BEAM'] == b]['index'])
            maxdex = max(gediDF[gediDF['BEAM'] == b]['index']) + 1
            shots = gedi[shot][mindex:maxdex]
        except ValueError:
            print(f"No intersecting shots found for {b}")
            continue

        for s in beamSDS:
            j += 1
            sName = s.split('/', 1)[-1].replace('/', '_')

            # Datasets with consistent structure as shots
            if gedi[s].shape == gedi[shot].shape:
                beamDF[sName] = gedi[s][mindex:maxdex]  # Subset by index
            
            # Datasets with a length of one 
            elif len(gedi[s][()]) == 1:
                beamDF[sName] = [gedi[s][()][0]] * len(shots) # create array of same single value
            
            # Multidimensional datasets
            elif len(gedi[s].shape) == 2 and 'surface_type' not in s: 
                allData = gedi[s][()][mindex:maxdex]
                
                # For each additional dimension, create a new output column to store those data
                for i in range(gedi[s].shape[1]):
                    step = []
                    for a in allData:
                        step.append(a[i])
                    beamDF[f"{sName}_{i}"] = step

            # Waveforms
            elif s.endswith('waveform') or s.endswith('pgap_theta_z'):
                waveform = []
                
                if s.endswith('waveform'):
                    # Use sample_count and sample_start_index to identify the location of each waveform
                    start = gedi[f'{b}/{s.split("/")[-1][:2]}_sample_start_index'][mindex:maxdex]
                    count = gedi[f'{b}/{s.split("/")[-1][:2]}_sample_count'][mindex:maxdex]
                
                # for pgap_theta_z, use rx sample start index and count to subset
                else:
                    # Use sample_count and sample_start_index to identify the location of each waveform
                    start = gedi[f'{b}/rx_sample_start_index'][mindex:maxdex]
                    count = gedi[f'{b}/rx_sample_count'][mindex:maxdex]
                wave = gedi[s][()]
                
                # in the dataframe, each waveform will be stored as a list of values
                for k in range(len(start)):
                    singleWF = wave[int(start[k] - 1): int(start[k] - 1 + count[k])]
                    waveform.append(','.join([str(q) for q in singleWF]))
                beamDF[sName] = waveform
            
            # Surface type 
            elif s.endswith('surface_type'):
                surfaces = ['land', 'ocean', 'sea_ice', 'land_ice', 'inland_water']
                allData = gedi[s][()]
                for i in range(gedi[s].shape[0]):
                    beamDF[f'{surfaces[i]}'] = allData[i][mindex:maxdex]
                del allData
            else:
                print(f"SDS: {s} not found")
            print(f"Processing {j} of {len(beamSDS) * len(beams)}: {s}")

        beamsDF = beamsDF.append(beamDF)
    del beamDF, beamSDS, beams, gedi, gediSDS, shots, sdsSubset
    return beamsDF

def gedi_open(g):
    gedi = h5py.File(g, 'r')      # Open file
    gediName = g.split('.h5')[0]  # Keep original filename
    gedi_objs = []            
    gedi.visit(gedi_objs.append)  # Retrieve list of datasets
    return gedi, gediName, gedi_objs

def gedi_load_metadata(g, l1bSubset,l2aSubset,l2bSubset,extra_layers,beamSubset):
    gedi, gediName, gedi_objs = gedi_open(g)
    gediSDS = [str(o) for o in gedi_objs if isinstance(gedi[o], h5py.Dataset)]
    sdsSubset = match_product_to_layers(g,l1bSubset,l2aSubset,l2bSubset)
    sdsSubset = add_extra_layer(extra_layers,sdsSubset)
    gediSDS = [c for c in gediSDS if any(c.endswith(d) for d in sdsSubset)]
    beams = get_all_unique_beams(gediSDS,beamSubset)
    del gedi_objs
    return gedi, gediName, sdsSubset, gediSDS,beams

# parallelize from here 
def add_beams_to_gedi_base_df(beams,gediSDS, sdsSubset, gedi, ROI): 
    gediDF = pd.DataFrame() # check this again for the interval between size []
    for b in beams:
        beamSDS = [s for s in gediSDS if b in s]
        geoDF = gedi_to_basedf(b, beamSDS, sdsSubset,gedi,ROI)
        gediDF = gediDF.append(geoDF)
        del geoDF
    gediDF = gpd.GeoDataFrame(gediDF)
    gediDF.crs = 'EPSG:4326'
    return gediDF

def gedi_to_basedf(b, beamSDS, sdsSubset,gedi,ROI):

  # Search for latitude, longitude, and shot number SDS
  lat = [l for l in beamSDS if sdsSubset[0] in l][0]  
  lon = [l for l in beamSDS if sdsSubset[1] in l][0]
  shot = f'{b}/shot_number'          

  # Open latitude, longitude, and shot number SDS
  shots = gedi[shot][()]
  lats = gedi[lat][()]
  lons = gedi[lon][()]

  # Append BEAM, shot number, latitude, longitude and an index to the GEDI dataframe
  geoDF = pd.DataFrame({'BEAM': len(shots) * [b], 
                        shot.split('/', 1)[-1].replace('/', '_'): shots,
                      'Latitude':lats, 'Longitude':lons, 
                      'index': np.arange(0, len(shots), 1)})
  geoDF = gpd.GeoDataFrame(geoDF, geometry=gpd.points_from_xy(geoDF.Longitude, geoDF.Latitude))
          
  # Clip to only include points within the user-defined bounding box
  geoDF = geoDF[geoDF['geometry'].within(ROI.envelope)] 
  return geoDF


def subset_beams(beams=None):
    if beams is not None:
        #beamSubset = args.beams.split(',')
        beamSubset = beams
        return beamSubset
    else:
        beamSubset = ['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']
        return beamSubset

def add_extra_layer(extra_layers,sdsSubset):
    if extra_layers is not None:
        [sdsSubset.append(y) for y in extra_layers]
    return sdsSubset

def get_all_unique_beams(gediSDS,beamSubset):
    beams = []
    for h in gediSDS:
        beam = h.split('/', 1)[0]
        if beam not in beams and beam in beamSubset:
            beams.append(beam)
    return beams

# Define subset of layers based on product
def match_product_to_layers(g,l1bSubset,l2aSubset,l2bSubset):
    if 'GEDI01_B' in g:
        sdsSubset = copy.deepcopy(l1bSubset)
    elif 'GEDI02_A' in g:
        sdsSubset = copy.deepcopy(l2aSubset)
    else:
        sdsSubset = copy.deepcopy(l2bSubset)
    return sdsSubset

def get_roi_gdf(ROI):
    if ROI.endswith('.geojson') or ROI.endswith('.shp') or ROI.endswith('.gpkg'):
        try:
            ROI = gpd.GeoDataFrame.from_file(ROI)
            ROI.crs = 'EPSG:4326'
            if len(ROI) > 1:
                print('Multi-feature polygon detected. Only the first feature will be used to subset the GEDI data.')
            ROI = ROI.geometry[0]
        except:
            print('error: unable to read input geojson file or the file was not found')
            sys.exit(2)
    else:
        ROI = ROI.replace("'", "")
        ROI = ROI.split(',')
        ROI = [float(r) for r in ROI]
        try:
            ROI = Polygon([(ROI[1], ROI[0]), (ROI[3], ROI[0]), (ROI[3], ROI[2]), (ROI[1], ROI[2])]) 
        except:
            print('error: unable to read input bounding box coordinates, the required format is: ul_lat,ul_lon,lr_lat,lr_lon')
            sys.exit(2)
    finalClip = gpd.GeoDataFrame(index=[0], geometry=[ROI], crs='EPSG:4326')  
    print (finalClip)
    return ROI, finalClip

def make_dirpath_system_agnostic(idir):
    if idir[-1] != '/' and idir[-1] != '\\':
        inDir = idir.strip("'").strip('"') + os.sep
    else:
        inDir = idir

    # Find input directory
    try:
        os.chdir(inDir)
    except FileNotFoundError:
        print('error: input directory (--dir) provided does not exist or was not found')
        sys.exit(2)
    return inDir

def create_directory(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    print(dirpath)