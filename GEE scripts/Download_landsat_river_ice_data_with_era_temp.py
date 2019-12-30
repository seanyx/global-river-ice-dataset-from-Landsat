# // This script samples pixel types along river centerlines in GRWL and export ice covered and total centerline pixel areas for each landsat image
# // by: Xiao Yang
# // 03/17/2018

import ee
import argparse
import math
from functions import *

parser = argparse.ArgumentParser(prog = 'Download_landsat_river_ice_data.py', description = "")

parser.add_argument('-sy', '--start_year', help = 'starting (inclusive) year of the export data', type = int, default = 2001)
parser.add_argument('-sm', '--start_month', help = 'starting (inclusive) month of the export data', type = int, default = 8)
parser.add_argument('-ey', '--end_year', help = 'end (inclusive) year of the export data', type = int, default = 2001)
parser.add_argument('-em', '--end_month', help = 'end (inclusive) month of the export data', type = int, default = 9)

parser.add_argument('-minc', '--minimal_cloud_score', help = 'minimal (inclusive) cloud score to include', type = int, default = 0)
parser.add_argument('-maxc', '--maximum_cloud_score', help = 'maximum (inclusive) cloud score to include', type = int, default = 25)

parser.add_argument('-o', '--output_dir', help = 'output directory in the Google drive', type = str, default = 'default')
parser.add_argument('-n', '--diag_note', help = 'note appended to the filename', type = str, default = '')

parser.add_argument('--restart_no', help = 'index to start the task from', type = int, default = 0)

args = parser.parse_args()
start_year = args.start_year
start_month = args.start_month
end_year = args.end_year
end_month = args.end_month
minc = args.minimal_cloud_score
maxc = args.maximum_cloud_score
output_dir = args.output_dir
rs = args.restart_no
diagnote = args.diag_note

ee.Initialize()

ls_data = merge_collections_std_bandnames_collection1tier1()
aoi = ee.Geometry.LineString([[180, 90], [-180, -90]])
fil = ee.Filter.geometry(geometry = aoi).Not()

n = (end_year - start_year) * 12 + (end_month - start_month) + start_month + 1

if output_dir == 'default':
    output_dir = 'cloudiness_landsat578_collection1_tier1_' + 'cs' + str(minc) + '-' + str(maxc) + '_' + str(start_year) + '-' + str(start_month) + '_' + str(end_year) + '-' + str(end_month)

for i in range(start_month + rs, n):

    yearOffset = math.floor((i - 1) / 12)
    year = int(start_year + yearOffset)
    month = int(i - yearOffset * 12)

    year_offset_next = math.floor(i / 12)
    year_next = int(start_year + year_offset_next)
    month_next = int(i + 1 - year_offset_next * 12)

    t_start = str(year) + '-' + format(month, '02d') + '-01'
    t_end = str(year_next) + '-' + format(month_next, '02d') + '-01'

    calculate_ice_fraction = calculate_ice_fraction_gen(0.9)

    result = ee.FeatureCollection(ls_data
    .filter(fil)
    .filterDate(t_start, t_end)
    .filterMetadata('CLOUD_COVER', 'not_greater_than', maxc)
    .filterMetadata('CLOUD_COVER', 'not_less_than', minc)
    # .limit(100)
    .map(calculate_ice_fraction))

    ## add header column to avoid missing columns in the exported file
    columns = ee.List([
    'LANDSAT_SCENE_ID',
    'base_count',
    'cloud_and_shadow_count',
    'cloud_and_shadow_mean',
    'snow_ice_count',
    'snow_ice_mean',
    'base_mean',
    'hillShadow_count',
    'hillShadow_mean',
    'pre30T2mMean_mean',
    'pre30T2mMean_count',
    'cloudscore'
    ])

    header = ee.Feature(None, ee.Dictionary.fromLists(columns, ee.List.repeat(-9999, columns.size())))
    result = ee.FeatureCollection([header]).merge(result)

    fn = 'river_ice_temp_cs_' + format(minc, '02d') + '-' + format(maxc, '02d') + '_' + t_start + '_' + diagnote

    task = (ee.batch.Export.table.toDrive(
            collection = result,
            description = fn,
            folder = output_dir,
            fileNamePrefix = fn,
            fileFormat = 'CSV'))

    task.start()

    maximum_no_of_tasks(10, 120)

    print('task', i - start_month + 1, 'of', n - start_month, ', (', t_start, t_end, ')', 'has submitted.')
