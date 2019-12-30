import ee

def merge_collections_std_bandnames_collection1tier1():
    """merge landsat 5, 7, 8 collection 1 tier 1 imageCollections and standardize band names
    """
    ## standardize band names
    bn8 = ['B2', 'B3', 'B4', 'B6', 'BQA', 'B5']
    bn7 = ['B1', 'B2', 'B3', 'B5', 'BQA', 'B4']
    bn5 = ['B1', 'B2', 'B3', 'B5', 'BQA', 'B4']
    bns = ['Blue', 'Green', 'Red', 'Swir1', 'BQA', 'Nir']

    # create a merged collection from landsat 5, 7, and 8
    ls5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_TOA").select(bn5, bns)

    ls7 = (ee.ImageCollection("LANDSAT/LE07/C01/T1_RT_TOA")
           .filterDate('1999-01-01', '2003-01-01')
           .select(bn7, bns))

    ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT_TOA").select(bn8, bns)

    merged = ee.ImageCollection(ls5.merge(ls7).merge(ls8))

    return(merged)

def Unpack(bitBand, startingBit, bitWidth):
    # unpacking bit bands
    # see: https://groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
    return (ee.Image(bitBand)
            .rightShift(startingBit)
            .bitwiseAnd(ee.Number(2).pow(ee.Number(bitWidth)).subtract(ee.Number(1)).int()))

bitAffected = {
    'Cloud': [4, 1],
    'CloudConfidence': [5, 2],
    'CloudShadowConfidence': [7, 2],
    'SnowIceConfidence': [9, 2]
}

def UnpackAll(bitBand, bitInfo):
    unpackedImage = ee.Image.cat([Unpack(bitBand, bitInfo[key][0], bitInfo[key][1]).rename([key]) for key in bitInfo])
    return unpackedImage

def genAffectedArea(bitInfo):
    def affectedArea(quality_band):
        bqa = UnpackAll(quality_band, bitInfo)
        return bqa
    return affectedArea

affectedArea = genAffectedArea(bitAffected)

def bufferGrwl(f):
    """buffer GRWL by 15 meters to extract information along the centerline"""
    return f.buffer(15)

def addAffectedPercentage(image):
    grwl = ee.FeatureCollection("users/eeProject/GRWL_summaryStats")
    # grwl = ee.FeatureCollection("users/eeProject/grwl")
    aoi = image.select('BQA').geometry()

    river_geometry = (grwl
                      .filterBounds(aoi)
                      .filterMetadata('width_mean', 'greater_than', 60)
                      .filterMetadata('lakeFlag', 'equals', 0)
                      .map(bufferGrwl)
                      .union()
                      .geometry()
                      .intersection(aoi))

    qa = affectedArea(ee.Image(image.select(['BQA'])))
    cloud = qa.select(['Cloud'])
    cloudConfidence = qa.select(['CloudConfidence'])
    cloudShadowConfidence = qa.select(['CloudShadowConfidence'])
    snowIceConfidence = qa.select(['SnowIceConfidence'])

    mndwi = image.normalizedDifference(['Green', 'Swir1'])
    clear = cloud.Or(cloudShadowConfidence.eq(3)).Or(snowIceConfidence.eq(3)).Not()

    # add grwl river mask
    grwl_river_mask = ee.Image("users/eeProject/grwl_river_mask")
    resolution = 30
    shift_units = 2.5
    disp = ee.Image(shift_units * resolution).addBands(ee.Image(-shift_units * resolution))
    grwl_river_mask_shifted = grwl_river_mask.displace(disp, None)

    areaBand = ee.Image.pixelArea()
    affected = (areaBand.mask(cloud).rename(['cloud'])
                .addBands(areaBand.mask(cloudConfidence.eq(0)).rename(['cc0']))
                .addBands(areaBand.mask(cloudConfidence.eq(1)).rename(['cc1']))
                .addBands(areaBand.mask(cloudConfidence.eq(2)).rename(['cc2']))
                .addBands(areaBand.mask(cloudConfidence.eq(3)).rename(['cc3']))
                .addBands(areaBand.mask(snowIceConfidence.eq(0)).rename(['sic0']))
                .addBands(areaBand.mask(snowIceConfidence.eq(1)).rename(['sic1']))
                .addBands(areaBand.mask(snowIceConfidence.eq(2)).rename(['sic2']))
                .addBands(areaBand.mask(snowIceConfidence.eq(3)).rename(['sic3']))
                .addBands(areaBand.mask(cloudShadowConfidence.eq(0)).rename(['csc0']))
                .addBands(areaBand.mask(cloudShadowConfidence.eq(1)).rename(['csc1']))
                .addBands(areaBand.mask(cloudShadowConfidence.eq(2)).rename(['csc2']))
                .addBands(areaBand.mask(cloudShadowConfidence.eq(3)).rename(['csc3']))
                .addBands(areaBand.rename(['area']))
                .addBands(clear.rename(['n_clear_pixels']))
                .addBands(grwl_river_mask_shifted.rename(['river_type'])))

    affected = (affected.updateMask(grwl_river_mask_shifted.gt(0))
                .updateMask(mndwi.gte(0)))

    condition = (affected.reduceRegion(
        reducer = ee.Reducer.sum(),
        crs = qa.select(['Cloud']).projection(),
        scale = 30,
        geometry = river_geometry,
        bestEffort = False,
        tileScale = 1))

    condition = condition.set('LANDSAT_SCENE_ID', image.get('LANDSAT_SCENE_ID'))

    return ee.Feature(None, condition)

def calculate_ice_fraction_gen(h_occ):
    jrc = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
    occ_mask = jrc.select(['occurrence']).gte(h_occ)
    dem = ee.Image("users/eeProject/MERIT").select(['elevation'])
    grwlf = ee.FeatureCollection("users/eeProject/grwl")
    grwls = ee.FeatureCollection("users/eeProject/GRWL_summaryStats")

    era5t2m1 = ee.ImageCollection('users/yxiao/ERA5_DAILY_T2M')
    era5t2m2 = ee.ImageCollection('users/yx/ERA5_DAILY_T2M')
    era5t2m = era5t2m1.merge(era5t2m2)

    def CalcPre30T2mMean(i):
        endDate = ee.Date(i.get('system:time_start'))
        startDate = endDate.advance(-30, 'day')
        mean30T2m = (era5t2m
        .filterDate(startDate, endDate)
        .mean()
        .add(-273.15))

        return(ee.Image(mean30T2m).rename(['pre30T2mMean']))

    grwlfil = ee.FeatureCollection('users/eeProject/GRWL_JRC_occGte90_widthGte90_chn1_lake0')

    def cif(image):

        bqa = ee.Image(image.select(['BQA']))
        aoi = bqa.geometry()

        river_geometry = (grwlfil
        .filterBounds(aoi)
        .union()
        .geometry())

        qa = affectedArea(bqa)
        crs = qa.select(['Cloud']).projection()
        cloud = qa.select(['Cloud'])
        cloudShadow = qa.select(['CloudShadowConfidence']).eq(3)
        snowIce = qa.select(['SnowIceConfidence']).eq(3)

        base = bqa.gte(0)

        # // hill shadow
        hillShadow = (ee.Terrain.hillShadow(dem, ee.Number(image.get('SUN_AZIMUTH')),
        ee.Number(90).subtract(image.get('SUN_ELEVATION')), 500).rename(['hillShadow']))

        affected = cloud.Or(cloudShadow)
        notAffected = affected.Not();

        # imageHsv = (image.select(['Red', 'Green', 'Blue']).rgbToHsv().multiply(360))

        input = (base.rename(['base'])
        .addBands(affected.mask(affected).rename(['cloud_and_shadow']))
        .addBands(snowIce.mask(snowIce).rename(['snow_ice']))
        .addBands(hillShadow)
        .addBands(CalcPre30T2mMean(image)))
        # .addBands(imageHsv.mask(snowIce).rename(['snow_ice_h', 'snow_ice_s', 'snow_ice_v']))
        # .addBands(imageHsv.mask(affected.Or(snowIce).Not()).rename(['water_h', 'water_s', 'water_v'])))

        condition = (input.reduceRegion(
        reducer = ee.Reducer.count().combine(ee.Reducer.mean(), None, True),
        crs = crs,
        scale = 30,
        geometry = river_geometry,
        bestEffort = False,
        maxPixels = 300000000,
        tileScale = 1))

        condition = (condition.set('LANDSAT_SCENE_ID', image.get('LANDSAT_SCENE_ID')))

        return ee.Feature(None, condition).set('cloudscore', image.get('CLOUD_COVER'))

    return cif

# def export_per_pixel_ice_status_gen():
#     # jrc = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
#     # occ_mask = jrc.select(['occurrence']).gte(h_occ)
#     alos = ee.Image("JAXA/ALOS/AW3D30_V1_1")
#     dem = alos.select(['MED'])
#     gridmet = ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")
#     ncep = ee.ImageCollection("NCEP_RE/surface_temp")
#
#     subGRWL = ee.FeatureCollection("users/eeProject/sampled_1percent_grwl")
#
#     def cifp(image):
#
#         bqa = ee.Image(image.select(['BQA']))
#         aoi = bqa.geometry()
#
#         river_geometry = (subGRWL.filterBounds(aoi))
#
#         timage = ee.Date(image.get('system:time_start'))
#
#         gridmet_temp_30 = (gridmet.select(['tmmn', 'tmmx'])
#         .filterDate(timage.advance(-29, 'day'), timage.advance(1, 'day')).reduce(ee.Reducer.mean())
#         .rename(['tmmn_30_mean', 'tmmx_30_mean'])
#         .unmask(-999))
#
#         gridmet_temp_1515 = (gridmet.select(['tmmn', 'tmmx'])
#         .filterDate(timage.advance(-15, 'day'), timage.advance(16, 'day')).reduce(ee.Reducer.mean())
#         .rename(['tmmn_1515_mean', 'tmmx_1515_mean'])
#         .unmask(-999))
#
#         ncep_temp_30 = (ncep.select(['air']).filterDate(timage.advance(-29, 'day'), timage.advance(1, 'day')).reduce(ee.Reducer.mean())).rename('ncep30mean')
#
#         ncep_temp_15 = (ncep.select(['air']).filterDate(timage.advance(-14, 'day'), timage.advance(1, 'day')).reduce(ee.Reducer.mean())).rename('ncep15mean')
#
#         ncep_temp_1515 = (ncep.select(['air']).filterDate(timage.advance(-15, 'day'), timage.advance(16, 'day')).reduce(ee.Reducer.mean())).rename('ncep1515mean')
#
#         # // var occ_mask = jrc.select(['occurrence']).gte(h_occ);
#         qa = affectedArea(ee.Image(image.select(['BQA'])))
#         crs = qa.select(['Cloud']).projection()
#         cloud = qa.select(['Cloud'])
#         cloudConfidence = qa.select(['CloudConfidence'])
#         cloudShadow = qa.select(['CloudShadowConfidence']).eq(3)
#         snowIce = qa.select(['SnowIceConfidence']).eq(3)
#
#         base = image.select(['BQA']).gte(0)
#
#         # // hill shade
#         hillShade = (ee.Terrain.hillshade(dem, ee.Number(image.get('SUN_AZIMUTH')).add(360), image.get('SUN_ELEVATION')).rename(['hillShade']))
#
#         affected = cloud.Or(cloudShadow)
#         notAffected = affected.Not()
#
#         imageHsv = image.select(['Red', 'Green', 'Blue']).rgbToHsv().multiply(360)
#
#         input = (base.rename(['base'])
#         .addBands(ee.Image.pixelLonLat())
#         .addBands(affected.rename(['cloud_and_shadow']))
#         .addBands(snowIce.rename(['snow_ice']))
#         .addBands(hillShade)
#         .addBands(imageHsv.rename(['hue', 'satuation', 'value']))
#         .addBands(gridmet_temp_30)
#         .addBands(gridmet_temp_1515)
#         .addBands(ncep_temp_30)
#         .addBands(ncep_temp_1515)
#         .addBands(ncep_temp_15))
#
#         condition = (input.sampleRegions(
#         collection = river_geometry,
#         scale = 30,
#         projection = 'EPSG:4326'))
#
#         def addSceneId(f):
#             return f.set('LANDSAT_SCENE_ID', image.get('LANDSAT_SCENE_ID'))
#
#         condition = condition.map(addSceneId)
#
#         return condition
#
#     return cifp

def maximum_no_of_tasks(MaxNActive, waitingPeriod):
	"""maintain a maximum number of active tasks
	"""
	import time
	import ee
	ee.Initialize()

	time.sleep(2)
	## initialize submitting jobs
	ts = list(ee.batch.Task.list())

	NActive = 0
	for task in ts:
		if ('RUNNING' in str(task) or 'READY' in str(task)):
			NActive += 1
	## wait if the number of current active tasks reach the maximum number
	## defined in MaxNActive
	while (NActive >= MaxNActive):
		time.sleep(waitingPeriod) # if reach or over maximum no. of active tasks, wait for 2min and check again
		ts = list(ee.batch.Task.list())
		NActive = 0
		for task in ts:
			if ('RUNNING' in str(task) or 'READY' in str(task)):
				NActive += 1
	return()
