var occ = ee.Image("JRC/GSW1_0/GlobalSurfaceWater"),
    grwl = ee.FeatureCollection("users/eeProject/grwl"),
    mergedStations = ee.FeatureCollection("users/eeProject/fmaskValidation/merged_val_stations");

var merge_collections_std_bandnames_collection1tier1 = function() {
  // merge landsat 5, 7, 8 collection 1 tier 1 imageCollections and standardize band names
  // refer and define standard band names
  var bn8 = ['B2', 'B3', 'B4', 'B6', 'BQA', 'B5'];
  var bn7 = ['B1', 'B2', 'B3', 'B5', 'BQA', 'B4'];
  var bn5 = ['B1', 'B2', 'B3', 'B5', 'BQA', 'B4'];
  var bns = ['Blue', 'Green', 'Red', 'Swir1', 'BQA', 'Nir'];

  // create a merged collection from landsat 5, 7, and 8
  var ls5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_TOA").select(bn5, bns);
  var ls7 = (ee.ImageCollection("LANDSAT/LE07/C01/T1_RT_TOA")
  .filterDate('1999-04-15', '2003-05-30') // exclude LE7 images that affected by failure of Scan Line Corrector (SLC)
  .select(bn7, bns));
  var ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT_TOA").select(bn8, bns);
  var merged = ls5.merge(ls7).merge(ls8);

  return(merged);
};

/* functions to add quality band (Fmask) */
var Unpack = function(qualityBand, startingBit, bitWidth) {
  // unpacking bit information from the quality band given the bit position and width
  // see: https://groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
  return(qualityBand
  .rightShift(startingBit)
  .bitwiseAnd(ee.Number(2).pow(bitWidth).subtract(1).int()));
};
var UnpackAll = function(bitBand) {
  // apply Unpack function for multiple pixel qualities
  var bitInfoTOA = {
    'Cloud': [4, 1],
    'CloudShadowConfidence': [7, 2],
    'SnowIceConfidence': [9, 2]
  };
  var unpackedImage = ee.Image();
  for (var key in bitInfoTOA) {
    unpackedImage = ee.Image.cat(
      unpackedImage, Unpack(bitBand, bitInfoTOA[key][0], bitInfoTOA[key][1])
      .rename(key));
  }
  return(unpackedImage.select(Object.keys(bitInfoTOA)));
};
var ClassifyWaterFmask = function(image) {
  // apply water classification from Fmask (see Zhu and Woodcock, 2012)
  var ndvi = image.normalizedDifference(['Nir', 'Red']);
  var nir = image.select(['Nir']);
  var fwater = ndvi.lt(0.01).and(nir.lt(0.11)).or(ndvi.lt(0.1).and(nir.lt(0.05)));

  return(fwater);
};
var AddFmask = function(image) {
  // add fmask as a separate band to the input image
  var temp = UnpackAll(image.select(['BQA']));
  var fwater = ClassifyWaterFmask(image);

  // construct the fmask
  var fmask = (fwater.rename(['fmask'])
  .where(temp.select(['SnowIceConfidence']).eq(3), ee.Image(3))
  .where(temp.select(['CloudShadowConfidence']).eq(3), ee.Image(2))
  .where(temp.select(['Cloud']), ee.Image(4)))
  .mask(temp.select(['Cloud']).gte(0)); // mask the fmask so that it has the same footprint as the quality (BQA) band

  return(image.addBands(fmask));
};

/* functions to calculate hill shadow */
var CalcHillShadow = function(image) {
  var dem = ee.Image("users/eeProject/MERIT").clip(image.geometry().buffer(9000).bounds());
  var SOLAR_AZIMUTH_ANGLE = ee.Number(image.get('SUN_AZIMUTH'));
  var SOLAR_ZENITH_ANGLE = ee.Number(90).subtract(image.get('SUN_ELEVATION'));
  return(ee.Terrain.hillShadow(dem, SOLAR_AZIMUTH_ANGLE, SOLAR_ZENITH_ANGLE, 100, true).rename(['hillshadow']));
};


var ls = merge_collections_std_bandnames_collection1tier1().filterMetadata('CLOUD_COVER', 'less_than', 25);

print(mergedStations.first());
print(ls.first());

var occ90 = occ.select('occurrence').gte(90);

var grwlFil = grwl
      .filterMetadata('width_m', 'greater_than', 90)
      .filterMetadata('lakeFlag', 'equals', 0)
      .filterMetadata('nchannels', 'equals', 1);

// implment

var riverIce = mergedStations.map(function(f) {

  var getFmaskStatus = function(i) {

    // add fmask and mask that with 90% occurrence of water
    var fmask = AddFmask(i).select('fmask').mask(occ90);

    // add hill shadow
    fmask = fmask.addBands(CalcHillShadow(i));

    var fmaskSampled = fmask.reduceRegion({
      reducer: ee.Reducer.toList(),
      geometry: grwlFil
      .filterBounds(f.geometry().buffer(1500))
      .union(),
      scale: 30});

    return(f
    .setGeometry(null)
    .set(fmaskSampled)
    .copyProperties(i, ['CLOUD_COVER', 'LANDSAT_SCENE_ID', 'system:time_start'])
    );
  };

  // only calculate images (1) from the years when there were in situ observations (2) contains the guages
  var fil = ee.Filter.contains({
    leftField: '.geo',
    rightValue: f.geometry().buffer(1500)
  });

  var ice_fraction = ls
  .filterDate(ee.String(ee.Number(f.get('date1')).int()).cat('-01-01'), ee.String(ee.Number(f.get('date2')).int().add(1)).cat('-12-31'))
  .filter(fil)
  .map(getFmaskStatus)

  return(ice_fraction);

}).flatten();

print(riverIce.toList(3));

Map.addLayer(mergedStations, {color: 'red'});

Export.table.toDrive({
  collection: riverIce,
  description: 'river_ice_validation',
  folder: 'Validating-Fmask-snow-ice-classification',
  fileNamePrefix: 'river_ice_validation',
  fileFormat: 'csv'});
