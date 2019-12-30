var cmip5 = ee.ImageCollection("NASA/NEX-GDDP"),
    land = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017");

var changeBetweenTwoPeriod = function(model, scenario, aoi, period1, period2, spatialReducer) {
  var cmip5 = ee.ImageCollection("NASA/NEX-GDDP");
  var dailyTas1 = cmip5
  .filterMetadata('model', 'equals', model)
  .filterMetadata('scenario', 'equals', scenario)
  .filterDate(period1[0], period1[1])
  .select(['tasmin', 'tasmax'])
  .reduce(ee.Reducer.mean());

  var dailyTas2 = cmip5
  .filterMetadata('model', 'equals', model)
  .filterMetadata('scenario', 'equals', scenario)
  .filterDate(period2[0], period2[1])
  .select(['tasmin', 'tasmax'])
  .reduce(ee.Reducer.mean());

  var diff = dailyTas2.subtract(dailyTas1);

  var diffDailyMean = diff.select(['tasmin_mean'])
  .add(diff.select(['tasmax_mean']))
  .divide(2)
  .rename('dailyMeanDiff');

  var result = diffDailyMean.reduceRegion({
    reducer: spatialReducer,
    geometry: aoi,
    scale: 300000,
    crs: 'EPSG:4326',
    bestEffort: true
  });

  return(ee.Feature(null, result).setMulti({
    model: model,
    scenario: scenario,
    basePeriod: period1,
    targetPeriod: period2,
    aoi: aoi,
    diffImage: diffDailyMean
  }));
};



// simple test
var period1 = ['2006-01-01', '2007-01-01'];
var period2 = ['2069-12-31', '2070-12-31'];
var model = 'CESM1-BGC';
var scenario = 'rcp85';
var aoi = ee.Geometry.Rectangle([-179, -90, 179, 90], 'EPSG:4326', false, true).aside(Map.addLayer);


var combinedReducer = ee.Reducer.mean()
.combine(ee.Reducer.median(), null, true)
.combine(ee.Reducer.count(), null, true);

// var test = changeBetweenTwoPeriod(model, scenario, aoi, period1, period2, combinedReducer);

// print(test);

// Map.addLayer(ee.Image(test.get('diffImage')), {palette: ['blue', 'white', 'red'], min: -8, max: 8});

// Export.table.toDrive({
//   collection: ee.FeatureCollection([test]),
//   description: 'simpleTest',
//   folder: '',
//   fileNamePrefix: 'tempdiffclimatemodel',
//   fileFormat: 'csv',
//   selectors: ['dailyMeanDiff_mean', 'dailyMeanDiff_median', 'dailyMeanDiff_count', 'model', 'scenario',
//   'basePeriod', 'targetPeriod', 'aoi']
// });


// // 30 year test
// var period1 = ['2006-01-01', '2036-01-01'];
// var period2 = ['2069-12-31', '2099-12-31'];
// var model = 'CESM1-BGC';
// var scenario = 'rcp85';
// var aoi = ee.Geometry.Rectangle([-179, -90, 179, 90], 'EPSG:4326', false, true).aside(Map.addLayer);

// var combinedReducer = ee.Reducer.mean()
// .combine(ee.Reducer.median(), null, true)
// .combine(ee.Reducer.count(), null, true);

// var test = changeBetweenTwoPeriod(model, scenario, aoi, period1, period2, combinedReducer);

// // print(test);

// // Map.addLayer(ee.Image(test.get('diffImage')), {palette: ['blue', 'white', 'red'], min: -8, max: 8});

// Export.table.toDrive({
//   collection: ee.FeatureCollection([test]),
//   description: '30_year_test',
//   folder: '',
//   fileNamePrefix: 'tempdiffclimatemodel',
//   fileFormat: 'csv',
//   selectors: ['dailyMeanDiff_mean', 'dailyMeanDiff_median', 'dailyMeanDiff_count', 'model', 'scenario',
//   'basePeriod', 'targetPeriod', 'aoi']
// });


// all model comparison

var period1 = ['2006-01-01', '2036-01-01'];
var period2 = ['2069-12-31', '2099-12-31'];

// var period1 = ['2006-01-01', '2007-01-01'];
// var period2 = ['2069-12-31', '2070-12-31'];
// var model = 'CESM1-BGC';
var scenario = 'rcp45';
var aoi = ee.Geometry.Rectangle([-179, -90, 179, 90], 'EPSG:4326', false, true).aside(Map.addLayer);

var models = ee.List(['ACCESS1-0', 'bcc-csm1-1', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MIROC5', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'NorESM1-M']);

var combinedReducer = ee.Reducer.mean()
.combine(ee.Reducer.median(), null, true)
.combine(ee.Reducer.count(), null, true);

var result = models.map(function(l) {
  var temp = changeBetweenTwoPeriod(l, scenario, aoi, period1, period2, combinedReducer);
  return(temp);
});

print(result);

// Map.addLayer(ee.Image(test.get('diffImage')), {palette: ['blue', 'white', 'red'], min: -8, max: 8});

Export.table.toDrive({
  collection: ee.FeatureCollection(result),
  description: '30_year_all_models',
  folder: '',
  fileNamePrefix: 'tempdiffclimatemodel',
  fileFormat: 'csv',
  selectors: ['dailyMeanDiff_mean', 'dailyMeanDiff_median', 'dailyMeanDiff_count', 'model', 'scenario',
  'basePeriod', 'targetPeriod', 'aoi']
});
