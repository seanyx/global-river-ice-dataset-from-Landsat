var cmip5 = ee.ImageCollection("NASA/NEX-GDDP"),
    grwl = ee.FeatureCollection("users/eeProject/grwl");

var model = 'MIROC-ESM'; //['CESM1-BGC', 'GFDL-ESM2M', 'MIROC-ESM']
var scenario = 'rcp85';
var date1 = '2007-12-01';
var date2 = '2099-12-31';

var dailyTas = cmip5
.filterMetadata('model', 'equals', model)
.filterMetadata('scenario', 'equals', scenario)
.filterDate(date1, date2)
.select('tasmin', 'tasmax');

print(dailyTas.first());

var calculate30DayPriorMeanTempGen = function(cmip5) {
  var temp = function(i) {
    var endDate = ee.Date(i.get('system:time_start'));
    var startDate = endDate.advance(-30, 'day');
    var mean30Tas = cmip5
    .filterDate(startDate, endDate)
    .map(function(i) {
      return(i.select(['tasmin']).add(i.select(['tasmax'])).divide(2).rename(['tasmean']));
    })
    .mean()
    .add(-273.15)
    .setMulti({
      'firstDate': startDate,
      'lastDate': endDate,
      'year': endDate.get('year')
    })
    .copyProperties(i, ['system:time_start']);

    return(mean30Tas);
  };
  return(temp);
};

var calculate30DayPriorMeanTemp = calculate30DayPriorMeanTempGen(dailyTas);

// var tempImage = calculate30DayPriorMeanTemp(ee.Image(dailyTas.toList(1, 40).get(0)));

// print(tempImage);

var appendPeriodBand = function(i) {
  var doy = ee.Number.parse(ee.Date(i.get('system:time_start')).format('D'));
  var period = ee.Algorithms.If(doy.gt(220).or(doy.lt(40)), 'freeze-up', 'breakup');
  var periodValue = ee.Algorithms.If(doy.gt(220).or(doy.lt(40)), 0, 1);
  return(ee.Image(i.setMulti({'period': period, 'periodValue': periodValue})));
};

// var tempImagePeriod = appendPeriodBand(tempImage);
// print(tempImagePeriod);

// Map.addLayer(tempImagePeriod);

var calculateRiverIce = function(i) {
  var iceFraction = ee.Image(1).subtract(ee.Image(1).divide(i.multiply(ee.Number(-0.32233).add(ee.Number(i.get('periodValue')).multiply(-0.04515))).add(-0.82245).exp().add(1)));
  return(i.addBands(iceFraction.rename(['iceFraction'])));
};

// var tempIceFraction = calculateRiverIce(tempImagePeriod).aside(print);

// Map.addLayer(tempIceFraction.select(['iceFraction']), {min: 0, max: 1});

var calculateIceFractionDuration = function(i) {
  var temp = calculate30DayPriorMeanTemp(i);
  temp = appendPeriodBand(temp);
  var iceFraction = calculateRiverIce(temp);
  var iceStatus = iceFraction.select(['iceFraction']).gte(0.5).rename(['iceStatus']);
  return(iceFraction.addBands(iceStatus));
};

// var iceOneday = calculateIceFraction(ee.Image(dailyTas.toList(1, 40).get(0)));

// print(iceOneday);
// Map.addLayer(iceOneday);

var calculateIceFractionDurationPeriod = function(collection, startDate, endDate) {
  var riverIceStats = collection
  .filterDate(startDate, endDate)
  .map(calculateIceFractionDuration)
  .mean();

  return(riverIceStats);
};

// var iceOneYear = calculateIceFractionPeriod(dailyTas, '2019-01-01', '2020-01-01');

// print(iceOneYear);

// Map.addLayer(iceOneYear);

print(grwl
  .filterMetadata('width_m', 'greater_than', 90)
  .filterMetadata('nchannels', 'equals', 1)
  .filterMetadata('lakeFlag', 'equals', 0)
  // .filterMetadata('lat', 'greater_than', 0)
  // .filterMetadata('lat', 'not_greater_than', 23.5)
  .size(), 'all npoints');

var partitionGRWL = function() {
  var grwlFil = grwl
  .filterMetadata('width_m', 'greater_than', 90)
  .filterMetadata('nchannels', 'equals', 1)
  .filterMetadata('lakeFlag', 'equals', 0);
  // .filterMetadata('lat', 'greater_than', 23.5);

  var step = 15;
  var lonStart = ee.List.sequence(-180, 180 - step, step);

  var grwlFilCollection = lonStart.map(function(f) {
    var a = grwlFil
    .filterMetadata('lon', 'not_less_than', f)
    .filterMetadata('lon', 'less_than', ee.Number(f).add(step));
    var tmp = ee.Feature(a.union()
    .first()).setMulti({'npoints': a.size(),
      'lonStart': f
    });

    return(tmp);
  });

  return(ee.FeatureCollection(grwlFilCollection));
};

var grwlFilCollection = partitionGRWL();

ee.Feature(grwlFilCollection.first()).get('npoints').aside(print)

// print(grwlFilCollection.aggregate_array('npoints'));

var weights = grwlFilCollection
.map(function(f) {return(f.setGeometry(null));});

// Export.table.toAsset(weights, 'grwlByLonBands', 'river_ice/grwlByLonBandsNpoints');
Export.table.toDrive({
  collection: weights,
  description: 'grwlByLonBands',
  folder: '',
  fileNamePrefix: 'grwlByLonBandsNpoints',
  fileFormat: 'csv',
  selectors: ['npoints', 'lonStart']
});

// Map.addLayer(ee.Feature(grwlFilCollection.first().aside(print)));

var CalculateOneYearIceMeanGen = function(geometry) {
  return(function(l) {
    var firstDate = ee.Date(ee.Number(l).format('%4d').cat('-01-01')).advance(29, 'day');
    var lastDate = ee.Date(ee.Number(l).add(1).format('%4d').cat('-01-01'));

    var tempIceMean = calculateIceFractionDurationPeriod(dailyTas, firstDate, lastDate);

    tempIceMean = tempIceMean
    .select(['iceStatus'])
    .gt(0)
    .rename(['iceAffected'])
    .addBands(tempIceMean);

    var annualTempIceGRWL = tempIceMean.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geometry,
      crs: 'EPSG:4326',
      crsTransform: [0.25, 0, 0, 0, 0.25, 0],
      tileScale: 1,
      maxPixels: 100000000000
    });

    return(ee.Feature(null, annualTempIceGRWL).setMulti({'year': l}));
  });

};

for (var i = 0; i <= 23; i++) {
  // var grwlFilGeometry = ee.Feature(grwlFil.union().first()).geometry();

  var CalculateOneYearIceMean = CalculateOneYearIceMeanGen(ee.Feature(grwlFilCollection.toList(1, i).get(0)).geometry());

  // print(ee.Feature(grwlFilCollection.toList(1, i).get(0)).get('npoints'));

  var annualMeanIce = ee.List.sequence(2008, 2099, 1)
  .map(CalculateOneYearIceMean).flatten();

  Export.table.toDrive({
    collection: ee.FeatureCollection(annualMeanIce),
    description: model + '_' + scenario + '_' + i,
    folder: 'predicted_annual_river_ice_global_0899',
    fileNamePrefix: model + '_' + scenario + '_' + i,
    fileFormat: 'csv'});
}


// calculate annual mean temperature

var calculateMeanTemp = function(i) {
      return(i.select(['tasmin']).add(i.select(['tasmax'])).divide(2).rename(['tasmean']));
};

var annualMeanTemp = function(l) {
    var firstDate = ee.Date(ee.Number(l).format('%4d').cat('-01-01')).advance(29, 'day');
    var lastDate = ee.Date(ee.Number(l).add(1).format('%4d').cat('-01-01'));

    var tempMean = calculateMeanTemp(dailyTas
    .filterDate(firstDate, lastDate)
    .mean())
    .add(-273.15);

    var coords = [[-180.0, -90.0], [180,  -90.0], [180, 90], [-180,90]];
    var geometry = ee.Geometry.Polygon(coords, 'EPSG:4326', false);

    // Map.addLayer(geometry);

    var annualMeanTemp = tempMean.reduceRegion({
      reducer: ee.Reducer.mean().combine(ee.Reducer.count(), null, true),
      geometry: geometry,
      crs: 'EPSG:4326',
      crsTransform: [0.25, 0, 0, 0, 0.25, 0],
      tileScale: 1,
      maxPixels: 100000000000,
      bestEffort: false
    });

    return(ee.Feature(null, annualMeanTemp).setMulti({'year': l}));
};

print(annualMeanTemp(2019));

var annualMeanTemp = ee.List.sequence(2008, 2099, 1)
  .map(annualMeanTemp);

Export.table.toDrive({
    collection: ee.FeatureCollection(annualMeanTemp),
    description: model + '_' + scenario + '_global_mean_temp',
    folder: 'predicted_global_mean_SAT',
    fileNamePrefix: model + '_' + scenario + '_global_mean_temp',
    fileFormat: 'csv'});
