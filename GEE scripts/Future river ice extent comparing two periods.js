var cmip5 = ee.ImageCollection("NASA/NEX-GDDP"),
    grwl = ee.FeatureCollection("users/eeProject/grwl");

var model = 'CESM1-BGC'; //['CESM1-BGC', 'GFDL-ESM2M', 'MIROC-ESM']
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

var appendPeriodBand = function(i) {
  var doy = ee.Number.parse(ee.Date(i.get('system:time_start')).format('D'));
  var period = ee.Algorithms.If(doy.gt(220).or(doy.lt(40)), 'freeze-up', 'breakup');
  var periodValue = ee.Algorithms.If(doy.gt(220).or(doy.lt(40)), 0, 1);
  return(ee.Image(i.setMulti({'period': period, 'periodValue': periodValue})));
};

var calculateRiverIce = function(i) {
  var iceFraction = ee.Image(1).subtract(ee.Image(1).divide(i.multiply(ee.Number(-0.32233).add(ee.Number(i.get('periodValue')).multiply(-0.04515))).add(-0.82245).exp().add(1)));
  return(i.addBands(iceFraction.rename(['iceFraction'])));
};

var calculateIceFraction = function(i) {
  var temp = calculate30DayPriorMeanTemp(i);
  temp = appendPeriodBand(temp);
  var iceFraction = calculateRiverIce(temp);
  return(iceFraction);
};

var calculateIceFractionPeriod = function(collection, startDate, endDate) {
  var riverIceStats = collection
  .filterDate(startDate, endDate)
  .map(calculateIceFraction)
  .mean();

  return(riverIceStats);
};

var iceOneYear = calculateIceFractionPeriod(dailyTas, '2019-01-01', '2020-01-01').aside(print);

for (var i = 1; i <= 12; i++) {
  var period1 = calculateIceFractionPeriod(dailyTas.filterMetadata('month', 'equals', i), '2008-12-31', '2028-12-31');
  var period2 = calculateIceFractionPeriod(dailyTas.filterMetadata('month', 'equals', i), '2079-12-31', '2099-12-31');

  var diff = period2.subtract(period1).set('month', i);

  Map.addLayer(diff.select(['iceFraction']), {min: -0.5, max: 0.5, palette: ['red', 'grey', 'blue']}, 'month' + i, false);

  Export.image.toDrive({
    image: diff,
    description: 'diff_month_' + model + '_' + scenario + '_' + i,
    folder: 'monthly_map_river_ice_future_change_global',
    fileNamePrefix: 'diff_month_' + model + '_' + scenario + '_' + i,
    dimensions: '1440x720',
    crs: 'EPSG:4326',
    crsTransform: [0.25, 0, -180, 0, 0.25, -90],
    fileFormat: 'GeoTIFF'});
}

// calculate grwl stats

var partitionGRWL = function() {
  var grwlFil = grwl
  .filterMetadata('width_m', 'greater_than', 90)
  .filterMetadata('nchannels', 'equals', 1)
  .filterMetadata('lakeFlag', 'equals', 0)
  .filterMetadata('lat', 'greater_than', 23.5);

  var grwlFil1 = grwlFil
  .filterMetadata('lon', 'greater_than', -180)
  .filterMetadata('lon', 'not_greater_than', -135);

  var grwlFil2 = grwlFil
  .filterMetadata('lon', 'greater_than', -135)
  .filterMetadata('lon', 'not_greater_than', -90);

  var grwlFil3 = grwlFil
  .filterMetadata('lon', 'greater_than', -90)
  .filterMetadata('lon', 'not_greater_than', -45);

  var grwlFil4 = grwlFil
  .filterMetadata('lon', 'greater_than', -45)
  .filterMetadata('lon', 'not_greater_than', 0);

  var grwlFil5 = grwlFil
  .filterMetadata('lon', 'greater_than', 0)
  .filterMetadata('lon', 'not_greater_than', 45);

  var grwlFil6 = grwlFil
  .filterMetadata('lon', 'greater_than', 45)
  .filterMetadata('lon', 'not_greater_than', 90);

  var grwlFil7 = grwlFil
  .filterMetadata('lon', 'greater_than', 90)
  .filterMetadata('lon', 'not_greater_than', 135);

  var grwlFil8 = grwlFil
  .filterMetadata('lon', 'greater_than', 135)
  .filterMetadata('lon', 'not_greater_than', 180);

  var grwlFilCollection = ee.FeatureCollection([
    ee.Feature(grwlFil1.union().first()).set('npoints', grwlFil1.size()),
    ee.Feature(grwlFil2.union().first()).set('npoints', grwlFil2.size()),
    ee.Feature(grwlFil3.union().first()).set('npoints', grwlFil3.size()),
    ee.Feature(grwlFil4.union().first()).set('npoints', grwlFil4.size()),
    ee.Feature(grwlFil5.union().first()).set('npoints', grwlFil5.size()),
    ee.Feature(grwlFil6.union().first()).set('npoints', grwlFil6.size()),
    ee.Feature(grwlFil7.union().first()).set('npoints', grwlFil7.size()),
    ee.Feature(grwlFil8.union().first()).set('npoints', grwlFil8.size())
    ]);

  return(grwlFilCollection);
};

var grwlFilCollection = partitionGRWL();

// var month = 1;

var tempGen = function(i) {
  var temp = function(l) {
    var period1 = calculateIceFractionPeriod(dailyTas.filterMetadata('month', 'equals', l), '2008-12-31', '2028-12-31').select(['iceFraction']);
    var period2 = calculateIceFractionPeriod(dailyTas.filterMetadata('month', 'equals', l), '2079-12-31', '2099-12-31').select(['iceFraction']);

    var diff = period2.subtract(period1).set('month', l);

    var geometry = ee.Feature(grwlFilCollection.toList(1, i).get(0)).geometry();

    var predictedRasters = period1.rename(['period1'])
    .addBands(period2.rename(['period2']))
    .addBands(diff.rename(['diff']));

    var grwlStats = predictedRasters.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geometry,
        crs: 'EPSG:4326',
        crsTransform: [0.25, 0, 0, 0, 0.25, 0],
        tileScale: 1,
        maxPixels: 100000000000
      });

    return(ee.Feature(null, grwlStats).set('month', l));
  };
  return(temp);
};

for (var i = 0; i <= 7; i++) {
  var predictedMonthlyStats = ee.List.sequence(1, 12, 1).map(tempGen(i));

  Export.table.toDrive({
    collection: ee.FeatureCollection(predictedMonthlyStats),
    description: 'monthlyStats' + '_' + model + '_' + scenario + '_' + i,
    folder: 'predicted_monthly_river_ice_change_number',
    fileNamePrefix: 'monthlyStats' + '_' + model + '_' + scenario + '_' + i,
    fileFormat: 'csv'});

}

print(ee.FeatureCollection(predictedMonthlyStats).first());


//  calculate river ice duration map
var calculateIceStatusPeriod = function(collection, startDate, endDate) {
  var riverIceStats = collection
  .filterDate(startDate, endDate)
  .map(function(i) {
    var temp = calculateIceFraction(i);
    return(temp.gte(0.5));
  })
  .mean();

  return(riverIceStats);
};


var duration1 = calculateIceStatusPeriod(dailyTas, '2008-12-31', '2028-12-31')
.select(['iceFraction']).multiply(365).rename(['duration']);
var duration2 = calculateIceStatusPeriod(dailyTas, '2079-12-31', '2099-12-31')
.select(['iceFraction']).multiply(365).rename(['duration']);

var duration_diff = duration2.subtract(duration1);

print(duration1);
Map.addLayer(duration1, {min: 0, max: 365, palette: ['orange', 'white', 'cyan']});

Export.image.toDrive({
    image: duration1,
    description: 'diff_duration_' + model + '_' + scenario + '_' + '08-28',
    folder: 'duration_change_maps',
    fileNamePrefix: 'diff_duration_' + model + '_' + scenario + '_' + '08-28',
    dimensions: '1440x720',
    crs: 'EPSG:4326',
    crsTransform: [0.25, 0, -180, 0, 0.25, -90],
    fileFormat: 'GeoTIFF'});

Export.image.toDrive({
    image: duration2,
    description: 'diff_duration_' + model + '_' + scenario + '_' + '79-99',
    folder: 'duration_change_maps',
    fileNamePrefix: 'diff_duration_' + model + '_' + scenario + '_' + '79-99',
    dimensions: '1440x720',
    crs: 'EPSG:4326',
    crsTransform: [0.25, 0, -180, 0, 0.25, -90],
    fileFormat: 'GeoTIFF'});

Export.image.toDrive({
    image: duration_diff,
    description: 'diff_duration_' + model + '_' + scenario + '_' + 'diff',
    folder: 'duration_change_maps',
    fileNamePrefix: 'diff_duration_' + model + '_' + scenario + '_' + 'diff',
    dimensions: '1440x720',
    crs: 'EPSG:4326',
    crsTransform: [0.25, 0, -180, 0, 0.25, -90],
    fileFormat: 'GeoTIFF'});

// calculate river ice duration grwl stats

var tempGen = function(i) {
  var temp = function() {
    var period1 = calculateIceStatusPeriod(dailyTas, '2008-12-31', '2028-12-31')
    .select(['iceFraction']).multiply(365).rename(['duration1']);
    var period2 = calculateIceStatusPeriod(dailyTas, '2079-12-31', '2099-12-31')
    .select(['iceFraction']).multiply(365).rename(['duration2']);

    var diff = period2.subtract(period1);

    diff = diff.mask(period1.gt(0));

    var geometry = ee.Feature(grwlFilCollection.toList(1, i).get(0)).geometry();

    var predictedRasters = period1.rename(['duration1'])
    .addBands(period2.rename(['duration2']))
    .addBands(diff.rename(['diff_duration_all']))
    .addBands(diff.mask(period1.gt(0)).rename(['diff_duration_d1nonZero']))
    .addBands(diff.mask(period1.gt(0).and(period2.gt(0))).rename(['diff_duration_allnonZero']));

    var grwlStats = predictedRasters.reduceRegion({
        reducer: ee.Reducer.mean().combine(ee.Reducer.count(), null, true),
        geometry: geometry,
        crs: 'EPSG:4326',
        crsTransform: [0.25, 0, 0, 0, 0.25, 0],
        tileScale: 1,
        maxPixels: 100000000000
      });

    return(ee.Feature(null, grwlStats));
  };
  return(temp);
};

for (var i = 0; i <= 7; i++) {
  var predictedMonthlyStats = tempGen(i)();

  Export.table.toDrive({
    collection: ee.FeatureCollection(predictedMonthlyStats),
    description: 'monthlyDurationMultipleStats' + '_' + model + '_' + scenario + '_' + i,
    folder: 'predicted_monthly_river_ice_change_number',
    fileNamePrefix: 'monthlyDurationMultipleStats' + '_' + model + '_' + scenario + '_' + i,
    fileFormat: 'csv'});

}
