


/***
 * Surface water change detection (copied from http://aqua-monitor.appspot.com)
 * 
 * Citation: Donchyts, Gennadii, et al. "Earth's surface water change over the past 30 years." Nature Climate Change 6.9 (2016): 810-813.
 *
 * License: GPL
 * 
 */

function renderLandsatMosaic(options, dateIntervalIndex) {
  var percentile = options.percentile;
  var start = options.dateIntervals[dateIntervalIndex][0];
  var stop = options.dateIntervals[dateIntervalIndex][1];
  var sharpen = options.sharpen;
  var smoothen = options.smoothen;
  var filterCount = options.filterCount;
  
  var bands = ['swir1', 'nir', 'green'];
  var l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterDate(start, stop).select(['B6', 'B5', 'B3'], bands);
  var l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').filterDate(start, stop).select(['B5', 'B4', 'B2'], bands);
  var l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA').filterDate(start, stop).select(['B5', 'B4', 'B2'], bands);
  var l4 = ee.ImageCollection('LANDSAT/LT04/C01/T1_TOA').filterDate(start, stop).select(['B5', 'B4', 'B2'], bands);

  var images = l8.merge(l7).merge(l5).merge(l4)
  
  //images = ee.ImageCollection(images.limit(100))
  
  //var images = ee.ImageCollection(l7.merge(l5).merge(l4))

  if(smoothen) {
    images = images.map(function(i) { return i.resample('bicubic'); })
  }

  var image = images
    //.filterMetadata('SUN_AZIMUTH', 'greater_than', 5) // almost empty
    .reduce(ee.Reducer.percentile([percentile]))
    .rename(bands)

  if(filterCount > 0) {
    image = image.mask(images.select(0).count().gt(filterCount));
  }
  
  // Map.addLayer(images.select(0).count(), {min:filterCount, max:200, palette:['ff0000', '00ff00']}, 'count', false)

  if(sharpen) {
    // LoG
    image = image.subtract(image.convolve(ee.Kernel.gaussian(30, 20, 'meters')).convolve(ee.Kernel.laplacian8(0.4)))
  }

  return image.visualize({min: 0.05, max: [0.5, 0.5, 0.6], gamma: 1.4})
}

// A helper to apply an expression and linearly rescale the output.
var rescale = function (img, thresholds) {
  return img.subtract(thresholds[0]).divide(ee.Number(thresholds[1]).subtract(thresholds[0]))
  .copyProperties(img)
  .copyProperties(img, ['system:time_start']);
};

function pad(n, width, z) {
  z = z || '0';
  n = n + '';
  return n.length >= width ? n : new Array(width - n.length + 1).join(z) + n;
}

var generateGrid = function(xmin, ymin, xmax, ymax, dx, dy) {
  var xx = ee.List.sequence(xmin, ee.Number(xmax).subtract(dx), dx);
  var yy = ee.List.sequence(ymin, ee.Number(ymax).subtract(dx), dy);
  var polys = xx.map(function(x) {F
    return yy.map(function(y) {
      var x1 = ee.Number(x)
      var x2 = ee.Number(x).add(dx)
      var y1 = ee.Number(y)
      var y2 = ee.Number(y).add(dy)
      
      var coords = ee.List([x1, y1, x2, y2]);

      return ee.Feature(ee.Algorithms.GeometryConstructors.Rectangle(coords));
    })
  }).flatten()

  return ee.FeatureCollection(polys);
}

function getIntersection(left, right) {
  var spatialFilter = ee.Filter.intersects({leftField: '.geo', rightField: '.geo', maxError: 1000});
  var saveAllJoin = ee.Join.saveAll({matchesKey: 'match'});
  var intersectJoined = saveAllJoin.apply(left, right, spatialFilter);

  return intersectJoined.map(function(f) { 
    var match = ee.List(f.get('match'));
    return f.set('count', match.length())
  }).filter(ee.Filter.gt('count', 0))
}

function getEdge(i) {
  var canny = ee.Algorithms.CannyEdgeDetector(i, 0.99, 0);
  canny = canny.mask(canny)
  return canny;
}

function createTimeBand(img) {
  var date = ee.Date(img.get('system:time_start'));
  var year = date.get('year').subtract(1970);
  return ee.Image(year).byte().addBands(img)
}

// TODO: split it into smaller functions
function renderWaterTrend(options) {
  var dateIntervals = options.dateIntervals;
  var percentile = options.percentile;
  
  var slopeThreshold = options.slopeThreshold;
  var slopeThresholdRatio = options.slopeThresholdRatio;
  var refine = options.refine;
  var slopeThresholdRefined = options.slopeThresholdRefined;
  
  var ndviFilter = options.ndviFilter;
  var filterCount = options.filterCount;

  var showEdges = options.showEdges;
  var smoothen = options.smoothen;
  var includeBackgroundSlope = options.includeBackgroundSlope;
  var backgroundSlopeOpacity = options.backgroundSlopeOpacity;
  var refineFactor = options.refineFactor;
  
  var ndwiMaxLand = options.ndwiMaxLand;
  var ndwiMinWater = options.ndwiMinWater;
  var alwaysWaterThreshold = options.alwaysWaterThreshold;

  var bands = ['green', 'swir1'];
  var bands8 = ['B3', 'B6'];
  var bands7 = ['B2', 'B5'];

  if(ndviFilter > -1) {
    bands = ['green', 'swir1', 'nir', 'red'];
    bands8 = ['B3', 'B6', 'B5', 'B4'];
    bands7 = ['B2', 'B5', 'B4', 'B3'];
  }
  
  var images = new ee.ImageCollection([])

  var images_l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').select(bands8, bands);
  images = images.merge(images_l8);

  var images_l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').select(bands7, bands);
  images = images.merge(images_l7);

  var images_l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA').select(bands7, bands);
  images = images.merge(images_l5);

  var images_l4 = ee.ImageCollection('LANDSAT/LT04/C01/T1_TOA').select(bands7, bands);
  images = images.merge(images_l4);
  
  // images = ee.ImageCollection(images.limit(100))

  var list = ee.List(dateIntervals);
  
  // add percentile images for debugging
  if(options.debugMapLayers) {
    list.getInfo().map(function (i) {
      var start = ee.Date(i[0].value); 
      var stop = ee.Date(i[1].value);
  
      var filtered = images.filterDate(start, stop)
  
      var percentiles = ee.List.sequence(0, 100, 1)
  
      var result = filtered
          .reduce(ee.Reducer.percentile(percentiles))
          .set('system:time_start', start)
          
      Map.addLayer(result, {}, 'all percentiles ' + start.format('YYYY-MM-dd').getInfo(), false)
    });
  }
  
  
  // compute a single annual percentile
  var annualPercentile = ee.ImageCollection(list.map(function (i) {
    var l = ee.List(i);
    var start = l.get(0); 
    var stop = l.get(1);

    var filtered = images.filterDate(start, stop)

    if(smoothen) {
      filtered = filtered.map(function(i) { return i.resample('bicubic'); })
    }

    var image = filtered
        .reduce(ee.Reducer.percentile([percentile])).rename(bands)

    var result = image
        .normalizedDifference(['green', 'swir1']).rename('water')
        .set('system:time_start', start);

    if(ndviFilter > -1) {
      var ndvi = image.normalizedDifference(['nir', 'red']).rename('ndvi');
      result = result.addBands(ndvi)
    }

    if(filterCount > 0) {
        result = result.addBands(filtered.select(0).count().rename('count'));
    }

    return result
  }));
  
  var mndwi = annualPercentile.select('water')

  if(ndviFilter > -1) {
    var ndvi = annualPercentile.select('ndvi')
  }

  var fit = mndwi
    .map(function(img) { 
      return rescale(img, [-0.6, 0.6]);
    })
    .map(createTimeBand)
    .reduce(ee.Reducer.linearFit().unweighted());

  var scale = fit.select('scale')
  
  var scaleMask = scale.mask();

  if(options.debugMapLayers) {
    Map.addLayer(scaleMask, {}, 'scale original mask', false)

    Map.addLayer(scale, {
        min: -slopeThreshold * slopeThresholdRatio,
        max: slopeThreshold * slopeThresholdRatio,
        palette: ['00ff00', '000000', '00d8ff'],
      }, 'scale', false)
  }

  var mndwiMin = mndwi.min();
  var mndwiMax = mndwi.max();
  
  if(ndviFilter > -1) {
    var ndviMin = ndvi.min();
  }
  
  if(options.debugMapLayers) {
    Map.addLayer(mndwiMin, {}, 'mndwi min (raw)', false)
    Map.addLayer(mndwiMax, {}, 'mndwi max (raw)', false)
    if(ndviFilter > -1) {
      Map.addLayer(ndvi.min(), {}, 'ndvi min (raw)', false)
      Map.addLayer(ndvi.max(), {}, 'ndvi max (raw)', false)
    }
  }

  if(options.useSwbdMask) {
    var swbd = ee.Image('MODIS/MOD44W/MOD44W_005_2000_02_24').select('water_mask')
    var swbdMask = swbd.unmask().not()
      .focal_max(10000, 'square', 'meters').reproject('EPSG:4326', null, 1000)
  }  

  // computes a mask representing a surface water change using a given slope (linear fit scale) threshold
  function computeSlopeMask(threshold) {
    var minWaterMask = mndwiMax.gt(ndwiMinWater) // maximum looks like water
    var maxLandMask = mndwiMin.lt(ndwiMaxLand) // minimum looks like land
    
    var mask = scale.unmask().abs().gt(threshold)
      .multiply(minWaterMask) 
      .multiply(maxLandMask) 

    if(ndviFilter > -1) {
      var ndviMask = ndviMin.lt(ndviFilter)
      mask = mask.multiply(ndviMask) // deforestation?
    }

    if(filterCount > 0) {
      var countMask = annualPercentile.select('count').min().gt(filterCount)
      mask = mask.multiply(countMask);
    }
    
    if(options.useSwbdMask) {
      mask = mask.multiply(swbdMask)
    }
    
    // add eroded original scale mask (small scale-friendly erosion, avoid kernel too large)
    // var erodedScaleMask = scaleMask
    //   .focal_min(10000, 'square', 'meters').reproject('EPSG:4326', null, 1000)

    // mask = mask.multiply(erodedScaleMask)

    if(options.debugMapLayers) {
      Map.addLayer(minWaterMask.not().mask(minWaterMask.not()), {}, 'min water mask', false)
      Map.addLayer(maxLandMask.not().mask(maxLandMask.not()), {}, 'max land mask', false)
      
      if(ndviFilter > -1) {
        Map.addLayer(ndviMask.not().mask(ndviMask.not()), {}, 'ndvi mask', false)
      }
      if(filterCount > 0) {
        Map.addLayer(countMask.not().mask(countMask.not()), {}, 'count mask', false)
      }
      
      if(options.useSwbdMask) {
        Map.addLayer(swbdMask.not().mask(swbdMask.not()), {}, 'swbd mask', false)
      }
      
      Map.addLayer(erodedScaleMask.not().mask(erodedScaleMask.not()), {}, 'scale original mask (eroded)', false)
    }

    return mask;
  }

  print('slope threshold: ', slopeThresholdRatio * slopeThreshold)
  var mask = computeSlopeMask(slopeThresholdRatio * slopeThreshold);

  if(refine) {
    // this should be easier, maybe use SRTM projection
    mask = mask.reproject('EPSG:4326', null, 30)
    var prj = mask.projection();

    // more a buffer around larger change
    var maskBuffer = mask
      .reduceResolution(ee.Reducer.max(), true)
      .focal_max(refineFactor)
      .reproject(prj.scale(refineFactor, refineFactor))
      .focal_mode(ee.Number(refineFactor).multiply(30), 'circle', 'meters')

    print('slope threshold (refined): ', slopeThresholdRefined * slopeThresholdRatio)
    var maskRefined = computeSlopeMask(slopeThresholdRefined * slopeThresholdRatio).mask(maskBuffer)

    if(options.debugMapLayers) {
      Map.addLayer(mask.mask(mask), {}, 'mask (raw)', false)
      Map.addLayer(maskBuffer.mask(maskBuffer), {}, 'mask buffer (raw)', false)
      Map.addLayer(maskRefined.mask(maskRefined), {}, 'mask refined (raw)', false)
    }
    
    // smoothen scale and mask
    if(smoothen) {
      scale = scale
        .focal_median(25, 'circle', 'meters', 3);

      mask = mask
        .focal_mode(35, 'circle', 'meters', 3)
    }
  }

  if(options.debugMapLayers) {
    Map.addLayer(scale, {}, 'scale (raw)', false)
  }

  var results = [];

  // background
  var bg = ee.Image(1).toInt().visualize({palette: '000000', opacity: 0.4});
  
  if(includeBackgroundSlope) {
    bg = scale.visualize({
        min: -slopeThreshold * slopeThresholdRatio,
        max: slopeThreshold * slopeThresholdRatio,
        palette: ['00ff00', '000000', '00d8ff'], opacity: backgroundSlopeOpacity
      });

    // exclude when both are water
    if(alwaysWaterThreshold) {
      bg = bg.mask(ee.Image(backgroundSlopeOpacity).toFloat().multiply(mndwiMin.gt(alwaysWaterThreshold).focal_mode(1).not())) 
    }
  } 

  if(filterCount > 0) {
    bg = bg.multiply(annualPercentile.select('count').min().gt(filterCount));
  }

  if(options.useSwbdMask) {
    bg = bg.multiply(swbdMask.gt(0))
  }

  results.push(bg);


  // surface water change
  if(refine) {
    if(options.debug) {
      var maskBufferVis = maskBuffer.mask(maskBuffer).visualize({palette:['ffffff', '000000'], opacity:0.5})
      results.push(maskBufferVis);
    }
    
    var edgeWater = getEdge(mask.mask(scale.gt(0))).visualize({palette: '00d8ff'})
    var edgeLand = getEdge(mask.mask(scale.lt(0))).visualize({palette: '00ff00'})

    scale = scale.mask(maskRefined)

    var scaleRefined = scale.visualize({
      min: -slopeThreshold * slopeThresholdRatio,
      max: slopeThreshold * slopeThresholdRatio,
      palette: ['00ff00', '000000', '00d8ff'],
      opacity: showEdges ? 0.3 : 1.0
    })

    results.push(scaleRefined)
    
    if(showEdges) {
      results.push(edgeWater, edgeLand)
    }
    
  } else {
    scale = scale.mask(mask)

    var change = scale.visualize({
      min: -slopeThreshold * slopeThresholdRatio,
      max: slopeThreshold * slopeThresholdRatio,
      palette: ['00ff00', '000000', '00d8ff'],
    })

    results.push(change);
  }

  return {changeVis: ee.ImageCollection.fromImages(results).mosaic(), change: scale.toFloat()};
}

function computeAggregatedSurfaceWaterChangeArea(scale, options) {
    // add aggregated version of change
    var changeAggregatedWater = scale.gt(0).multiply(ee.Image.pixelArea().divide(1000000))
      .reproject('EPSG:4326', null, 30)
      .reduceResolution(ee.Reducer.sum(), false, 100)
      .reproject('EPSG:4326', null, 300)

var statswater = changeAggregatedWater.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: geometry2, 
  scale: 30,
  maxPixels: 1e13

});

//print(statswater, 'water')

    var changeAggregatedLand = scale.lt(0).multiply(ee.Image.pixelArea().divide(1000000))
      .reproject('EPSG:4326', null, 30)
      .reduceResolution(ee.Reducer.sum(), false, 100)
      .reproject('EPSG:4326', null, 300)

var statsland = changeAggregatedLand.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: geometry2,
  scale: 30,
  maxPixels: 1e13

});

//print(statsland, 'land')

    var maxArea = 300
    var changeAggregatedWaterVis = changeAggregatedWater
      .visualize({
        min: 0,
        max: maxArea,
        palette: ['000000', '00d8ff'],
      })

    var changeAggregatedLandVis = changeAggregatedLand
      .visualize({
        min: 0,
        max: maxArea,
        palette: ['000000', '00ff00'],
      })

    Map.addLayer(changeAggregatedWater, {}, 'scale aggregated water (raw)', false)
    Map.addLayer(changeAggregatedLand, {}, 'scale aggregated land (raw)', false)
    
    Map.addLayer(changeAggregatedWaterVis.mask(changeAggregatedWater.divide(maxArea)), {}, 'change aggregated (water => land)', false)
    Map.addLayer(changeAggregatedLandVis.mask(changeAggregatedLand.divide(maxArea)), {}, 'change aggregated (land => water)', false)
    
    return {water: changeAggregatedWater.toFloat(), land: changeAggregatedLand.toFloat()};
}

function getWaterTrendChangeRatio(start, stop) {
  return (15 / (stop - start));  // empiricaly found ratio
}

// ======================================================= PARAMETERS AND SCRIPT

// start / stop times and averaging periods (in months)
var time0 = [ee.Date.fromYMD(1990, 1, 1), 24]; // [ee.Date.fromYMD(2000, 1, 1), 24];
var time1 = [ee.Date.fromYMD(2020, 1, 1), 24]; //[ee.Date.fromYMD(2015, 1, 1), 24];

// larger periods, more robust, slower, may contain less changes (used to compute results reported in the paper)
// var time0 = [ee.Date.fromYMD(1984, 1, 1), 240];
// var time1 = [ee.Date.fromYMD(2013, 1, 1), 48];

var options = {
  // intervals used for averaging and linear regression (the web site may use multiple intervals here)  
  dateIntervals: [
    [time0[0], time0[0].advance(time0[1], 'month')],
    [time1[0], time1[0].advance(time1[1], 'month')]
  ],

  percentile: 15,

  slopeThreshold: 0.025,
  slopeThresholdRatio: getWaterTrendChangeRatio(1990, 2020),

  slopeThresholdRefined: 0.015,

  // refine: true, // more expensive
  refine: false,
  refineFactor: 5,

  ndviFilter: 0.08, // the highest NDVI value for water
  //ndviFilter: -1,

  ndwiMinWater: -0.05, // minimum value of NDWI to assume as water
  ndwiMaxLand: 0.5, // maximum value of NDWI to assume as land
  
  alwaysWaterThreshold: 0.7,

  filterCount: 10,
  
  //useSwbdMask: true,
  useSwbdMask: false,
  
  showEdges: true,
  // showEdges: false,

  //includeBackgroundSlope: false,
  includeBackgroundSlope: true,

  backgroundSlopeOpacity: 0.2,

  //smoothen: false,
  smoothen: true,

  //debug: false, // shows a buffer used to refine changes
  debug: false,
  //debugMapLayers: true,
  debugMapLayers: false,

  //sharpen: true,
  sharpen: false,
};

print(options.dateIntervals[0][0].format('YYYY-MM-dd').getInfo() + ' - ' + options.dateIntervals[0][1].format('YYYY-MM-dd').getInfo());
print(options.dateIntervals[1][0].format('YYYY-MM-dd').getInfo() + ' - ' + options.dateIntervals[1][1].format('YYYY-MM-dd').getInfo());

// background
Map.addLayer(ee.Image(1).toInt(), {palette:['000000']}, 'bg (black)', false);
Map.addLayer(ee.Image(1).toInt(), {palette:['ffffff']}, 'bg (white)', false);

// average images
var timeCount = options.dateIntervals.length;
Map.addLayer(renderLandsatMosaic(options, 0), {}, options.dateIntervals[0][0].format('YYYY-MM-dd').getInfo(), true);
Map.addLayer(renderLandsatMosaic(options, timeCount - 1), {}, options.dateIntervals[timeCount - 1][0].format('YYYY-MM-dd').getInfo(), true);

// country boundaries
Map.addLayer(countries.map(function(f) { return f.buffer(15000) }), {}, 'countries', false);

// GLCF water
var water = glcf.map(function(i){return i.eq(2)}).mosaic();
Map.addLayer(water.mask(water), {palette:['2020aa'], opacity: 0.5}, 'GLCF water', false);

// surface water change trend
var trend1 = renderWaterTrend(options);
Map.addLayer(trend1.changeVis, {}, 'water change', true);

Export.image.toDrive({
  // image: trend1.change,
  image: trend1.changeVis, 
  // description: , 
  folder: 'Earth Engine folder', 
  fileNamePrefix: '1987_2015_waterchange', 
  // dimensions: , 
  region: geometry2, 
  scale: 30, 
  // crs: , 
  // crsTransform:, 
  // maxPixels:, 
  // shardSize:, 
  // fileDimensions:
});

/*
options.refine = false;
options.debugMapLayers = false;
var trend1 = renderWaterTrend(options)[0];
Map.addLayer(trend1.changeVis, {}, '1987 - 2015 (water change, no refine)', true)
*/

// temporary compute aggregated version here
var trend1Aggregated = computeAggregatedSurfaceWaterChangeArea(trend1.change, options);

// HAND
var handThreshold = 50;
var handBuffer = 150;
Map.addLayer(ee.Image(1).mask(hand.mosaic().gt(handThreshold)
  .focal_max(handBuffer, 'circle', 'meters')
  .focal_min(handBuffer, 'circle', 'meters')
  ), 
  {palette:['000000']}, 'HAND > ' + handThreshold + ' (+' + handBuffer + 'm closing)', false);

Map.addLayer(ee.Image(1).mask(hand.mosaic().gt(50)), {palette:['000000']}, 'HAND > ' + handThreshold + 'm', false);

// set map options
Map.setOptions('SATELLITE')

// center to a specific location
print(Map.getCenter())



/*
Copyright (c) 2018 Gennadii Donchyts. All rights reserved.

This work is licensed under the terms of the MIT license.  
For a copy, see <https://opensource.org/licenses/MIT>.
*/

var g = require('users/gena/packages:grid')

// generateGridForGeometry(g, dx, dy, opt_marginx, opt_marginy, opt_proj)

// WEB
var proj = 'EPSG:3857'
var dx = 400000
var dy = 400000
var marginx = 0, marginy = 0
var grid1 = g.generateGridForGeometry(geometry2, dx, dy, marginx, marginy, proj)
//Map.addLayer(grid1, { color: 'brown' }, 'Web')

// GEO
var proj = 'EPSG:4326'
var dx = 0.02
var dy = 0.02
var grid2 = g.generateGridForGeometry(geometry2, dx, dy, marginx, marginy, proj)
Map.addLayer(grid2, { color: 'black' }, 'GEO')

function erode(img, distance) {
  var d = (img.not().unmask(1)
       .fastDistanceTransform(256).sqrt()
       .multiply(ee.Image.pixelArea().sqrt()))
  return img.updateMask(d.gt(distance))
}

function dilate(img, distance) {
  var d = (img.fastDistanceTransform(256).sqrt()
       .multiply(ee.Image.pixelArea().sqrt()))
  return d.lt(distance)
}

var land = ee.Image("MERIT/DEM/v1_0_3").select('dem').unmask(0).gt(0)
var landMask = erode(dilate(land, 3000), 10000).mask().eq(1)
Map.addLayer(landMask.updateMask(landMask), {palette: "beige"}, "land 10km")
 
var oceanMask = erode(land.not(), 10000).mask().eq(1)
Map.addLayer(oceanMask.updateMask(oceanMask), {palette: "aqua"}, "ocean 10km")

