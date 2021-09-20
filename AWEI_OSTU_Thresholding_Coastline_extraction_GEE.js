//Gaining or losing ground? Tracking Asia's hunger for ‘new’ coastal land in the era of sea level rise. 
// Coastline extraction. Modofied by Dhritiraj Sengupta- Original code credit- GEE forum
//GEE code also used in the s tudy '21st-century global coastal reclamation is driven by trade, real-estate growth, and prestige.'
var table =Shanghai;

var cld = require('users/fitoprincipe/geetools:cloud_masks');

var Options = { 
  bandSelect: ['green', 'swir1', 'swir2', 'nir'],
  bands8: ['B3', 'B6', 'B7', 'B5'],
  bands7: ['B2', 'B5', 'B7', 'B4']
};
var landsatFunctions = {
   applyAWEI: function(image) {
    var awei = image.expression("4*(b('green')-b('swir1'))-(0.25*b('nir')+2.75*b('swir2'))");
    return awei.select([0], ['awei']);
  }
};
function FindLandMask(im) {
  var thresh = ee.Number(0.75);
  var sigma = ee.Number(1);
  var scale = ee.Number(400);
  
  im = ee.Image(im).clip(Shanghai);
  var sf = im.get('system:footprint');
  var edge = ee.Algorithms.CannyEdgeDetector(im,thresh,sigma).focal_max(1);
  var hist_info = ee.Dictionary(im.mask(edge).reduceRegion(
      ee.Reducer.histogram(510), table, scale).get('awei'));
  var hist = hist_info.get('histogram');
  var buckets = hist_info.get('bucketMeans');
  var erf = ee.Array(hist).erf();
  var len = ee.List(hist).length();
  var counter = ee.List.sequence(1,len.subtract(1));
  
  // OTSU to build binary image
  function OTSU(i) {
      
      var len = ee.List(hist).length();
      var total = ee.Array(hist).accum(0).get([-1]);
      var q_L = ee.Array(ee.List(hist).slice(0,ee.Number(i))).divide(total).accum(0).get([-1]);
      var q_H = ee.Array(ee.List(hist).slice(ee.Number(i))).divide(total).accum(0).get([-1]);
      var miu_L = ee.Array(ee.List(hist).slice(0,ee.Number(i))).divide(total)
                    .multiply(ee.Array(ee.List.sequence(1,ee.Number(i))))
                    .accum(0).get([-1]).divide(q_L);
      var miu_H = ee.Array(ee.List(hist).slice(ee.Number(i))).divide(total)
                    .multiply(ee.Array(ee.List.sequence(ee.Number(i),ee.Number(len).subtract(1))))
                    .accum(0).get([-1]).divide(q_H);
      var sigma = q_L.multiply(q_H).multiply((miu_L.subtract(miu_H)).pow(2));
      return sigma;
    
  }
  var out = counter.map(OTSU);
  var maxxer = out.reduce(ee.Reducer.max());
  var threshold_index = ee.Number(out.indexOf(maxxer)).subtract(1);
  var threshold = ee.List(buckets).get(threshold_index.round()); 
  threshold = ee.Number(threshold);
  
  // Find thresholds for low and high ranges
  var hist_l = ee.List(hist).slice(threshold_index.round());
  var buckets_l = ee.List(buckets).slice(threshold_index.round());
  var threshold_l = ee.Array(buckets_l).mask(ee.Array(hist_l).erf().gte(0.99)).reduce(ee.Reducer.median(), [0]).get([0]);
  var hist_w = ee.List(hist).slice(0, threshold_index.round());
  var buckets_w = ee.List(buckets).slice(0, threshold_index.round());
  var threshold_w = ee.Array(buckets_w).mask(ee.Array(hist_w).erf().gte(0.99)).reduce(ee.Reducer.median(), [0]).get([0]);
  
  // Creat land-water image
  var land = im.where(im.gte(threshold),1).where(im.lt(threshold),0);
  print(threshold)
  
  return ee.Image(land)
    .set({
      'thresh_mid':threshold,
      'thresh_low':threshold_w,
      'thresh_high':threshold_l,
      'histogram':hist_info,
      'system:footprint': sf
    });
}
var clipImage=function(image){image=image.clip(Shanghai);return image;}

function generateLandsatCollection(start,stop){
  var L4collection = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')
      .filterDate(start, stop)
      .map(cld.landsatSR()).filterBounds(table).map(clipImage)
      .select(Options.bands7, Options.bandSelect);
  var L5collection = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
      .filterDate(start, stop)
      .map(cld.landsatSR()).filterBounds(table).map(clipImage)
      .select(Options.bands7, Options.bandSelect);
  var L7collection = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
          .filterDate(start, stop)
      .map(cld.landsatSR()).filterBounds(table).map(clipImage)
      .select(Options.bands7, Options.bandSelect);
  var L8collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
          .filterDate(start, stop)
      .map(cld.landsatSR()).filterBounds(table).map(clipImage)
      .select(Options.bands8, Options.bandSelect);
  var collectionFull = ee.ImageCollection(L4collection
      .merge(L5collection)
      .merge(L7collection)
      .merge(L8collection))
      
  return collectionFull;
}

var years = ee.List.sequence(2000,2010);

years.evaluate(function(y) {
  y.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var stop =  ee.Date.fromYMD(y, 12, 31);
  var collection = generateLandsatCollection(start,stop);
  var AWEI = collection.map(landsatFunctions.applyAWEI).median();
   
 
  //---Extract  coastline----//
  AWEI=FindLandMask(AWEI.clip(Shanghai)).int()

  //Remove inland water
  var coastline= AWEI.addBands(AWEI)
   .reduceConnectedComponents(ee.Reducer.median(), 'awei', 100).unmask(99).eq(99).and(AWEI.neq(0));
  
  //Remove small islands  
 coastline= coastline.addBands(coastline)
  .reduceConnectedComponents(ee.Reducer.median(), 'awei_1', 1).unmask(99).eq(99).and(coastline.neq(1));    
    
  coastline=coastline.set({name:ee.String(ee.Number(y).int())});

  Export.image.toDrive({
    image: coastline,
    region:table,
    scale:30,
    maxPixels:1000000000000,
    description: y.toString()
    });
  }); 
});
