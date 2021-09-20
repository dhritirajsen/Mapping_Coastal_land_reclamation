// SOurce code:- https://developers.google.com/earth-engine/tutorials/community/extract-raster-values-for-points
function bufferPoints(radius, bounds) {
  return function(pt) {
    pt = ee.Feature(pt);
    return bounds ? pt.buffer(radius).bounds() : pt.buffer(radius);
  };
}

function zonalStats(ic, fc, params) {
  // Initialize internal params dictionary.
  var _params = {
    reducer: ee.Reducer.mean(),
    scale: null,
    crs: null,
    bands: null,
    bandsRename: null,
    imgProps: null,
    imgPropsRename: null,
    datetimeName: 'datetime',
    datetimeFormat: 'YYYY-MM-dd HH:MM:ss'
  };

  // Replace initialized params with provided params.
  if (params) {
    for (var param in params) {
      _params[param] = params[param] || _params[param];
    }
  }

  // Set default parameters based on an image representative.
  var imgRep = ic.first();
  var nonSystemImgProps = ee.Feature(null)
    .copyProperties(imgRep).propertyNames();
  if (!_params.bands) _params.bands = imgRep.bandNames();
  if (!_params.bandsRename) _params.bandsRename = _params.bands;
  if (!_params.imgProps) _params.imgProps = nonSystemImgProps;
  if (!_params.imgPropsRename) _params.imgPropsRename = _params.imgProps;

  // Map the reduceRegions function over the image collection.
  var results = ic.map(function(img) {
    // Select bands (optionally rename), set a datetime & timestamp property.
    img = ee.Image(img.select(_params.bands, _params.bandsRename))
      .set(_params.datetimeName, img.date().format(_params.datetimeFormat))
      .set('timestamp', img.get('system:time_start'));

    // Define final image property dictionary to set in output features.
    var propsFrom = ee.List(_params.imgProps)
      .cat(ee.List([_params.datetimeName, 'timestamp']));
    var propsTo = ee.List(_params.imgPropsRename)
      .cat(ee.List([_params.datetimeName, 'timestamp']));
    var imgProps = img.toDictionary(propsFrom).rename(propsFrom, propsTo);

    // Subset points that intersect the given image.
    var fcSub = fc.filterBounds(img.geometry());

    // Reduce the image by regions.
    return img.reduceRegions({
      collection: fcSub,
      reducer: _params.reducer,
      scale: _params.scale,
      crs: _params.crs
    })
    // Add metadata to each feature.
    .map(function(f) {
      return f.set(imgProps);
    });
  }).flatten().filter(ee.Filter.notNull(_params.bandsRename));

  return results;
}
var pts = ee.FeatureCollection([
  ee.Feature(Shanghai, {plot_id: 1}),
  ee.Feature(SNG, {plot_id: 2}),
   ee.Feature(Mombasa, {plot_id: 3}),
   ee.Feature(DaresSalaam, {plot_id: 4}),
   ee.Feature(Luanda, {plot_id: 5}),
   ee.Feature(PointeNoire, {plot_id: 6}),
   ee.Feature(Algiers, {plot_id: 7}),
   ee.Feature(Alexandria, {plot_id: 8}),
   ee.Feature(Tripoli, {plot_id: 9}),
   ee.Feature(Casablanca, {plot_id: 10}),
   ee.Feature(Rabat, {plot_id: 11}),
   ee.Feature(Tanger, {plot_id: 12}),
   ee.Feature(Tunis, {plot_id: 13}),
   ee.Feature(PE, {plot_id: 14}),
   ee.Feature(Abidjan, {plot_id: 15}),
   ee.Feature(Accra, {plot_id: 16}),
   ee.Feature(Conakry, {plot_id: 17}),
   ee.Feature(Nouakchott, {plot_id: 18}),
   ee.Feature(Lagos, {plot_id: 19}),
   ee.Feature(Dakar, {plot_id: 20}),
   ee.Feature(Freetown, {plot_id: 21}),
   ee.Feature(Lome, {plot_id: 22}),
   ee.Feature(Shantou, {plot_id: 23}),
   ee.Feature(Dalian, {plot_id: 24}),
   ee.Feature(Dalian, {plot_id: 25}),
   ee.Feature(Dongying, {plot_id: 26}),
   ee.Feature(Guangzhou, {plot_id: 27}),
   ee.Feature(Haikou, {plot_id: 28}),
   ee.Feature(Jiaxing, {plot_id: 29}),
   ee.Feature(Ningbo, {plot_id: 30}),
   ee.Feature(Putian, {plot_id: 31}),
   ee.Feature(Qingdao, {plot_id: 32}),
   ee.Feature(Qinhuangdao, {plot_id: 33}),
   ee.Feature(Rizhao, {plot_id: 34}),
   ee.Feature(Shenzhen, {plot_id: 35}),
    ee.Feature(Tangshan, {plot_id: 36}),
    ee.Feature(Tianjin, {plot_id: 37}),
    ee.Feature(Weihai, {plot_id: 38}),
     ee.Feature(Wenzhou, {plot_id: 39}),
       ee.Feature(Xiamen, {plot_id: 40}),
       ee.Feature(Yantai, {plot_id: 41}),
       ee.Feature(Yingkou, {plot_id: 42}),
       ee.Feature(Zhuhai, {plot_id: 43}),
       ee.Feature(HK, {plot_id: 44}),
       ee.Feature(Kaohsiung, {plot_id: 45}),
       ee.Feature(Taipei, {plot_id: 46}),
       ee.Feature(Taichung, {plot_id: 47}),
       ee.Feature(Tokoname, {plot_id: 48}),
       ee.Feature(Hiroshima, {plot_id: 49}),
       ee.Feature(Osaka, {plot_id: 50}),
       ee.Feature(Kitakyushu, {plot_id: 51}),
       ee.Feature(Tokyo, {plot_id: 52}),
       ee.Feature(Busan, {plot_id: 53}),
       ee.Feature(Incheon, {plot_id: 54}),
       ee.Feature(Chittagong, {plot_id: 55}),
       ee.Feature(Chennai, {plot_id: 56}),
       ee.Feature(Kochi, {plot_id: 57}),
       ee.Feature(Mumbai, {plot_id: 58}),
       ee.Feature(Thiruvanthapuram, {plot_id: 59}),
       ee.Feature(Visakhapatnam, {plot_id: 60}),
       ee.Feature(Karachi, {plot_id: 61}),
       ee.Feature(Bandar_Lampung, {plot_id: 62}),
       ee.Feature(Jakarta, {plot_id: 63}),
       ee.Feature(Makassar, {plot_id: 64}),
       ee.Feature(Samarang, {plot_id: 65}),
       ee.Feature(Surabaya, {plot_id: 66}),
       ee.Feature(Manila, {plot_id: 67}),
       ee.Feature(Da_Nang, {plot_id: 68}),
       ee.Feature(Haifa, {plot_id: 69}),
       ee.Feature(Kuwait_city, {plot_id: 70}),
       ee.Feature(Beirut, {plot_id: 71}),
       ee.Feature(Muscat, {plot_id: 72}),
       ee.Feature(Ad_Dammam, {plot_id: 73}),
       ee.Feature(Istanbul, {plot_id: 74}),
        ee.Feature(Abu_dhabi, {plot_id: 75}),
         ee.Feature(Sharjah, {plot_id: 76}),
         ee.Feature(Dubai, {plot_id: 77}),
         ee.Feature(St_petersburg, {plot_id: 78}),
         ee.Feature(Copenhagen, {plot_id: 79}),
         ee.Feature(Helsinki, {plot_id: 80}),
         ee.Feature(Athens, {plot_id: 81}),
         ee.Feature(Naples, {plot_id: 82}),
         ee.Feature(Porto, {plot_id: 83}),
         ee.Feature(Barcelona, {plot_id: 84}),
         ee.Feature(Amsterdam, {plot_id: 85}),
         ee.Feature(Rotterdam, {plot_id: 86}),
         ee.Feature(Santo_Domingo, {plot_id: 87}),
         ee.Feature(Port_au_Prince, {plot_id: 88}),
         ee.Feature(San_Juan, {plot_id: 89}),
         ee.Feature(Panama_city, {plot_id: 90}),
         ee.Feature(Florianoplis, {plot_id: 91}),
         ee.Feature(Fortaleza, {plot_id: 92}),
         ee.Feature(Rio_de_Janerio, {plot_id: 93}),
         ee.Feature(Cartagena, {plot_id: 94}),
         ee.Feature(Montevideo, {plot_id: 95}),
         ee.Feature(LA, {plot_id: 96}),
         ee.Feature(Adelaide, {plot_id: 97}),
        ee.Feature(Sydney, {plot_id: 98}),
             ee.Feature(Auckland, {plot_id: 99}),
                 ee.Feature(Buenos_Aires, {plot_id: 100}),  
    ee.Feature(Jiddah, {plot_id: 101}), 
     ee.Feature(Zhangjiang, {plot_id: 102}), 
      ee.Feature(Brisbane, {plot_id: 103}), 
      ee.Feature(Lima, {plot_id: 104}), 
      ee.Feature(Maputo, {plot_id: 105}),  
      ee.Feature(Qinghuangdao, {plot_id: 106}),
       ee.Feature(Zhongshan, {plot_id: 107}),
 
]);
Map.addLayer(pts)
var ptsTopo = pts.map(bufferPoints(90, false));
// Import the MERIT global elevation dataset.
var elev = ee.Image('MERIT/DEM/v1_0_3').select('dem');
Map.addLayer(elev, {min:0, max:10}, 'DEM');
// Calculate slope from the DEM.
var slope = ee.Terrain.slope(elev);

// Concatenate elevation and slope as two bands of an image.
var topo = ee.Image.cat(elev, slope)
  // Computed images do not have a 'system:time_start' property; add one based
  // on when the data were collected.
  .set('system:time_start', ee.Date('2000-01-01').millis());

// Wrap the single image in an ImageCollection for use in the zonalStats function.
var topoCol = ee.ImageCollection([topo]);
// Define parameters for the zonalStats function.
var params = {
  bands: [0, 1],
  bandsRename: ['elevation', 'slope']
};

// Extract zonal statistics per point per image.
var ptsTopoStats = zonalStats(topoCol, ptsTopo, params);
print(ptsTopoStats);

Export.table.toDrive({
  collection: ptsTopoStats,
  description: 'your_summary_table_name_here',
  fileFormat: 'CSV'
});