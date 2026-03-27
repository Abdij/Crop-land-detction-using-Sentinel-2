/************************************************************
 * FARMLAND DETECTION - SOMALIA AGRICULTURAL AREAS
 * Using Sentinel-2 Only with Somalia Boundary
 ************************************************************/

// ============================================================================
// 1. LOAD SOMALIA BOUNDARY
// ============================================================================

var somalia = ee.FeatureCollection('FAO/GAUL/2015/level0')
  .filter(ee.Filter.eq('ADM0_NAME', 'Somalia'));

var somaliaLSIB = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
  .filter(ee.Filter.eq('country_na', 'Somalia'));

// Use GAUL as primary
var somaliaBoundary = somalia;
var studyArea = somaliaBoundary.geometry();

print('Somalia Boundary:', somaliaBoundary);
Map.centerObject(somaliaBoundary, 6);
Map.addLayer(somaliaBoundary, {color: 'red'}, 'Somalia Boundary', true);

// ============================================================================
// 2. LOAD SENTINEL-2 IMAGERY
// ============================================================================

print('==================================================');
print('LOADING SENTINEL-2 IMAGERY FOR SOMALIA');
print('==================================================');

var sentinel2_2024 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2024-01-01', '2024-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

var sentinel2_2023 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2023-01-01', '2023-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

var sceneCount2024 = sentinel2_2024.size();
var sceneCount2023 = sentinel2_2023.size();

print('Sentinel-2 scenes (2024):', sceneCount2024);
print('Sentinel-2 scenes (2023):', sceneCount2023);

// Use 2024 if available, otherwise 2023
var useCollection = ee.ImageCollection(
  ee.Algorithms.If(sceneCount2024.gt(0), sentinel2_2024, sentinel2_2023)
);

var selectedYear = ee.String(
  ee.Algorithms.If(sceneCount2024.gt(0), '2024', '2023')
);

var selectedSceneCount = ee.Number(
  ee.Algorithms.If(sceneCount2024.gt(0), sceneCount2024, sceneCount2023)
);

// ============================================================================
// 3. CLOUD MASKING AND PREPROCESSING
// ============================================================================

function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000).copyProperties(image, image.propertyNames());
}

var processedS2 = useCollection.map(maskS2clouds);
var composite = processedS2.median().clip(studyArea);

// ============================================================================
// 4. CALCULATE VEGETATION INDICES
// ============================================================================

var ndvi = composite.normalizedDifference(['B8', 'B4']).rename('NDVI');

var evi = composite.expression(
  '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
  {
    NIR: composite.select('B8'),
    RED: composite.select('B4'),
    BLUE: composite.select('B2')
  }
).rename('EVI');

var savi = composite.expression(
  '((NIR - RED) / (NIR + RED + 0.5)) * 1.5',
  {
    NIR: composite.select('B8'),
    RED: composite.select('B4')
  }
).rename('SAVI');

var ndmi = composite.normalizedDifference(['B8', 'B11']).rename('NDMI');

var bsi = composite.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))',
  {
    SWIR1: composite.select('B11'),
    RED: composite.select('B4'),
    NIR: composite.select('B8'),
    BLUE: composite.select('B2')
  }
).rename('BSI');

var msavi = composite.expression(
  '(2 * NIR + 1 - sqrt(pow((2 * NIR + 1), 2) - 8 * (NIR - RED))) / 2',
  {
    NIR: composite.select('B8'),
    RED: composite.select('B4')
  }
).rename('MSAVI');

// ============================================================================
// 5. FARMLAND CLASSIFICATION
// ============================================================================

var vegetated = ndvi.gt(0.2).rename('vegetated');
var notForest = ndvi.lt(0.7).rename('notForest');
var soilAdjusted = msavi.gt(0.1).rename('soilAdjusted');
var hasMoisture = ndmi.gt(-0.1).rename('hasMoisture');
var notBare = bsi.lt(0).rename('notBare');
var agriVegetation = evi.gt(0.1).rename('agriVegetation');

var farmland = vegetated
  .and(notForest)
  .and(soilAdjusted)
  .and(hasMoisture)
  .and(notBare)
  .and(agriVegetation)
  .rename('farmland');

// ============================================================================
// 6. AREA FILTERING AND SMOOTHING
// ============================================================================

var minAreaHa = 0.5;
var minAreaPixels = (minAreaHa * 10000) / (10 * 10); // 50 pixels at 10 m

var pixelCount = farmland.selfMask().connectedPixelCount({
  maxSize: 256,
  eightConnected: true
});

var filteredFarmland = farmland.updateMask(pixelCount.gte(minAreaPixels));

var smoothedFarmland = filteredFarmland.focal_mode({
  radius: 1,
  kernelType: 'square',
  iterations: 1
}).rename('farmland');

// ============================================================================
// 7. CREATE COMPOSITES
// ============================================================================

var naturalColor = composite.select(['B4', 'B3', 'B2']);
var falseColor = composite.select(['B8', 'B4', 'B3']);
var agriColor = composite.select(['B11', 'B8', 'B4']);

// ============================================================================
// 8. VISUALIZATION
// ============================================================================

Map.addLayer(naturalColor, {min: 0, max: 0.3}, 'Natural Color (RGB)', true);
Map.addLayer(falseColor, {min: 0, max: 0.5}, 'False Color (Vegetation=Red)', false);
Map.addLayer(agriColor, {min: 0, max: 0.5}, 'Agriculture Composite', false);

Map.addLayer(
  ndvi,
  {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
  'NDVI',
  false
);

Map.addLayer(
  msavi,
  {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
  'MSAVI',
  false
);

Map.addLayer(
  ndmi,
  {min: -0.3, max: 0.4, palette: ['brown', 'yellow', 'green', 'blue']},
  'NDMI',
  false
);

Map.addLayer(
  smoothedFarmland.selfMask(),
  {palette: ['limegreen'], opacity: 0.7},
  'Farmland Areas',
  true
);

// ============================================================================
// 9. LEGEND
// ============================================================================

var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '12px 16px',
    backgroundColor: 'white',
    border: '2px solid #333333',
    width: '260px'
  }
});

legend.add(ui.Label({
  value: 'SOMALIA FARMLAND DETECTION',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    color: '#2c3e50',
    margin: '0 0 10px 0'
  }
}));

legend.add(ui.Label({
  value: 'FARMLAND CLASSES',
  style: {
    fontWeight: 'bold',
    fontSize: '11px',
    color: '#34495e',
    margin: '12px 0 4px 0'
  }
}));

var classItems = [
  {color: '#32cd32', name: 'Farmland Areas', desc: 'Cultivated crops detected'},
  {color: '#ffff00', name: 'Potential Farmland', desc: 'Low confidence areas'},
  {color: '#ff0000', name: 'Non-Farmland', desc: 'Urban/bare soil/water/forest'}
];

classItems.forEach(function(item) {
  var row = ui.Panel({
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {margin: '4px 0', padding: '2px 4px', backgroundColor: '#f8f9fa'}
  });

  var colorBox = ui.Label({
    value: '■',
    style: {color: item.color, fontSize: '16px', margin: '0 8px 0 0'}
  });

  var textPanel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical')
  });

  textPanel.add(ui.Label({
    value: item.name,
    style: {fontWeight: 'bold', fontSize: '11px', color: '#2c3e50', margin: '0'}
  }));

  textPanel.add(ui.Label({
    value: item.desc,
    style: {fontSize: '9px', color: '#7f8c8d', margin: '0'}
  }));

  row.add(colorBox);
  row.add(textPanel);
  legend.add(row);
});

legend.add(ui.Label({
  value: 'DETECTION CRITERIA',
  style: {
    fontWeight: 'bold',
    fontSize: '11px',
    color: '#34495e',
    margin: '12px 0 4px 0'
  }
}));

[
  '• NDVI > 0.2 (vegetated)',
  '• NDVI < 0.7 (not forest)',
  '• MSAVI > 0.1 (soil-adjusted)',
  '• NDMI > -0.1 (has moisture)',
  '• BSI < 0 (not bare soil)',
  '• EVI > 0.1 (agricultural vegetation)',
  '• Min area: 0.5 hectares'
].forEach(function(c) {
  legend.add(ui.Label({
    value: c,
    style: {fontSize: '9px', color: '#34495e', margin: '2px 0'}
  }));
});

var footer = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {margin: '12px 0 0 0', padding: '8px 0 0 0'}
});

footer.add(ui.Label({value: 'Data: Sentinel-2 SR Harmonized', style: {fontSize: '9px', color: '#7f8c8d'}}));
footer.add(ui.Label({value: 'Resolution: 10m', style: {fontSize: '9px', color: '#7f8c8d'}}));
footer.add(ui.Label({value: 'Boundary: FAO GAUL 2015', style: {fontSize: '9px', color: '#7f8c8d'}}));
footer.add(ui.Label({value: 'For agricultural planning only', style: {fontSize: '9px', color: '#e74c3c', margin: '4px 0 0 0'}}));

legend.add(footer);
Map.add(legend);

// ============================================================================
// 10. CALCULATE FARMLAND STATISTICS
// ============================================================================

var pixelArea = ee.Image.pixelArea();

var farmlandArea = smoothedFarmland.selfMask()
  .multiply(pixelArea)
  .rename('farmland_area_m2')
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: studyArea,
    scale: 10,
    maxPixels: 1e13,
    bestEffort: true
  });

var totalArea = pixelArea.rename('area_m2').reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: studyArea,
  scale: 10,
  maxPixels: 1e13,
  bestEffort: true
});

var farmlandHa = ee.Number(farmlandArea.get('farmland_area_m2')).divide(10000);
var totalHa = ee.Number(totalArea.get('area_m2')).divide(10000);
var percentFarmland = farmlandHa.divide(totalHa).multiply(100);

print('==================================================');
print('SOMALIA FARMLAND STATISTICS');
print('==================================================');
print('Total Somalia Area (ha):', totalHa);
print('Farmland Area (ha):', farmlandHa);
print('Percentage Farmland:', percentFarmland);

// ============================================================================
// 11. STATS PANEL
// ============================================================================

var statsPanel = ui.Panel({
  style: {
    position: 'top-right',
    padding: '12px 16px',
    backgroundColor: 'white',
    border: '2px solid #333333',
    width: '240px'
  }
});

statsPanel.add(ui.Label({
  value: 'FARMLAND STATISTICS',
  style: {fontWeight: 'bold', fontSize: '13px', color: '#2c3e50', margin: '0 0 8px 0'}
}));

var totalLabel = ui.Label('Total Somalia Area: loading...', {fontSize: '10px', margin: '4px 0'});
var farmlandLabel = ui.Label('Farmland Detected: loading...', {fontSize: '10px', margin: '4px 0', fontWeight: 'bold', color: '#32cd32'});
var percentLabel = ui.Label('Percentage: loading...', {fontSize: '10px', margin: '4px 0', fontWeight: 'bold'});
var scenesLabel = ui.Label('Sentinel-2 Scenes: loading...', {fontSize: '9px', color: '#7f8c8d', margin: '8px 0 2px 0'});
var yearLabel = ui.Label('Year used: loading...', {fontSize: '9px', color: '#7f8c8d', margin: '2px 0'});
var resLabel = ui.Label('Resolution: 10m', {fontSize: '9px', color: '#7f8c8d', margin: '2px 0'});

statsPanel.add(totalLabel);
statsPanel.add(farmlandLabel);
statsPanel.add(percentLabel);
statsPanel.add(scenesLabel);
statsPanel.add(yearLabel);
statsPanel.add(resLabel);

Map.add(statsPanel);

// Evaluate server-side values before putting them in UI strings
totalHa.format('%.0f').evaluate(function(v) {
  totalLabel.setValue('Total Somalia Area: ' + v + ' ha');
});

farmlandHa.format('%.0f').evaluate(function(v) {
  farmlandLabel.setValue('Farmland Detected: ' + v + ' ha');
});

percentFarmland.format('%.1f').evaluate(function(v) {
  percentLabel.setValue('Percentage: ' + v + '%');
});

selectedSceneCount.evaluate(function(v) {
  scenesLabel.setValue('Sentinel-2 Scenes: ' + v);
});

selectedYear.evaluate(function(v) {
  yearLabel.setValue('Year used: ' + v);
});

// ============================================================================
// 12. EXPORTS
// ============================================================================

Export.image.toDrive({
  image: smoothedFarmland.selfMask(),
  description: 'Somalia_Farmland_Sentinel2',
  folder: 'GEE_Exports',
  fileNamePrefix: 'somalia_farmland',
  region: studyArea,
  scale: 10,
  maxPixels: 1e13
});

var statsFeature = ee.Feature(null, {
  total_area_ha: totalHa,
  farmland_area_ha: farmlandHa,
  farmland_percentage: percentFarmland,
  sentinel2_scenes: selectedSceneCount,
  year_used: selectedYear,
  date_processed: ee.Date(Date.now()).format('YYYY-MM-dd HH:mm:ss')
});

Export.table.toDrive({
  collection: ee.FeatureCollection([statsFeature]),
  description: 'Somalia_Farmland_Statistics',
  folder: 'GEE_Exports',
  fileNamePrefix: 'somalia_farmland_stats',
  fileFormat: 'CSV'
});

print('Export tasks created. Check the Tasks tab to run exports.');
