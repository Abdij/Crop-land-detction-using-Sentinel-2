/************************************************************
 * FARMLAND DETECTION - SOMALIA AGRICULTURAL AREAS
 * Using Sentinel-2 Only with Somalia Boundary
 * 
 * Purpose: Identify and visualize agricultural areas across
 * Somalia using Sentinel-2 imagery and official boundary
 ************************************************************/

// ============================================================================
// 1. LOAD SOMALIA BOUNDARY
// ============================================================================

// Load Somalia country boundary from FAO GAUL dataset
var somalia = ee.FeatureCollection('FAO/GAUL/2015/level0')
  .filter(ee.Filter.eq('ADM0_NAME', 'Somalia'));

// Alternative: Load from Large Scale International Boundary (LSIB)
var somaliaLSIB = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
  .filter(ee.Filter.eq('country_na', 'Somalia'));

// Use GAUL boundary as primary (more detailed)
var somaliaBoundary = somalia;

// Check if boundary loaded properly
print('Somalia Boundary:', somaliaBoundary);
Map.centerObject(somaliaBoundary, 6);

// Add boundary to map
Map.addLayer(somaliaBoundary, {color: 'red', width: 2}, 'Somalia Boundary', false);

// Define study area as Somalia boundary
var studyArea = somaliaBoundary;

// ============================================================================
// 2. LOAD SENTINEL-2 IMAGERY
// ============================================================================

print('==================================================');
print('LOADING SENTINEL-2 IMAGERY FOR SOMALIA');
print('==================================================');

// Load Sentinel-2 imagery for the main growing seasons
var sentinel2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2025-01-01', '2025-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

var sceneCount = sentinel2.size();
print('Sentinel-2 scenes (2025):', sceneCount);

// If no scenes found, try 2023
var sentinel2_2023 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2024-01-01', '2024-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

var sceneCount2023 = sentinel2_2023.size();
print('Sentinel-2 scenes (2024):', sceneCount2023);

// Use the collection with more scenes
var useCollection = ee.ImageCollection(
  ee.Algorithms.If(
    sceneCount.gt(0),
    sentinel2,
    sentinel2_2023
  )
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
  return image.updateMask(mask).divide(10000);
}

var processedS2 = useCollection.map(maskS2clouds);

// Create median composite for Somalia
var composite = processedS2.median().clip(studyArea);

// ============================================================================
// 4. CALCULATE VEGETATION INDICES FOR FARMLAND DETECTION
// ============================================================================

// Normalized Difference Vegetation Index (NDVI)
var ndvi = composite.normalizedDifference(['B8', 'B4']).rename('NDVI');

// Enhanced Vegetation Index (EVI) - better for dense vegetation
var evi = composite.expression(
  '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
    'NIR': composite.select('B8'),
    'RED': composite.select('B4'),
    'BLUE': composite.select('B2')
}).rename('EVI');

// Soil-Adjusted Vegetation Index (SAVI) - better for arid regions like Somalia
var savi = composite.expression(
  '((NIR - RED) / (NIR + RED + 0.5)) * 1.5', {
    'NIR': composite.select('B8'),
    'RED': composite.select('B4')
}).rename('SAVI');

// Normalized Difference Moisture Index (NDMI) - for irrigation detection
var ndmi = composite.normalizedDifference(['B8', 'B11']).rename('NDMI');

// Bare Soil Index (BSI)
var bsi = composite.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
    'SWIR1': composite.select('B11'),
    'RED': composite.select('B4'),
    'NIR': composite.select('B8'),
    'BLUE': composite.select('B2')
}).rename('BSI');

// Modified Soil Adjusted Vegetation Index (MSAVI) - even better for arid regions
var msavi = composite.expression(
  '(2 * NIR + 1 - sqrt(pow((2 * NIR + 1), 2) - 8 * (NIR - RED))) / 2', {
    'NIR': composite.select('B8'),
    'RED': composite.select('B4')
}).rename('MSAVI');

// ============================================================================
// 5. FARMLAND CLASSIFICATION USING MULTIPLE INDICES
// ============================================================================

// Step 1: Identify vegetated areas (NDVI > 0.2 - adjusted for arid Somalia)
var vegetated = ndvi.gt(0.2).rename('vegetated');

// Step 2: Filter out very dense vegetation (natural forests) - NDVI < 0.7
var notForest = ndvi.lt(0.7).rename('notForest');

// Step 3: Soil-adjusted filter for arid regions (using MSAVI for better results)
var soilAdjusted = msavi.gt(0.1).rename('soilAdjusted');

// Step 4: Moisture filter (irrigated crops have higher moisture)
var hasMoisture = ndmi.gt(-0.1).rename('hasMoisture');

// Step 5: Not bare soil
var notBare = bsi.lt(0).rename('notBare');

// Step 6: Additional check for agricultural areas (EVI > 0.1)
var agriVegetation = evi.gt(0.1).rename('agriVegetation');

// Combine all criteria for farmland detection
var farmland = vegetated
  .and(notForest)
  .and(soilAdjusted)
  .and(hasMoisture)
  .and(notBare)
  .and(agriVegetation)
  .rename('farmland');

// ============================================================================
// 6. APPLY AREA FILTERING AND SMOOTHING
// ============================================================================

var minAreaHa = 0.5; // Minimum 0.5 hectare to be considered farmland
var minAreaPixels = (minAreaHa * 10000) / (10 * 10); // ~50 pixels at 10m

// Count connected pixels
var pixelCount = farmland.selfMask()
  .connectedPixelCount({
    maxSize: 256,
    eightConnected: true
  });

// Filter out small patches
var filteredFarmland = farmland.updateMask(pixelCount.gte(minAreaPixels));

// Apply focal mode for smoothing (remove speckle)
var smoothedFarmland = filteredFarmland.focal_mode({
  radius: 1,
  kernelType: 'square',
  iterations: 1
});

// ============================================================================
// 7. CREATE RGB AND FALSE COLOR COMPOSITES
// ============================================================================

// Natural color composite (RGB: B4, B3, B2)
var naturalColor = composite.select(['B4', 'B3', 'B2']).rename(['red', 'green', 'blue']);

// False color (vegetation in red - NIR, Red, Green)
var falseColor = composite.select(['B8', 'B4', 'B3']).rename(['nir', 'red', 'green']);

// Agriculture composite (SWIR, NIR, Red) - highlights crops
var agriColor = composite.select(['B11', 'B8', 'B4']).rename(['swir', 'nir', 'red']);

// ============================================================================
// 8. VISUALIZATION
// ============================================================================

// Add background composites
Map.addLayer(
  naturalColor,
  {min: 0, max: 0.3},
  'Natural Color (RGB)',
  true
);

Map.addLayer(
  falseColor,
  {min: 0, max: 0.5},
  'False Color (Vegetation=Red)',
  false
);

Map.addLayer(
  agriColor,
  {min: 0, max: 0.5},
  'Agriculture Composite',
  false
);

// Add individual indices
Map.addLayer(
  ndvi,
  {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
  'NDVI',
  false
);

Map.addLayer(
  msavi,
  {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
  'MSAVI (Modified Soil Adjusted)',
  false
);

Map.addLayer(
  ndmi,
  {min: -0.3, max: 0.4, palette: ['brown', 'yellow', 'green', 'blue']},
  'NDMI (Moisture)',
  false
);

// Add final farmland layer (VISIBLE BY DEFAULT)
Map.addLayer(
  smoothedFarmland.selfMask(),
  {palette: ['limegreen'], opacity: 0.7},
  'Farmland Areas',
  true
);

// ============================================================================
// 9. CREATE LEGEND
// ============================================================================

// Create a legend panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '12px 16px',
    backgroundColor: 'white',
    border: '2px solid #333333',
    borderRadius: '8px',
    fontFamily: 'Arial, sans-serif',
    width: '260px'
  }
});

// Title
var title = ui.Label({
  value: 'SOMALIA FARMLAND DETECTION',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    color: '#2c3e50',
    margin: '0 0 10px 0',
    textAlign: 'center'
  }
});
legend.add(title);

// Farmland Classes Section
var classesTitle = ui.Label({
  value: 'FARMLAND CLASSES',
  style: {fontWeight: 'bold', fontSize: '11px', color: '#34495e', margin: '12px 0 4px 0'}
});
legend.add(classesTitle);

var classItems = [
  {color: '#32cd32', name: 'Farmland Areas', desc: 'Cultivated crops detected'},
  {color: '#ffff00', name: 'Potential Farmland', desc: 'Low confidence areas'},
  {color: '#ff0000', name: 'Non-Farmland', desc: 'Urban/bare soil/water/forest'}
];

classItems.forEach(function(item) {
  var row = ui.Panel({
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {margin: '4px 0', padding: '2px 4px', backgroundColor: '#f8f9fa', borderRadius: '4px'}
  });
  
  var colorBox = ui.Label({
    value: '●',
    style: {color: item.color, fontSize: '18px', margin: '0 8px 0 0', padding: '0'}
  });
  
  var textPanel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style: {margin: '0', padding: '0'}
  });
  
  var className = ui.Label({
    value: item.name,
    style: {fontWeight: 'bold', fontSize: '11px', color: '#2c3e50', margin: '0', padding: '0'}
  });
  
  var classDesc = ui.Label({
    value: item.desc,
    style: {fontSize: '9px', color: '#7f8c8d', margin: '0', padding: '0'}
  });
  
  textPanel.add(className);
  textPanel.add(classDesc);
  row.add(colorBox);
  row.add(textPanel);
  legend.add(row);
});

// Detection Criteria
var criteriaTitle = ui.Label({
  value: 'DETECTION CRITERIA',
  style: {fontWeight: 'bold', fontSize: '11px', color: '#34495e', margin: '12px 0 4px 0'}
});
legend.add(criteriaTitle);

var criteria = [
  '• NDVI > 0.2 (vegetated)',
  '• NDVI < 0.7 (not forest)',
  '• MSAVI > 0.1 (soil-adjusted)',
  '• NDMI > -0.1 (has moisture)',
  '• BSI < 0 (not bare soil)',
  '• EVI > 0.1 (agricultural vegetation)',
  '• Min area: 0.5 hectares'
];

criteria.forEach(function(c) {
  legend.add(ui.Label(c, {fontSize: '9px', color: '#34495e', margin: '2px 0'}));
});

// Footer with date and data source
var footer = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {margin: '12px 0 0 0', padding: '8px 0 0 0'}
});

footer.add(ui.Label('Data: 2023-2024 Sentinel-2', {fontSize: '9px', color: '#7f8c8d'}));
footer.add(ui.Label('Resolution: 10m', {fontSize: '9px', color: '#7f8c8d'}));
footer.add(ui.Label('Boundary: FAO GAUL 2015', {fontSize: '9px', color: '#7f8c8d'}));
footer.add(ui.Label('For agricultural planning only', {fontSize: '9px', color: '#e74c3c', margin: '4px 0 0 0'}));

legend.add(footer);

// Add legend to map
Map.add(legend);

// ============================================================================
// 10. CALCULATE FARMLAND STATISTICS FOR SOMALIA
// ============================================================================

var pixelArea = ee.Image.pixelArea();

// Calculate farmland area within Somalia boundary
var farmlandArea = smoothedFarmland.selfMask()
  .multiply(pixelArea)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: studyArea,
    scale: 10,
    maxPixels: 1e13,
    bestEffort: true
  });

var totalArea = pixelArea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: studyArea,
  scale: 10,
  bestEffort: true
});

var farmlandHa = ee.Number(farmlandArea.get('farmland')).divide(10000);
var totalHa = ee.Number(totalArea.get('area')).divide(10000);
var percentFarmland = farmlandHa.divide(totalHa).multiply(100);

// Print statistics
print('==================================================');
print('SOMALIA FARMLAND STATISTICS');
print('==================================================');
print('Total Somalia Area (ha):', totalHa);
print('Farmland Area (ha):', farmlandHa);
print('Percentage Farmland:', percentFarmland);

// Create statistics panel
var statsPanel = ui.Panel({
  style: {
    position: 'top-right',
    padding: '12px 16px',
    backgroundColor: 'white',
    border: '2px solid #333333',
    borderRadius: '8px',
    fontFamily: 'Arial, sans-serif',
    width: '240px'
  }
});

statsPanel.add(ui.Label({
  value: 'FARMLAND STATISTICS',
  style: {fontWeight: 'bold', fontSize: '13px', color: '#2c3e50', margin: '0 0 8px 0'}
}));

statsPanel.add(ui.Label({
  value: 'Total Somalia Area: ' + totalHa.format('%.0f') + ' ha',
  style: {fontSize: '10px', margin: '4px 0'}
}));

statsPanel.add(ui.Label({
  value: 'Farmland Detected: ' + farmlandHa.format('%.0f') + ' ha',
  style: {fontSize: '10px', margin: '4px 0', fontWeight: 'bold', color: '#32cd32'}
}));

statsPanel.add(ui.Label({
  value: 'Percentage: ' + percentFarmland.format('%.1f') + '%',
  style: {fontSize: '10px', margin: '4px 0', fontWeight: 'bold'}
}));

statsPanel.add(ui.Label({
  value: 'Sentinel-2 Scenes: ' + sceneCount,
  style: {fontSize: '9px', color: '#7f8c8d', margin: '8px 0 2px 0'}
}));

statsPanel.add(ui.Label({
  value: 'Resolution: 10m',
  style: {fontSize: '9px', color: '#7f8c8d', margin: '2px 0'}
}));

Map.add(statsPanel);

// ============================================================================
// 11. EXPORT FARMLAND MAP
// ============================================================================

Export.image.toDrive({
  image: smoothedFarmland.selfMask(),
  description: 'Somalia_Farmland_Sentinel2',
  folder: 'GEE_Exports',
  fileNamePrefix: 'somalia_farmland_2024',
  region: studyArea,
  scale: 10,
  maxPixels: 1e13
});

// Optional: Export statistics as CSV
var statsFeature = ee.Feature(null, {
  'total_area_ha': totalHa,
  'farmland_area_ha': farmlandHa,
  'farmland_percentage': percentFarmland,
  'sentinel2_scenes': sceneCount,
  'date_processed': ee.Date(Date.now()).format()
});

Export.table.toDrive({
  collection: ee.FeatureCollection([statsFeature]),
  description: 'Somalia_Farmland_Statistics',
  folder: 'GEE_Exports',
  fileNamePrefix: 'somalia_farmland_stats',
  fileFormat: 'CSV'
});

print('Export tasks created. Check Tasks tab to run exports.');
