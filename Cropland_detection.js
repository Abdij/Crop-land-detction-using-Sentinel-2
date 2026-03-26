/************************************************************
 * FARMLAND DETECTION - SOMALIA AGRICULTURAL AREAS
 * Using Sentinel-2 Only
 * 
 * Purpose: Identify and visualize agricultural areas across
 * major farming regions in Somalia using Sentinel-2 imagery
 ************************************************************/

// ============================================================================
// 1. USER CONFIGURATION
// ============================================================================

// Define Somalia's major agricultural regions
var somaliaAOI = ee.Geometry.Rectangle([41.0, -1.7, 51.0, 11.5]);

// Specific agricultural zones
var lowerShabelle = ee.Geometry.Rectangle([43.5, 1.5, 45.5, 3.0]);
var middleShabelle = ee.Geometry.Rectangle([45.0, 2.0, 46.5, 4.5]);
var bayRegion = ee.Geometry.Rectangle([43.0, 2.5, 44.5, 4.0]);
var jubaValley = ee.Geometry.Rectangle([41.5, -0.5, 43.0, 2.0]);
var gabiley = ee.Geometry.Rectangle([43.5, 9.5, 44.5, 10.5]);

// Select which region to analyze (change as needed)
var studyArea = lowerShabelle; // Change to any of the above

// Center map on selected area
Map.centerObject(studyArea, 10);
Map.addLayer(ee.Image().paint(studyArea, 1, 2), {palette: ['yellow']}, 'Study Area', true);

// ============================================================================
// 2. LOAD SENTINEL-2 IMAGERY (GUARANTEED TO WORK)
// ============================================================================

print('==================================================');
print('LOADING SENTINEL-2 IMAGERY');
print('==================================================');

// Load Sentinel-2 imagery for the main growing seasons
var sentinel2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2025-01-01', '2025-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

var sceneCount = sentinel2.size();
print('Sentinel-2 scenes available:', sceneCount);

// If no scenes found, try 2023
var sentinel2_2024 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2024-01-01', '2024-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

var sceneCount2024 = sentinel2_2024.size();
print('Sentinel-2 scenes (2024):', sceneCount2024);

// Use the collection with more scenes
var useCollection = ee.ImageCollection(
  ee.Algorithms.If(
    sceneCount.gt(0),
    sentinel2,
    sentinel2_2024
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

// Bare Soil Index
var bsi = composite.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
    'SWIR1': composite.select('B11'),
    'RED': composite.select('B4'),
    'NIR': composite.select('B8'),
    'BLUE': composite.select('B2')
}).rename('BSI');

// ============================================================================
// 5. FARMLAND CLASSIFICATION USING MULTIPLE INDICES
// ============================================================================

// Create composite image with all indices
var indices = ee.Image.cat([ndvi, evi, savi, ndmi, bsi]);

// Step 1: Identify vegetated areas (NDVI > 0.2 - adjusted for arid regions)
var vegetated = ndvi.gt(0.2).rename('vegetated');

// Step 2: Filter out very dense vegetation (forests) - NDVI < 0.7
var notForest = ndvi.lt(0.7).rename('notForest');

// Step 3: Soil-adjusted filter for arid regions
var soilAdjusted = savi.gt(0.15).rename('soilAdjusted');

// Step 4: Moisture filter (irrigated crops have higher moisture)
var hasMoisture = ndmi.gt(-0.1).rename('hasMoisture');

// Step 5: Not bare soil
var notBare = bsi.lt(0).rename('notBare');

// Combine all criteria for farmland detection
var farmland = vegetated
  .and(notForest)
  .and(soilAdjusted)
  .and(hasMoisture)
  .and(notBare)
  .rename('farmland');

// ============================================================================
// 6. APPLY AREA FILTERING (REMOVE ISOLATED PIXELS)
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
  false
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
  savi,
  {min: 0, max: 0.6, palette: ['brown', 'yellow', 'green']},
  'SAVI (Soil Adjusted)',
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
  filteredFarmland.selfMask(),
  {palette: ['limegreen']},
  'Farmland Areas',
  true
);

// ============================================================================
// 9. CREATE LEGEND (FIXED - REMOVED boxShadow)
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
  value: 'FARMLAND DETECTION LEGEND',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    color: '#2c3e50',
    margin: '0 0 10px 0',
    textAlign: 'center'
  }
});
legend.add(title);

// Study Area
var studyAreaLabel = ui.Label({
  value: 'Study Area: ' + 
    (studyArea == lowerShabelle ? 'Lower Shabelle' :
     studyArea == middleShabelle ? 'Middle Shabelle' :
     studyArea == bayRegion ? 'Bay Region' :
     studyArea == jubaValley ? 'Juba Valley' : 'Gabiley'),
  style: {fontSize: '11px', color: '#34495e', margin: '0 0 8px 0'}
});
legend.add(studyAreaLabel);

// RGB Composites Section
var rgbTitle = ui.Label({
  value: 'RGB COMPOSITES',
  style: {fontWeight: 'bold', fontSize: '11px', color: '#34495e', margin: '8px 0 4px 0'}
});
legend.add(rgbTitle);

var rgbItems = [
  {color: '#ffffff', name: 'Natural Color (4-3-2)', desc: 'True color'},
  {color: '#ffaaaa', name: 'False Color (8-4-3)', desc: 'Vegetation=Red'},
  {color: '#aaffaa', name: 'Agriculture (11-8-4)', desc: 'Highlights crops'}
];

rgbItems.forEach(function(item) {
  var row = ui.Panel({
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {margin: '2px 0'}
  });
  row.add(ui.Label('●', {color: item.color, fontSize: '14px', margin: '0 4px 0 0'}));
  row.add(ui.Label(item.name, {fontSize: '10px', margin: '0 4px 0 0'}));
  row.add(ui.Label('(' + item.desc + ')', {fontSize: '9px', color: '#7f8c8d'}));
  legend.add(row);
});

// Indices Section
var indicesTitle = ui.Label({
  value: 'VEGETATION INDICES',
  style: {fontWeight: 'bold', fontSize: '11px', color: '#34495e', margin: '12px 0 4px 0'}
});
legend.add(indicesTitle);

var indexItems = [
  {name: 'NDVI', palette: 'brown→yellow→green', desc: 'Vegetation health'},
  {name: 'SAVI', palette: 'brown→yellow→green', desc: 'Soil-adjusted'},
  {name: 'NDMI', palette: 'brown→blue', desc: 'Moisture content'}
];

indexItems.forEach(function(item) {
  var row = ui.Panel({
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {margin: '2px 0'}
  });
  row.add(ui.Label('•', {fontSize: '12px', margin: '0 4px 0 0'}));
  row.add(ui.Label(item.name + ':', {fontWeight: 'bold', fontSize: '10px', margin: '0 4px 0 0'}));
  row.add(ui.Label(item.palette, {fontSize: '9px', color: '#7f8c8d', margin: '0 4px 0 0'}));
  row.add(ui.Label('(' + item.desc + ')', {fontSize: '9px', color: '#7f8c8d'}));
  legend.add(row);
});

// Farmland Classes Section
var classesTitle = ui.Label({
  value: 'FARMLAND CLASSES',
  style: {fontWeight: 'bold', fontSize: '11px', color: '#34495e', margin: '12px 0 4px 0'}
});
legend.add(classesTitle);

var classItems = [
  {color: '#32cd32', name: 'Farmland Areas', desc: 'Cultivated crops'},
  {color: '#ffff00', name: 'Potential Farmland', desc: 'Low confidence'},
  {color: '#ff0000', name: 'Non-Farmland', desc: 'Urban/bare soil/water'}
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
  '• SAVI > 0.15 (soil-adjusted)',
  '• NDMI > -0.1 (has moisture)',
  '• BSI < 0 (not bare soil)',
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
footer.add(ui.Label('For agricultural planning only', {fontSize: '9px', color: '#e74c3c', margin: '4px 0 0 0'}));

legend.add(footer);

// Add legend to map
Map.add(legend);

// ============================================================================
// 10. CALCULATE STATISTICS
// ============================================================================

var pixelArea = ee.Image.pixelArea();

// Calculate farmland area
var farmlandArea = filteredFarmland.selfMask()
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

print('==================================================');
print('FARMLAND STATISTICS');
print('==================================================');
print('Study area:', studyArea == lowerShabelle ? 'Lower Shabelle' :
      studyArea == middleShabelle ? 'Middle Shabelle' :
      studyArea == bayRegion ? 'Bay Region' :
      studyArea == jubaValley ? 'Juba Valley' : 'Gabiley');
print('Total area (hectares):', totalHa);
print('Farmland detected (hectares):', farmlandHa);
print('Percentage farmland:', percentFarmland);
print('Number of Sentinel-2 scenes:', sceneCount);

// ============================================================================
// 11. ZONE COMPARISON (if using full Somalia)
// ============================================================================

function compareAllZones() {
  var zones = [
    {geom: lowerShabelle, name: 'Lower Shabelle'},
    {geom: middleShabelle, name: 'Middle Shabelle'},
    {geom: bayRegion, name: 'Bay Region'},
    {geom: jubaValley, name: 'Juba Valley'},
    {geom: gabiley, name: 'Gabiley'}
  ];
  
  print('==================================================');
  print('ALL AGRICULTURAL ZONES COMPARISON');
  print('==================================================');
  
  zones.forEach(function(z) {
    var zoneFarmland = filteredFarmland.selfMask()
      .multiply(pixelArea)
      .reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry: z.geom,
        scale: 10,
        bestEffort: true
      });
    
    var zoneArea = pixelArea.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: z.geom,
      scale: 10,
      bestEffort: true
    });
    
    zoneFarmland.evaluate(function(farm) {
      zoneArea.evaluate(function(total) {
        var farmHa = (farm && farm.farmland) ? farm.farmland / 10000 : 0;
        var totalHa = (total && total.area) ? total.area / 10000 : 0;
        var percent = totalHa > 0 ? (farmHa / totalHa * 100).toFixed(1) : 0;
        print(z.name + ': ' + Math.round(farmHa) + ' ha (' + percent + '%)');
      });
    });
  });
}

// Uncomment to compare all zones
// compareAllZones();

// ============================================================================
// 12. EXPORT
// ============================================================================

Export.image.toDrive({
  image: filteredFarmland.selfMask(),
  description: 'Somalia_Farmland_Sentinel2',
  folder: 'GEE_Exports',
  fileNamePrefix: 'somalia_farmland_s2',
  region: studyArea,
  scale: 10,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: ndvi,
  description: 'Somalia_NDVI',
  folder: 'GEE_Exports',
  fileNamePrefix: 'somalia_ndvi',
  region: studyArea,
  scale: 10,
  maxPixels: 1e13
});

print('==================================================');
print('FARMLAND DETECTION COMPLETE');
print('==================================================');
print('Legend added to map (bottom-right)');
print('Green areas: Detected farmland');
print('Check Console for statistics');
print('Export tasks available in Tasks tab');