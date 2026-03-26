/************************************************************
 * SOMALIA SEASONAL FARMLAND DETECTION
 * Interactive Gu vs Deyr comparison using Sentinel-2
 *
 * Features:
 * - Somalia AOI = entire Somalia
 * - No study area box
 * - Dropdown filter for SWALIM-informed agricultural areas
 * - Analysis updates automatically when zone changes
 * - Gu vs Deyr seasonal farmland comparison
 *
 * NOTE:
 * These agricultural zones are SWALIM-informed approximate geometries,
 * not official downloaded SWALIM polygons.
 ************************************************************/

// ============================================================================
// 1. USER CONFIGURATION
// ============================================================================

// Entire Somalia extent
var somaliaAOI = ee.Geometry.Rectangle([40.5, -1.8, 51.5, 12.2]);

// Analysis year
var year = 2024;
var fallbackYear = 2023;

// Seasonal windows aligned to crop peak greenness
var guStart = year + '-05-01';
var guEnd   = year + '-06-30';

var deyrStart = year + '-11-01';
var deyrEnd   = year + '-12-31';

var guFallbackStart = fallbackYear + '-05-01';
var guFallbackEnd   = fallbackYear + '-06-30';

var deyrFallbackStart = fallbackYear + '-11-01';
var deyrFallbackEnd   = fallbackYear + '-12-31';

// ============================================================================
// 2. SWALIM-INFORMED AGRICULTURAL ZONES (APPROXIMATE)
// ============================================================================

// 1) Shabelle riverine agriculture
var lowerShabelleRiverine = ee.Geometry.Rectangle([43.2, 1.6, 45.3, 3.2]);
var middleShabelleRiverine = ee.Geometry.Rectangle([45.0, 2.1, 46.7, 4.7]);

// 2) Juba riverine agriculture
var lowerJubaRiverine = ee.Geometry.Rectangle([41.0, -0.3, 42.6, 1.3]);
var middleJubaRiverine = ee.Geometry.Rectangle([42.0, 0.3, 43.6, 2.3]);

// 3) Southern rain-fed belt
var bayAgropastoral = ee.Geometry.Rectangle([43.0, 2.3, 44.8, 4.2]);
var bakoolAgropastoral = ee.Geometry.Rectangle([43.6, 3.0, 45.2, 4.9]);
var gedoAgropastoral = ee.Geometry.Rectangle([41.8, 2.0, 43.8, 4.2]);

// 4) Northwestern rain-fed pockets
var gebileyPocket = ee.Geometry.Rectangle([43.3, 9.5, 44.6, 10.4]);
var hargeisaPocket = ee.Geometry.Rectangle([44.0, 9.3, 45.0, 10.1]);
var boramaPocket = ee.Geometry.Rectangle([42.6, 9.7, 43.6, 10.5]);

// Dropdown zone names
var zoneNames = [
  'Lower Shabelle Riverine',
  'Middle Shabelle Riverine',
  'Lower Juba Riverine',
  'Middle Juba Riverine',
  'Bay Agro-pastoral',
  'Bakool Agro-pastoral',
  'Gedo Agro-pastoral',
  'Gebiley Pocket',
  'Hargeisa Pocket',
  'Borama Pocket',
  'Shabelle Belt',
  'Juba Belt',
  'NW Rainfed Pockets',
  'Southern Rainfed Belt',
  'All SWALIM Agricultural Areas'
];

// ============================================================================
// 3. HELPER FUNCTIONS
// ============================================================================

function getSelectedZoneGeometry(zoneName) {
  if (zoneName === 'Lower Shabelle Riverine') return lowerShabelleRiverine;
  if (zoneName === 'Middle Shabelle Riverine') return middleShabelleRiverine;
  if (zoneName === 'Lower Juba Riverine') return lowerJubaRiverine;
  if (zoneName === 'Middle Juba Riverine') return middleJubaRiverine;

  if (zoneName === 'Bay Agro-pastoral') return bayAgropastoral;
  if (zoneName === 'Bakool Agro-pastoral') return bakoolAgropastoral;
  if (zoneName === 'Gedo Agro-pastoral') return gedoAgropastoral;

  if (zoneName === 'Gebiley Pocket') return gebileyPocket;
  if (zoneName === 'Hargeisa Pocket') return hargeisaPocket;
  if (zoneName === 'Borama Pocket') return boramaPocket;

  if (zoneName === 'Shabelle Belt') {
    return lowerShabelleRiverine.union(middleShabelleRiverine);
  }

  if (zoneName === 'Juba Belt') {
    return lowerJubaRiverine.union(middleJubaRiverine);
  }

  if (zoneName === 'NW Rainfed Pockets') {
    return gebileyPocket.union(hargeisaPocket).union(boramaPocket);
  }

  if (zoneName === 'Southern Rainfed Belt') {
    return bayAgropastoral.union(bakoolAgropastoral).union(gedoAgropastoral);
  }

  if (zoneName === 'All SWALIM Agricultural Areas') {
    return lowerShabelleRiverine
      .union(middleShabelleRiverine)
      .union(lowerJubaRiverine)
      .union(middleJubaRiverine)
      .union(bayAgropastoral)
      .union(bakoolAgropastoral)
      .union(gedoAgropastoral)
      .union(gebileyPocket)
      .union(hargeisaPocket)
      .union(boramaPocket);
  }

  return lowerShabelleRiverine;
}

function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}

function getSeasonCollection(filterArea, startDate, endDate, fallbackStart, fallbackEnd) {
  var primary = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(somaliaAOI)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

  var fallback = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(somaliaAOI)
    .filterDate(fallbackStart, fallbackEnd)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

  var primaryCount = primary.size();

  return ee.ImageCollection(
    ee.Algorithms.If(primaryCount.gt(0), primary, fallback)
  )
  .map(maskS2clouds)
  .map(function(img) {
    return img.clip(filterArea);
  });
}

function addIndices(composite) {
  var ndvi = composite.normalizedDifference(['B8', 'B4']).rename('NDVI');

  var evi = composite.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': composite.select('B8'),
      'RED': composite.select('B4'),
      'BLUE': composite.select('B2')
    }).rename('EVI');

  var savi = composite.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5', {
      'NIR': composite.select('B8'),
      'RED': composite.select('B4')
    }).rename('SAVI');

  var ndmi = composite.normalizedDifference(['B8', 'B11']).rename('NDMI');

  var bsi = composite.expression(
    '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
      'SWIR1': composite.select('B11'),
      'RED': composite.select('B4'),
      'NIR': composite.select('B8'),
      'BLUE': composite.select('B2')
    }).rename('BSI');

  return composite.addBands([ndvi, evi, savi, ndmi, bsi]);
}

function classifyFarmland(image, seasonName, filterArea) {
  var ndvi = image.select('NDVI');
  var savi = image.select('SAVI');
  var ndmi = image.select('NDMI');
  var bsi = image.select('BSI');

  var vegetated = ndvi.gt(0.25);
  var notForest = ndvi.lt(0.75);
  var soilAdjusted = savi.gt(0.18);
  var hasMoisture = ndmi.gt(-0.05);
  var notBare = bsi.lt(0.0);

  var farmland = vegetated
    .and(notForest)
    .and(soilAdjusted)
    .and(hasMoisture)
    .and(notBare)
    .rename(seasonName + '_farmland');

  var minAreaHa = 0.5;
  var minAreaPixels = (minAreaHa * 10000) / (10 * 10);

  var pixelCount = farmland.selfMask().connectedPixelCount({
    maxSize: 256,
    eightConnected: true
  });

  return farmland.updateMask(pixelCount.gte(minAreaPixels)).clip(filterArea);
}

function areaStats(maskImage, geometry, name) {
  var pixelArea = ee.Image.pixelArea();

  var area = maskImage.selfMask()
    .multiply(pixelArea)
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry,
      scale: 10,
      maxPixels: 1e13,
      bestEffort: true
    });

  var totalArea = pixelArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geometry,
    scale: 10,
    maxPixels: 1e13,
    bestEffort: true
  });

  var bandName = ee.String(maskImage.bandNames().get(0));
  var areaHa = ee.Number(area.get(bandName)).divide(10000);
  var totalHa = ee.Number(totalArea.get('area')).divide(10000);
  var percent = areaHa.divide(totalHa).multiply(100);

  print(name + ' area (ha):', areaHa);
  print(name + ' percent:', percent);
}

// ============================================================================
// 4. UI SETUP
// ============================================================================

ui.root.clear();

var map = ui.Map();
map.setOptions('SATELLITE');
ui.root.add(map);

var controlPanel = ui.Panel({
  style: {
    position: 'top-left',
    width: '340px',
    padding: '12px',
    backgroundColor: 'white'
  }
});

var legendPanel = ui.Panel({
  style: {
    position: 'bottom-left',
    width: '300px',
    padding: '12px 16px',
    backgroundColor: 'white',
    border: '2px solid #333333',
    borderRadius: '8px'
  }
});

map.add(legendPanel);

controlPanel.add(ui.Label({
  value: 'Somalia Seasonal Farmland Detection',
  style: {
    fontWeight: 'bold',
    fontSize: '16px',
    margin: '0 0 8px 0'
  }
}));

controlPanel.add(ui.Label({
  value: 'Select a SWALIM-informed agricultural area to run Gu vs Deyr farmland analysis.',
  style: {
    fontSize: '11px',
    color: 'gray',
    margin: '0 0 10px 0'
  }
}));

var zoneSelect = ui.Select({
  items: zoneNames,
  value: 'Shabelle Belt',
  placeholder: 'Select agricultural zone'
});

controlPanel.add(ui.Label('Agricultural Area', {
  fontWeight: 'bold',
  fontSize: '12px'
}));
controlPanel.add(zoneSelect);

var statusLabel = ui.Label('Ready.', {
  fontSize: '11px',
  color: 'gray',
  margin: '10px 0 0 0'
});
controlPanel.add(statusLabel);

ui.root.add(controlPanel);

// ============================================================================
// 5. LEGEND
// ============================================================================

function drawLegend(selectedZone) {
  legendPanel.clear();

  legendPanel.add(ui.Label({
    value: 'SEASONAL FARMLAND LEGEND',
    style: {
      fontWeight: 'bold',
      fontSize: '14px',
      margin: '0 0 10px 0',
      textAlign: 'center'
    }
  }));

  legendPanel.add(ui.Label({
    value: 'Filter Area: ' + selectedZone,
    style: {fontSize: '11px', margin: '0 0 8px 0'}
  }));

  function makeLegendRow(color, name, desc) {
    var row = ui.Panel({
      layout: ui.Panel.Layout.flow('horizontal'),
      style: {margin: '4px 0'}
    });

    row.add(ui.Label({
      value: '■',
      style: {color: color, fontSize: '16px', margin: '0 8px 0 0'}
    }));

    var text = ui.Panel({
      layout: ui.Panel.Layout.flow('vertical'),
      style: {margin: '0', padding: '0'}
    });

    text.add(ui.Label({
      value: name,
      style: {fontWeight: 'bold', fontSize: '11px', margin: '0'}
    }));

    text.add(ui.Label({
      value: desc,
      style: {fontSize: '9px', color: 'gray', margin: '0'}
    }));

    row.add(text);
    return row;
  }

  legendPanel.add(ui.Label({
    value: 'Seasonal Classes',
    style: {fontWeight: 'bold', fontSize: '11px', margin: '10px 0 4px 0'}
  }));

  legendPanel.add(makeLegendRow('#00cc44', 'Gu Farmland', 'Detected in Gu peak season'));
  legendPanel.add(makeLegendRow('#0066ff', 'Deyr Farmland', 'Detected in Deyr peak season'));
  legendPanel.add(makeLegendRow('#ffff00', 'Persistent Farmland', 'Detected in both seasons'));

  legendPanel.add(ui.Label({
    value: 'Peak Windows',
    style: {fontWeight: 'bold', fontSize: '11px', margin: '10px 0 4px 0'}
  }));

  legendPanel.add(ui.Label('Gu peak: May-June', {fontSize: '10px'}));
  legendPanel.add(ui.Label('Deyr peak: Nov-Dec', {fontSize: '10px'}));

  legendPanel.add(ui.Label({
    value: 'Thresholds',
    style: {fontWeight: 'bold', fontSize: '11px', margin: '10px 0 4px 0'}
  }));

  legendPanel.add(ui.Label('NDVI > 0.25', {fontSize: '10px'}));
  legendPanel.add(ui.Label('NDVI < 0.75', {fontSize: '10px'}));
  legendPanel.add(ui.Label('SAVI > 0.18', {fontSize: '10px'}));
  legendPanel.add(ui.Label('NDMI > -0.05', {fontSize: '10px'}));
  legendPanel.add(ui.Label('BSI < 0', {fontSize: '10px'}));
  legendPanel.add(ui.Label('Min area: 0.5 ha', {fontSize: '10px'}));

  legendPanel.add(ui.Label({
    value: 'Sentinel-2 | 10m | SWALIM-informed zones',
    style: {fontSize: '9px', color: 'gray', margin: '10px 0 0 0'}
  }));
}

// ============================================================================
// 6. MAIN ANALYSIS
// ============================================================================

function runAnalysis(selectedZone) {
  statusLabel.setValue('Running analysis for: ' + selectedZone + ' ...');

  var filterArea = getSelectedZoneGeometry(selectedZone);

  map.layers().reset();
  map.centerObject(filterArea, selectedZone === 'All SWALIM Agricultural Areas' ? 6 : 8);

  // Optional Somalia outline
  map.addLayer(
    ee.Image().paint(somaliaAOI, 1, 2),
    {palette: ['999999']},
    'Somalia AOI',
    false
  );

  // Selected filter area
  map.addLayer(
    ee.Image().paint(filterArea, 1, 2),
    {palette: ['yellow']},
    'Selected Agricultural Filter Area',
    true
  );

  var guCollection = getSeasonCollection(
    filterArea, guStart, guEnd, guFallbackStart, guFallbackEnd
  );

  var deyrCollection = getSeasonCollection(
    filterArea, deyrStart, deyrEnd, deyrFallbackStart, deyrFallbackEnd
  );

  var guComposite = guCollection.median().clip(filterArea);
  var deyrComposite = deyrCollection.median().clip(filterArea);

  var guImage = addIndices(guComposite);
  var deyrImage = addIndices(deyrComposite);

  var guFarmland = classifyFarmland(guImage, 'Gu', filterArea);
  var deyrFarmland = classifyFarmland(deyrImage, 'Deyr', filterArea);

  var persistentFarmland = guFarmland.and(deyrFarmland).rename('Persistent_Farmland');
  var guOnly = guFarmland.and(deyrFarmland.not()).rename('Gu_Only');
  var deyrOnly = deyrFarmland.and(guFarmland.not()).rename('Deyr_Only');
  var anyFarmland = guFarmland.or(deyrFarmland).rename('Any_Farmland');

  var ndviDiff = guImage.select('NDVI')
    .subtract(deyrImage.select('NDVI'))
    .rename('NDVI_Gu_minus_Deyr')
    .clip(filterArea);

  // Basemap layers
  map.addLayer(
    guComposite.select(['B4', 'B3', 'B2']),
    {min: 0, max: 0.3},
    'Gu Natural Color',
    false
  );

  map.addLayer(
    deyrComposite.select(['B4', 'B3', 'B2']),
    {min: 0, max: 0.3},
    'Deyr Natural Color',
    false
  );

  map.addLayer(
    guImage.select('NDVI'),
    {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
    'Gu NDVI Peak',
    false
  );

  map.addLayer(
    deyrImage.select('NDVI'),
    {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
    'Deyr NDVI Peak',
    false
  );

  // Farmland layers
  map.addLayer(
    guFarmland.selfMask(),
    {palette: ['00cc44']},
    'Gu Farmland',
    true
  );

  map.addLayer(
    deyrFarmland.selfMask(),
    {palette: ['0066ff']},
    'Deyr Farmland',
    true
  );

  map.addLayer(
    persistentFarmland.selfMask(),
    {palette: ['ffff00']},
    'Persistent Farmland (Gu + Deyr)',
    true
  );

  map.addLayer(
    guOnly.selfMask(),
    {palette: ['00cc44']},
    'Gu Only Farmland',
    false
  );

  map.addLayer(
    deyrOnly.selfMask(),
    {palette: ['0066ff']},
    'Deyr Only Farmland',
    false
  );

  map.addLayer(
    anyFarmland.selfMask(),
    {palette: ['limegreen']},
    'Any Seasonal Farmland',
    false
  );

  map.addLayer(
    ndviDiff,
    {min: -0.4, max: 0.4, palette: ['blue', 'white', 'green']},
    'NDVI Difference (Gu - Deyr)',
    false
  );

  drawLegend(selectedZone);

  print('==================================================');
  print('SEASONAL IMAGE AVAILABILITY');
  print('==================================================');
  print('Selected zone:', selectedZone);
  print('Gu scenes:', guCollection.size());
  print('Deyr scenes:', deyrCollection.size());

  print('==================================================');
  print('SEASONAL FARMLAND STATISTICS');
  print('==================================================');
  print('Somalia AOI: Entire Somalia');
  print('Selected agricultural filter area:', selectedZone);

  areaStats(guFarmland, filterArea, 'Gu farmland');
  areaStats(deyrFarmland, filterArea, 'Deyr farmland');
  areaStats(persistentFarmland, filterArea, 'Persistent farmland');
  areaStats(anyFarmland, filterArea, 'Any seasonal farmland');

  statusLabel.setValue('Analysis complete for: ' + selectedZone);
}

// ============================================================================
// 7. DROPDOWN EVENT
// ============================================================================

zoneSelect.onChange(function(value) {
  runAnalysis(value);
});

// ============================================================================
// 8. INITIAL RUN
// ============================================================================

runAnalysis(zoneSelect.getValue());