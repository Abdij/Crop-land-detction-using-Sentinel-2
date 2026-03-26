/************************************************************
 * SOMALIA PHENOLOGY-BASED SEASONAL FARMLAND DETECTION
 * Interactive Gu vs Deyr comparison using Sentinel-2
 *
 * Features:
 * - Somalia AOI = entire Somalia
 * - No study area box
 * - Dropdown filter for SWALIM-informed agricultural areas
 * - Dynamic Gu and Deyr date sliders
 * - 5-day date step aligned to Sentinel-2 revisit
 * - Phenology-based seasonal farmland detection
 *
 * NOTE:
 * These agricultural zones are SWALIM-informed approximate
 * geometries, not official downloaded SWALIM polygons.
 ************************************************************/

// ============================================================================
// 1. SOMALIA AOI
// ============================================================================

var somaliaAOI = ee.Geometry.Rectangle([40.5, -1.8, 51.5, 12.2]);

// ============================================================================
// 2. SWALIM-INFORMED AGRICULTURAL ZONES (APPROXIMATE)
// ============================================================================

// Shabelle riverine agriculture
var lowerShabelleRiverine = ee.Geometry.Rectangle([43.2, 1.6, 45.3, 3.2]);
var middleShabelleRiverine = ee.Geometry.Rectangle([45.0, 2.1, 46.7, 4.7]);

// Juba riverine agriculture
var lowerJubaRiverine = ee.Geometry.Rectangle([41.0, -0.3, 42.6, 1.3]);
var middleJubaRiverine = ee.Geometry.Rectangle([42.0, 0.3, 43.6, 2.3]);

// Southern rain-fed belt
var bayAgropastoral = ee.Geometry.Rectangle([43.0, 2.3, 44.8, 4.2]);
var bakoolAgropastoral = ee.Geometry.Rectangle([43.6, 3.0, 45.2, 4.9]);
var gedoAgropastoral = ee.Geometry.Rectangle([41.8, 2.0, 43.8, 4.2]);

// Northwestern rain-fed pockets
var gebileyPocket = ee.Geometry.Rectangle([43.3, 9.5, 44.6, 10.4]);
var hargeisaPocket = ee.Geometry.Rectangle([44.0, 9.3, 45.0, 10.1]);
var boramaPocket = ee.Geometry.Rectangle([42.6, 9.7, 43.6, 10.5]);

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
// 3. HELPERS
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

function getDateRangeFromSlider(slider) {
  var range = slider.getValue();
  return {
    start: range[0],
    end: range[1]
  };
}

function formatDateForDisplay(d) {
  return ee.Date(d).format('YYYY-MM-dd');
}

function updateDateLabel(labelWidget, prefix, rangeObj) {
  var start = formatDateForDisplay(rangeObj.start);
  var end = formatDateForDisplay(rangeObj.end);

  start.evaluate(function(s) {
    end.evaluate(function(e) {
      labelWidget.setValue(prefix + ': ' + s + ' to ' + e);
    });
  });
}

function updateDateLabel(labelWidget, prefix, rangeObj) {
  labelWidget.setValue(
    prefix + ': ' + formatJsDate(rangeObj.start) + ' to ' + formatJsDate(rangeObj.end)
  );
}

function maskS2(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  var qaMask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  // Scene Classification Layer mask
  var scl = image.select('SCL');
  var sclMask = scl.neq(3)   // cloud shadow
    .and(scl.neq(8))         // cloud medium probability
    .and(scl.neq(9))         // cloud high probability
    .and(scl.neq(10))        // cirrus
    .and(scl.neq(11));       // snow/ice

  return image
    .updateMask(qaMask)
    .updateMask(sclMask)
    .divide(10000)
    .copyProperties(image, ['system:time_start']);
}

function getCollection(filterArea, startDate, endDate) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(filterArea)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .map(maskS2)
    .map(function(img) {
      return img.clip(filterArea);
    });
}

function safeBaseImage(filterArea) {
  return ee.Image.constant([0, 0, 0, 0])
    .rename(['B2', 'B4', 'B8', 'B11'])
    .clip(filterArea)
    .toFloat();
}

function getMedianComposite(collection, filterArea) {
  return ee.Image(
    ee.Algorithms.If(
      collection.size().gt(0),
      collection.median().select(['B2', 'B4', 'B8', 'B11']).clip(filterArea),
      safeBaseImage(filterArea)
    )
  );
}

function makeNdvi(image) {
  return image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    .copyProperties(image, ['system:time_start']);
}

function makeNdmi(image) {
  return image.normalizedDifference(['B8', 'B11']).rename('NDMI')
    .copyProperties(image, ['system:time_start']);
}

function makeSavi(image) {
  return image.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5', {
      'NIR': image.select('B8'),
      'RED': image.select('B4')
    }).rename('SAVI')
    .copyProperties(image, ['system:time_start']);
}

function makeBsi(image) {
  return image.expression(
    '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
      'SWIR1': image.select('B11'),
      'RED': image.select('B4'),
      'NIR': image.select('B8'),
      'BLUE': image.select('B2')
    }).rename('BSI')
    .copyProperties(image, ['system:time_start']);
}

function buildPhenologyMetrics(collection, prefix, filterArea) {
  var ndviCol = collection.map(makeNdvi);
  var ndmiCol = collection.map(makeNdmi);
  var saviCol = collection.map(makeSavi);
  var bsiCol = collection.map(makeBsi);

  var hasData = collection.size().gt(0);

  var ndviP90 = ee.Image(ee.Algorithms.If(
    hasData,
    ndviCol.reduce(ee.Reducer.percentile([90])).rename(prefix + '_NDVI_P90'),
    ee.Image.constant(0).rename(prefix + '_NDVI_P90').clip(filterArea)
  ));

  var ndviMin = ee.Image(ee.Algorithms.If(
    hasData,
    ndviCol.min().rename(prefix + '_NDVI_MIN'),
    ee.Image.constant(0).rename(prefix + '_NDVI_MIN').clip(filterArea)
  ));

  var ndviMax = ee.Image(ee.Algorithms.If(
    hasData,
    ndviCol.max().rename(prefix + '_NDVI_MAX'),
    ee.Image.constant(0).rename(prefix + '_NDVI_MAX').clip(filterArea)
  ));

  var ndviAmp = ndviMax.subtract(ndviMin).rename(prefix + '_NDVI_AMP');

  var saviMax = ee.Image(ee.Algorithms.If(
    hasData,
    saviCol.max().rename(prefix + '_SAVI_MAX'),
    ee.Image.constant(0).rename(prefix + '_SAVI_MAX').clip(filterArea)
  ));

  var ndmiMean = ee.Image(ee.Algorithms.If(
    hasData,
    ndmiCol.mean().rename(prefix + '_NDMI_MEAN'),
    ee.Image.constant(0).rename(prefix + '_NDMI_MEAN').clip(filterArea)
  ));

  var bsiMed = ee.Image(ee.Algorithms.If(
    hasData,
    bsiCol.median().rename(prefix + '_BSI_MED'),
    ee.Image.constant(0).rename(prefix + '_BSI_MED').clip(filterArea)
  ));

  return ee.Image.cat([
    ndviP90, ndviMin, ndviMax, ndviAmp, saviMax, ndmiMean, bsiMed
  ]).clip(filterArea);
}

function classifyPhenologyFarmland(metrics, prefix, filterArea) {
  var ndviP90 = metrics.select(prefix + '_NDVI_P90');
  var ndviMin = metrics.select(prefix + '_NDVI_MIN');
  var ndviAmp = metrics.select(prefix + '_NDVI_AMP');
  var saviMax = metrics.select(prefix + '_SAVI_MAX');
  var ndmiMean = metrics.select(prefix + '_NDMI_MEAN');
  var bsiMed = metrics.select(prefix + '_BSI_MED');

  // Phenology-oriented crop logic:
  // 1) season reaches a green peak
  // 2) season shows noticeable rise/fall amplitude
  // 3) baseline stays low enough to avoid evergreen/natural vegetation
  // 4) soil/moisture support crop signal
  var peakGreen = ndviP90.gt(0.32);
  var seasonalSwing = ndviAmp.gt(0.12);
  var lowBaseline = ndviMin.lt(0.45);
  var soilAdjusted = saviMax.gt(0.18);
  var someMoisture = ndmiMean.gt(-0.08);
  var notBare = bsiMed.lt(0.15);

  var farmland = peakGreen
    .and(seasonalSwing)
    .and(lowBaseline)
    .and(soilAdjusted)
    .and(someMoisture)
    .and(notBare)
    .rename(prefix + '_Farmland');

  var minAreaHa = 0.5;
  var minAreaPixels = (minAreaHa * 10000) / (10 * 10);

  var pixelCount = farmland.selfMask().connectedPixelCount({
    maxSize: 256,
    eightConnected: true
  });

  return farmland.updateMask(pixelCount.gte(minAreaPixels)).clip(filterArea);
}

function areaStats(maskImage, bandName, geometry, name) {
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
    width: '370px',
    padding: '12px',
    backgroundColor: 'white'
  }
});

var legendPanel = ui.Panel({
  style: {
    position: 'bottom-left',
    width: '320px',
    padding: '12px 16px',
    backgroundColor: 'white',
    border: '2px solid #333333',
    borderRadius: '8px'
  }
});

map.add(legendPanel);

controlPanel.add(ui.Label({
  value: 'Somalia Phenology-Based Farmland Detection',
  style: {
    fontWeight: 'bold',
    fontSize: '16px',
    margin: '0 0 8px 0'
  }
}));

controlPanel.add(ui.Label({
  value: 'Pick a SWALIM-informed agricultural zone and adjust Gu and Deyr time windows. The algorithm uses seasonal NDVI peak and amplitude, plus SAVI, NDMI, and BSI.',
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

// Sentinel-2 public archive starts in 2015, but practical L2A coverage is later.
// Start slider from 2019 for cleaner behavior.
var sliderStart = new Date('2019-01-01');
var sliderEnd = new Date();

var defaultGuStart = new Date('2025-05-01');
var defaultGuEnd = new Date('2025-06-30');
var defaultDeyrStart = new Date('2025-11-01');
var defaultDeyrEnd = new Date('2025-12-31');

if (sliderEnd < defaultGuEnd) {
  defaultGuStart = new Date('2024-05-01');
  defaultGuEnd = new Date('2024-06-30');
}
if (sliderEnd < defaultDeyrEnd) {
  defaultDeyrStart = new Date('2024-11-01');
  defaultDeyrEnd = new Date('2024-12-31');
}

controlPanel.add(ui.Label('Gu Date Range', {
  fontWeight: 'bold',
  fontSize: '12px',
  margin: '12px 0 4px 0'
}));

var guDateLabel = ui.Label('Gu Range:', {
  fontSize: '11px',
  color: 'gray',
  margin: '0 0 4px 0'
});
controlPanel.add(guDateLabel);

var guDateSlider = ui.DateSlider({
  start: sliderStart,
  end: sliderEnd,
  value: [defaultGuStart, defaultGuEnd],
  period: 5,
  style: {stretch: 'horizontal'}
});
controlPanel.add(guDateSlider);

controlPanel.add(ui.Label('Deyr Date Range', {
  fontWeight: 'bold',
  fontSize: '12px',
  margin: '12px 0 4px 0'
}));

var deyrDateLabel = ui.Label('Deyr Range:', {
  fontSize: '11px',
  color: 'gray',
  margin: '0 0 4px 0'
});
controlPanel.add(deyrDateLabel);

var deyrDateSlider = ui.DateSlider({
  start: sliderStart,
  end: sliderEnd,
  value: [defaultDeyrStart, defaultDeyrEnd],
  period: 5,
  style: {stretch: 'horizontal'}
});
controlPanel.add(deyrDateSlider);

var runButton = ui.Button({
  label: 'Run Analysis',
  style: {stretch: 'horizontal', margin: '12px 0 0 0'}
});
controlPanel.add(runButton);

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

function drawLegend(selectedZone, guRange, deyrRange) {
  legendPanel.clear();

  legendPanel.add(ui.Label({
    value: 'PHENOLOGY FARMLAND LEGEND',
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
    value: 'Classes',
    style: {fontWeight: 'bold', fontSize: '11px', margin: '10px 0 4px 0'}
  }));

  legendPanel.add(makeLegendRow('#00cc44', 'Gu Farmland', 'Detected by seasonal phenology in Gu window'));
  legendPanel.add(makeLegendRow('#0066ff', 'Deyr Farmland', 'Detected by seasonal phenology in Deyr window'));
  legendPanel.add(makeLegendRow('#ffff00', 'Persistent Farmland', 'Detected in both selected windows'));

  legendPanel.add(ui.Label({
    value: 'Selected Dates',
    style: {fontWeight: 'bold', fontSize: '11px', margin: '10px 0 4px 0'}
  }));

  legendPanel.add(ui.Label(
    'Gu: ' + formatJsDate(guRange.start) + ' to ' + formatJsDate(guRange.end),
    {fontSize: '10px'}
  ));
  legendPanel.add(ui.Label(
    'Deyr: ' + formatJsDate(deyrRange.start) + ' to ' + formatJsDate(deyrRange.end),
    {fontSize: '10px'}
  ));

  legendPanel.add(ui.Label({
    value: 'Phenology Rules',
    style: {fontWeight: 'bold', fontSize: '11px', margin: '10px 0 4px 0'}
  }));

  legendPanel.add(ui.Label('NDVI P90 > 0.32', {fontSize: '10px'}));
  legendPanel.add(ui.Label('NDVI amplitude > 0.12', {fontSize: '10px'}));
  legendPanel.add(ui.Label('NDVI minimum < 0.45', {fontSize: '10px'}));
  legendPanel.add(ui.Label('SAVI max > 0.18', {fontSize: '10px'}));
  legendPanel.add(ui.Label('NDMI mean > -0.08', {fontSize: '10px'}));
  legendPanel.add(ui.Label('BSI median < 0.15', {fontSize: '10px'}));
  legendPanel.add(ui.Label('Minimum patch size = 0.5 ha', {fontSize: '10px'}));

  legendPanel.add(ui.Label({
    value: 'Sentinel-2 | 10m | 5-day slider step',
    style: {fontSize: '9px', color: 'gray', margin: '10px 0 0 0'}
  }));
}

// ============================================================================
// 6. MAIN ANALYSIS
// ============================================================================

function runAnalysis() {
  var selectedZone = zoneSelect.getValue();
  var filterArea = getSelectedZoneGeometry(selectedZone);

  var guRange = getDateRangeFromSlider(guDateSlider);
  var deyrRange = getDateRangeFromSlider(deyrDateSlider);

  updateDateLabel(guDateLabel, 'Gu Range', guRange);
  updateDateLabel(deyrDateLabel, 'Deyr Range', deyrRange);

  statusLabel.setValue('Running analysis for: ' + selectedZone + ' ...');

  map.layers().reset();
  map.centerObject(filterArea, selectedZone === 'All SWALIM Agricultural Areas' ? 6 : 8);

  map.addLayer(
    ee.Image().paint(somaliaAOI, 1, 2),
    {palette: ['999999']},
    'Somalia AOI',
    false
  );

  map.addLayer(
    ee.Image().paint(filterArea, 1, 2),
    {palette: ['yellow']},
    'Selected Agricultural Filter Area',
    true
  );

  var guCollection = getCollection(
  filterArea,
  formatJsDate(guRange.start),
  formatJsDate(guRange.end)
);

var deyrCollection = getCollection(
  filterArea,
  formatJsDate(deyrRange.start),
  formatJsDate(deyrRange.end)
);

  var guComposite = getMedianComposite(guCollection, filterArea);
  var deyrComposite = getMedianComposite(deyrCollection, filterArea);

  var guMetrics = buildPhenologyMetrics(guCollection, 'Gu', filterArea);
  var deyrMetrics = buildPhenologyMetrics(deyrCollection, 'Deyr', filterArea);

  var guFarmland = classifyPhenologyFarmland(guMetrics, 'Gu', filterArea);
  var deyrFarmland = classifyPhenologyFarmland(deyrMetrics, 'Deyr', filterArea);

  var persistentFarmland = guFarmland.and(deyrFarmland).rename('Persistent_Farmland');
  var guOnly = guFarmland.and(deyrFarmland.not()).rename('Gu_Only');
  var deyrOnly = deyrFarmland.and(guFarmland.not()).rename('Deyr_Only');
  var anyFarmland = guFarmland.or(deyrFarmland).rename('Any_Farmland');

  var ndviPeakDiff = guMetrics.select('Gu_NDVI_P90')
    .subtract(deyrMetrics.select('Deyr_NDVI_P90'))
    .rename('NDVI_Peak_Diff')
    .clip(filterArea);

  var guAmp = guMetrics.select('Gu_NDVI_AMP');
  var deyrAmp = deyrMetrics.select('Deyr_NDVI_AMP');

  // Base imagery
  map.addLayer(
    guComposite.select(['B4', 'B8', 'B11']),
    {min: 0.02, max: 0.35},
    'Gu Agriculture Composite',
    false
  );

  map.addLayer(
    deyrComposite.select(['B4', 'B8', 'B11']),
    {min: 0.02, max: 0.35},
    'Deyr Agriculture Composite',
    false
  );

  // Phenology layers
  map.addLayer(
    guMetrics.select('Gu_NDVI_P90'),
    {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
    'Gu NDVI Peak (P90)',
    false
  );

  map.addLayer(
    deyrMetrics.select('Deyr_NDVI_P90'),
    {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green', 'darkgreen']},
    'Deyr NDVI Peak (P90)',
    false
  );

  map.addLayer(
    guAmp,
    {min: 0, max: 0.4, palette: ['white', 'yellow', 'orange', 'red']},
    'Gu NDVI Amplitude',
    false
  );

  map.addLayer(
    deyrAmp,
    {min: 0, max: 0.4, palette: ['white', 'yellow', 'orange', 'red']},
    'Deyr NDVI Amplitude',
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
    ndviPeakDiff,
    {min: -0.35, max: 0.35, palette: ['blue', 'white', 'green']},
    'NDVI Peak Difference (Gu - Deyr)',
    false
  );

  drawLegend(selectedZone, guRange, deyrRange);

  print('==================================================');
  print('PHENOLOGY-BASED IMAGE AVAILABILITY');
  print('==================================================');
  print('Selected zone:', selectedZone);
  print('Gu start:', formatJsDate(guRange.start));
  print('Gu end:', formatJsDate(guRange.end));
  print('Deyr start:', formatJsDate(deyrRange.start));
  print('Deyr end:', formatJsDate(deyrRange.end));
  print('Gu scenes:', guCollection.size());
  print('Deyr scenes:', deyrCollection.size());

  print('==================================================');
  print('PHENOLOGY-BASED FARMLAND STATISTICS');
  print('==================================================');
  print('Somalia AOI: Entire Somalia');
  print('Selected agricultural filter area:', selectedZone);

  areaStats(guFarmland, 'Gu_Farmland', filterArea, 'Gu farmland');
  areaStats(deyrFarmland, 'Deyr_Farmland', filterArea, 'Deyr farmland');
  areaStats(persistentFarmland, 'Persistent_Farmland', filterArea, 'Persistent farmland');
  areaStats(anyFarmland, 'Any_Farmland', filterArea, 'Any seasonal farmland');

  statusLabel.setValue('Analysis complete for: ' + selectedZone);
}

// ============================================================================
// 7. EVENTS
// ============================================================================

runButton.onClick(function() {
  runAnalysis();
});

zoneSelect.onChange(function() {
  updateDateLabel(guDateLabel, 'Gu Range', getDateRangeFromSlider(guDateSlider));
  updateDateLabel(deyrDateLabel, 'Deyr Range', getDateRangeFromSlider(deyrDateSlider));
});

guDateSlider.onChange(function() {
  updateDateLabel(guDateLabel, 'Gu Range', getDateRangeFromSlider(guDateSlider));
});

deyrDateSlider.onChange(function() {
  updateDateLabel(deyrDateLabel, 'Deyr Range', getDateRangeFromSlider(deyrDateSlider));
});

// ============================================================================
// 8. INITIAL LABELS + RUN
// ============================================================================

updateDateLabel(guDateLabel, 'Gu Range', getDateRangeFromSlider(guDateSlider));
updateDateLabel(deyrDateLabel, 'Deyr Range', getDateRangeFromSlider(deyrDateSlider));
runAnalysis();