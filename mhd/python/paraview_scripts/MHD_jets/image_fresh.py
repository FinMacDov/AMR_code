# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
jetB50V40_0 = FindSource('jetB50V40_0*')

# create a new 'Clip'
clip1 = Clip(Input=jetB50V40_0)
clip1.ClipType = 'Plane'
clip1.Scalars = ['CELLS', 'Alfv']
clip1.Value = 59146291.734375

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.0, 4499999744.0, 0.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [0.0, 1000000000.0, 0.0]
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [0.0, 1000000000.0, 0.0]
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [969, 390]

# show data in view
clip1Display = Show(clip1, renderView1)

# get color transfer function/color map for 'rho'
rhoLUT = GetColorTransferFunction('rho')

# get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['CELLS', 'rho']
clip1Display.LookupTable = rhoLUT
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 100000000.0
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 5000000.0
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.SelectionCellLabelFontFile = ''
clip1Display.SelectionPointLabelFontFile = ''
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = rhoPWF
clip1Display.ScalarOpacityUnitDistance = 25834862.904569723

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip1Display.DataAxesGrid.XTitleFontFile = ''
clip1Display.DataAxesGrid.YTitleFontFile = ''
clip1Display.DataAxesGrid.ZTitleFontFile = ''
clip1Display.DataAxesGrid.XLabelFontFile = ''
clip1Display.DataAxesGrid.YLabelFontFile = ''
clip1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(jetB50V40_0, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# hide data in view
Hide(clip1, renderView1)

# set active source
SetActiveSource(clip1)

# show data in view
clip1Display = Show(clip1, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
rhoLUT.ApplyPreset('Black-Body Radiation', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
rhoLUT.ApplyPreset('X Ray', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
rhoLUT.ApplyPreset('Black-Body Radiation', True)

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [950, 1025]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.StereoType = 0
renderView2.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView2.AxesGrid.Visibility = 1
renderView2.AxesGrid.YTitle = '     Y Axis'
renderView2.AxesGrid.XTitleFontFile = ''
renderView2.AxesGrid.XTitleBold = 1
renderView2.AxesGrid.XTitleFontSize = 15
renderView2.AxesGrid.YTitleFontFile = ''
renderView2.AxesGrid.YTitleBold = 1
renderView2.AxesGrid.YTitleFontSize = 15
renderView2.AxesGrid.ZTitleFontFile = ''
renderView2.AxesGrid.ZTitleBold = 1
renderView2.AxesGrid.ZTitleFontSize = 15
renderView2.AxesGrid.XLabelFontFile = ''
renderView2.AxesGrid.XLabelBold = 1
renderView2.AxesGrid.XLabelFontSize = 15
renderView2.AxesGrid.YLabelFontFile = ''
renderView2.AxesGrid.YLabelBold = 1
renderView2.AxesGrid.YLabelFontSize = 15
renderView2.AxesGrid.ZLabelFontFile = ''
renderView2.AxesGrid.ZLabelBold = 1
renderView2.AxesGrid.ZLabelFontSize = 15

# place view in the layout
layout1.AssignView(2, renderView2)

# split cell
layout1.SplitHorizontal(2, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.ViewSize = [470, 1025]
renderView3.AxesGrid = 'GridAxes3DActor'
renderView3.StereoType = 0
renderView3.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView3.AxesGrid.Visibility = 1
renderView3.AxesGrid.YTitle = '     Y Axis'
renderView3.AxesGrid.XTitleFontFile = ''
renderView3.AxesGrid.XTitleBold = 1
renderView3.AxesGrid.XTitleFontSize = 15
renderView3.AxesGrid.YTitleFontFile = ''
renderView3.AxesGrid.YTitleBold = 1
renderView3.AxesGrid.YTitleFontSize = 15
renderView3.AxesGrid.ZTitleFontFile = ''
renderView3.AxesGrid.ZTitleBold = 1
renderView3.AxesGrid.ZTitleFontSize = 15
renderView3.AxesGrid.XLabelFontFile = ''
renderView3.AxesGrid.XLabelBold = 1
renderView3.AxesGrid.XLabelFontSize = 15
renderView3.AxesGrid.YLabelFontFile = ''
renderView3.AxesGrid.YLabelBold = 1
renderView3.AxesGrid.YLabelFontSize = 15
renderView3.AxesGrid.ZLabelFontFile = ''
renderView3.AxesGrid.ZLabelBold = 1
renderView3.AxesGrid.ZLabelFontSize = 15

# place view in the layout
layout1.AssignView(6, renderView3)

# resize frame
layout1.SetSplitFraction(0, 0.405645582854)

# resize frame
layout1.SetSplitFraction(0, 0.36696288552)

# resize frame
layout1.SetSplitFraction(2, 0.53642384106)

# set active view
SetActiveView(renderView1)

# hide data in view
Hide(clip1, renderView1)

# show data in view
clip1Display = Show(clip1, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# get color legend/bar for rhoLUT in view renderView1
rhoLUTColorBar = GetScalarBar(rhoLUT, renderView1)

# change scalar bar placement
rhoLUTColorBar.Orientation = 'Horizontal'
rhoLUTColorBar.Position = [0.22659585938768312, 0.08885017421602807]
rhoLUTColorBar.ScalarBarLength = 0.6442857142857145

# set active view
SetActiveView(renderView2)

#change interaction mode for render view
renderView2.InteractionMode = '2D'

# set active view
SetActiveView(renderView3)

#change interaction mode for render view
renderView3.InteractionMode = '2D'

# set active view
SetActiveView(renderView1)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# set active view
SetActiveView(renderView2)

# Hide orientation axes
renderView2.OrientationAxesVisibility = 0

# set active view
SetActiveView(renderView3)

# Hide orientation axes
renderView3.OrientationAxesVisibility = 0

# set active view
SetActiveView(renderView1)

# change scalar bar placement
rhoLUTColorBar.ScalarBarLength = 0.7362397372742202

# change scalar bar placement
rhoLUTColorBar.Position = [0.1317682731807866, 0.10055749128919882]
rhoLUTColorBar.ScalarBarLength = 0.7362397372742202

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.HorizontalTitle = 1

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.AutoOrient = 0
rhoLUTColorBar.Orientation = 'Vertical'

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.Orientation = 'Horizontal'

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.TitleJustification = 'Right'

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.TitleJustification = 'Centered'

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.ComponentTitle = 'test'

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.ComponentTitle = ''

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.HorizontalTitle = 0

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.TitleFontSize = 15

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.TitleFontSize = 19

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.TitleFontSize = 20

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.LabelFontSize = 20

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.LabelFontSize = 15

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.ScalarBarThickness = 20

# Properties modified on rhoLUTColorBar
rhoLUTColorBar.ScalarBarThickness = 15

# Properties modified on clip1Display.DataAxesGrid
clip1Display.DataAxesGrid.XTitleFontSize = 15
clip1Display.DataAxesGrid.YTitleFontSize = 15
clip1Display.DataAxesGrid.ZTitleFontSize = 15

# Properties modified on clip1Display.DataAxesGrid
clip1Display.DataAxesGrid.XLabelFontSize = 15
clip1Display.DataAxesGrid.YLabelFontSize = 15
clip1Display.DataAxesGrid.ZLabelFontSize = 15

# Properties modified on clip1Display.DataAxesGrid
clip1Display.DataAxesGrid.ShowGrid = 1

# Properties modified on clip1Display.DataAxesGrid
clip1Display.DataAxesGrid.ShowGrid = 0

# Properties modified on clip1Display.DataAxesGrid
clip1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = 'Y Axis'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = '    Y Axis'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = '      Y Axis'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = '         Y Axis'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = 'Y Axis'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = '   Y Axis'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YTitle = '    Y Axis'

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XLabelFontSize = 21

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XLabelFontSize = 15

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YLabelFontSize = 18

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YLabelFontSize = 20

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.YLabelFontSize = 18

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XLabelFontSize = 18

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = ''
renderView1.AxesGrid.YTitle = ''
renderView1.AxesGrid.XTitleFontSize = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XLabelFontSize = 20
renderView1.AxesGrid.YLabelFontSize = 20

# Show orientation axes
renderView1.OrientationAxesVisibility = 1

# change scalar bar placement
rhoLUTColorBar.Position = [0.12458436513480964, 0.10641114982578412]
rhoLUTColorBar.ScalarBarLength = 0.7362397372742199

# change scalar bar placement
rhoLUTColorBar.Position = [0.14757287088193616, 0.8283623693379791]

# set active view
SetActiveView(renderView2)

# show data in view
clip1Display_1 = Show(clip1, renderView2)

# trace defaults for the display properties.
clip1Display_1.Representation = 'Surface'
clip1Display_1.ColorArrayName = [None, '']
clip1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display_1.SelectOrientationVectors = 'None'
clip1Display_1.ScaleFactor = 100000000.0
clip1Display_1.SelectScaleArray = 'None'
clip1Display_1.GlyphType = 'Arrow'
clip1Display_1.GlyphTableIndexArray = 'None'
clip1Display_1.GaussianRadius = 5000000.0
clip1Display_1.SetScaleArray = [None, '']
clip1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display_1.OpacityArray = [None, '']
clip1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display_1.DataAxesGrid = 'GridAxesRepresentation'
clip1Display_1.SelectionCellLabelFontFile = ''
clip1Display_1.SelectionPointLabelFontFile = ''
clip1Display_1.PolarAxes = 'PolarAxesRepresentation'
clip1Display_1.ScalarOpacityUnitDistance = 25834862.904569723

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip1Display_1.DataAxesGrid.XTitleFontFile = ''
clip1Display_1.DataAxesGrid.YTitleFontFile = ''
clip1Display_1.DataAxesGrid.ZTitleFontFile = ''
clip1Display_1.DataAxesGrid.XLabelFontFile = ''
clip1Display_1.DataAxesGrid.YLabelFontFile = ''
clip1Display_1.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip1Display_1.PolarAxes.PolarAxisTitleFontFile = ''
clip1Display_1.PolarAxes.PolarAxisLabelFontFile = ''
clip1Display_1.PolarAxes.LastRadialAxisTextFontFile = ''
clip1Display_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView2.ResetCamera()

# set scalar coloring
ColorBy(clip1Display_1, ('CELLS', 'Te'))

# rescale color and/or opacity maps used to include current data range
clip1Display_1.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display_1.SetScalarBarVisibility(renderView2, True)

# get color transfer function/color map for 'Te'
teLUT = GetColorTransferFunction('Te')

# get color legend/bar for teLUT in view renderView2
teLUTColorBar = GetScalarBar(teLUT, renderView2)

# change scalar bar placement
teLUTColorBar.Position = [0.12398011118701718, 0.8186062717770036]
teLUTColorBar.ScalarBarLength = 0.73623973727422

# change scalar bar placement
teLUTColorBar.ScalarBarLength = 0.7409126344704784

# change scalar bar placement
teLUTColorBar.Position = [0.1457869647695717, 0.8147038327526134]

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_1.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_1.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_1.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_1.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_1.RescaleTransferFunctionToDataRange(False, True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Cool to Warm', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Cool to Warm (Extended)', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Cool to Warm', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Cold and Hot', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Rainbow Desaturated', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Warm to Cool (Extended)', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Warm to Cool', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
teLUT.ApplyPreset('Cool to Warm', True)

# set active view
SetActiveView(renderView3)

# show data in view
clip1Display_2 = Show(clip1, renderView3)

# trace defaults for the display properties.
clip1Display_2.Representation = 'Surface'
clip1Display_2.ColorArrayName = [None, '']
clip1Display_2.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display_2.SelectOrientationVectors = 'None'
clip1Display_2.ScaleFactor = 100000000.0
clip1Display_2.SelectScaleArray = 'None'
clip1Display_2.GlyphType = 'Arrow'
clip1Display_2.GlyphTableIndexArray = 'None'
clip1Display_2.GaussianRadius = 5000000.0
clip1Display_2.SetScaleArray = [None, '']
clip1Display_2.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display_2.OpacityArray = [None, '']
clip1Display_2.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display_2.DataAxesGrid = 'GridAxesRepresentation'
clip1Display_2.SelectionCellLabelFontFile = ''
clip1Display_2.SelectionPointLabelFontFile = ''
clip1Display_2.PolarAxes = 'PolarAxesRepresentation'
clip1Display_2.ScalarOpacityUnitDistance = 25834862.904569723

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip1Display_2.DataAxesGrid.XTitleFontFile = ''
clip1Display_2.DataAxesGrid.YTitleFontFile = ''
clip1Display_2.DataAxesGrid.ZTitleFontFile = ''
clip1Display_2.DataAxesGrid.XLabelFontFile = ''
clip1Display_2.DataAxesGrid.YLabelFontFile = ''
clip1Display_2.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip1Display_2.PolarAxes.PolarAxisTitleFontFile = ''
clip1Display_2.PolarAxes.PolarAxisLabelFontFile = ''
clip1Display_2.PolarAxes.LastRadialAxisTextFontFile = ''
clip1Display_2.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView3.ResetCamera()

# set scalar coloring
ColorBy(clip1Display_2, ('CELLS', 'trp1'))

# rescale color and/or opacity maps used to include current data range
clip1Display_2.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display_2.SetScalarBarVisibility(renderView3, True)

# get color transfer function/color map for 'trp1'
trp1LUT = GetColorTransferFunction('trp1')

# set active view
SetActiveView(renderView1)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# Properties modified on clip1Display.DataAxesGrid
clip1Display.DataAxesGrid.XTitle = ''
clip1Display.DataAxesGrid.YTitle = ''
clip1Display.DataAxesGrid.ZTitle = ''
clip1Display.DataAxesGrid.XTitleFontSize = 18
clip1Display.DataAxesGrid.YTitleFontSize = 18
clip1Display.DataAxesGrid.ZTitleFontSize = 18
clip1Display.DataAxesGrid.XLabelFontSize = 18
clip1Display.DataAxesGrid.YLabelFontSize = 18
clip1Display.DataAxesGrid.ZLabelFontSize = 18

# set active view
SetActiveView(renderView2)

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.XTitle = ''
renderView2.AxesGrid.YTitle = ''
renderView2.AxesGrid.ZTitle = ''
renderView2.AxesGrid.XTitleFontSize = 18
renderView2.AxesGrid.YTitleFontSize = 18
renderView2.AxesGrid.ZTitleFontSize = 18
renderView2.AxesGrid.XLabelFontSize = 18
renderView2.AxesGrid.YLabelFontSize = 18
renderView2.AxesGrid.ZLabelFontSize = 18

# hide data in view
Hide(clip1, renderView2)

# show data in view
clip1Display_1 = Show(clip1, renderView2)

# show color bar/color legend
clip1Display_1.SetScalarBarVisibility(renderView2, True)

# reset view to fit data
renderView2.ResetCamera()

# hide data in view
Hide(clip1, renderView2)

# show data in view
clip1Display_1 = Show(clip1, renderView2)

# show color bar/color legend
clip1Display_1.SetScalarBarVisibility(renderView2, True)

# reset view to fit data
renderView2.ResetCamera()

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.Visibility = 1

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.Visibility = 1

# Properties modified on renderView2.AxesGrid
renderView2.AxesGrid.XTitle = ''
renderView2.AxesGrid.YTitle = ''
renderView2.AxesGrid.ZTitle = ''
renderView2.AxesGrid.XLabelBold = 1
renderView2.AxesGrid.XLabelFontSize = 18
renderView2.AxesGrid.YLabelBold = 1
renderView2.AxesGrid.YLabelFontSize = 18

# set active view
SetActiveView(renderView3)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
trp1LUT.ApplyPreset('Black-Body Radiation', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
trp1LUT.ApplyPreset('X Ray', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
trp1LUT.ApplyPreset('Grayscale', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
trp1LUT.ApplyPreset('X Ray', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
trp1LUT.ApplyPreset('Black, Orange and White', True)

# Properties modified on renderView3.AxesGrid
renderView3.AxesGrid.XTitle = ''
renderView3.AxesGrid.YTitle = ''
renderView3.AxesGrid.ZTitle = ''
renderView3.AxesGrid.XLabelFontSize = 18
renderView3.AxesGrid.YLabelFontSize = 18

# get color legend/bar for trp1LUT in view renderView3
trp1LUTColorBar = GetScalarBar(trp1LUT, renderView3)

# change scalar bar placement
trp1LUTColorBar.Position = [0.07220148617717653, 0.8029965156794427]
trp1LUTColorBar.ScalarBarLength = 0.7362397372742202

# change scalar bar placement
trp1LUTColorBar.ScalarBarLength = 0.8986946109204288

# change scalar bar placement
trp1LUTColorBar.Position = [0.05956610711580462, 0.8078745644599306]

# change scalar bar placement
trp1LUTColorBar.ScalarBarLength = 0.8553733112814401

# change scalar bar placement
trp1LUTColorBar.Position = [0.05956610711580462, 0.8147038327526135]

# set active view
SetActiveView(renderView1)

# change scalar bar placement
rhoLUTColorBar.ScalarBarLength = 0.6873891625615762

# change scalar bar placement
rhoLUTColorBar.Position = [0.16050390536469472, 0.826411149825784]

# change scalar bar placement
rhoLUTColorBar.ScalarBarLength = 0.7146880131362888

# change scalar bar placement
rhoLUTColorBar.Position = [0.14038896283595903, 0.826411149825784]

# set active view
SetActiveView(renderView2)

# change scalar bar placement
teLUTColorBar.Position = [0.13644117037704834, 0.8205574912891987]
teLUTColorBar.ScalarBarLength = 0.7409126344704786

# set active view
SetActiveView(renderView3)

# change scalar bar placement
trp1LUTColorBar.Position = [0.09386213599667106, 0.8166550522648086]
trp1LUTColorBar.ScalarBarLength = 0.8553733112814405

# set active view
SetActiveView(renderView2)

# set active view
SetActiveView(renderView3)

# set active view
SetActiveView(renderView2)

# Show orientation axes
renderView2.OrientationAxesVisibility = 1

# set active view
SetActiveView(renderView3)

# Show orientation axes
renderView3.OrientationAxesVisibility = 1

#### saving camera placements for all active views

# current camera placement for renderView2
renderView2.InteractionMode = '2D'
renderView2.CameraPosition = [0.0, 488031543.4814339, 3868654631.9850993]
renderView2.CameraFocalPoint = [0.0, 488031543.4814339, 0.0]
renderView2.CameraParallelScale = 965511381.4789301

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 500000000.0, 3577058526.585485]
renderView1.CameraFocalPoint = [0.0, 500000000.0, 0.0]
renderView1.CameraParallelScale = 942988842.4380271

# current camera placement for renderView3
renderView3.InteractionMode = '2D'
renderView3.CameraPosition = [0.0, 496179178.1857369, 4467462039.9058485]
renderView3.CameraFocalPoint = [0.0, 496179178.1857369, 0.0]
renderView3.CameraParallelScale = 979085589.9049242

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
