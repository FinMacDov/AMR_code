# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
clip1 = FindSource('Clip1')

# set active source
SetActiveSource(clip1)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [696, 1025]

# get display properties
clip1Display = GetDisplayProperties(clip1, view=renderView1)

# get color transfer function/color map for 'rho'
rhoLUT = GetColorTransferFunction('rho')

# set scalar coloring
ColorBy(clip1Display, ('CELLS', 'v1'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(rhoLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'v1'
v1LUT = GetColorTransferFunction('v1')

# find view
renderView2 = FindViewOrCreate('RenderView2', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView2.ViewSize = [642, 1025]

# set active view
SetActiveView(renderView2)

# get display properties
clip1Display_1 = GetDisplayProperties(clip1, view=renderView2)

# get color transfer function/color map for 'Te'
teLUT = GetColorTransferFunction('Te')

# set active view
SetActiveView(renderView1)

# get color legend/bar for v1LUT in view renderView1
v1LUTColorBar = GetScalarBar(v1LUT, renderView1)

# change scalar bar placement
v1LUTColorBar.Orientation = 'Horizontal'
v1LUTColorBar.Position = [0.20901656314699774, 0.7863827283762905]
v1LUTColorBar.ScalarBarLength = 0.6442857142857152

# change scalar bar placement
v1LUTColorBar.Position = [0.19177518383665287, 0.8195534600836075]

# change scalar bar placement
v1LUTColorBar.ScalarBarLength = 0.6916995073891627

# change scalar bar placement
v1LUTColorBar.Position = [0.1357407010780322, 0.81857785032751]

# Properties modified on v1LUTColorBar
v1LUTColorBar.HorizontalTitle = 1

# Properties modified on v1LUTColorBar
v1LUTColorBar.TitleFontSize = 18

# Properties modified on v1LUTColorBar
v1LUTColorBar.LabelFontSize = 18

# change scalar bar placement
v1LUTColorBar.Position = [0.15729242521596326, 0.81857785032751]

# get animation scene
animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

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

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# set active view
SetActiveView(renderView2)

# set scalar coloring
ColorBy(clip1Display_1, ('CELLS', 'v2'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(teLUT, renderView2)

# rescale color and/or opacity maps used to include current data range
clip1Display_1.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display_1.SetScalarBarVisibility(renderView2, True)

# get color transfer function/color map for 'v2'
v2LUT = GetColorTransferFunction('v2')

# get color legend/bar for v2LUT in view renderView2
v2LUTColorBar = GetScalarBar(v2LUT, renderView2)

# change scalar bar placement
v2LUTColorBar.Position = [0.1489022295670795, 0.8234843205574914]
v2LUTColorBar.ScalarBarLength = 0.73623973727422

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
v2LUT.ApplyPreset('Cool to Warm', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
v2LUT.ApplyPreset('Cool to Warm (Extended)', True)

# set active view
SetActiveView(renderView1)

# Properties modified on v1LUTColorBar
v1LUTColorBar.ScalarBarThickness = 15

# find view
renderView3 = FindViewOrCreate('RenderView3', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView3.ViewSize = [554, 1025]

# set active view
SetActiveView(renderView3)

# get display properties
clip1Display_2 = GetDisplayProperties(clip1, view=renderView3)

# get color transfer function/color map for 'trp1'
trp1LUT = GetColorTransferFunction('trp1')

# set scalar coloring
ColorBy(clip1Display_2, ('CELLS', 'schrho'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(trp1LUT, renderView3)

# rescale color and/or opacity maps used to include current data range
clip1Display_2.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display_2.SetScalarBarVisibility(renderView3, True)

# get color transfer function/color map for 'schrho'
schrhoLUT = GetColorTransferFunction('schrho')

# get color legend/bar for schrhoLUT in view renderView3
schrhoLUTColorBar = GetScalarBar(schrhoLUT, renderView3)

# change scalar bar placement
schrhoLUTColorBar.Position = [0.11010762336129204, 0.8390940766550524]
schrhoLUTColorBar.ScalarBarLength = 0.73623973727422

# change scalar bar placement
schrhoLUTColorBar.ScalarBarLength = 0.8481530946749436

# change scalar bar placement
schrhoLUTColorBar.Position = [0.0866419193901729, 0.824459930313589]

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_2.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_2.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display_2.RescaleTransferFunctionToDataRange(False, True)

# change scalar bar placement
schrhoLUTColorBar.Position = [0.09386213599667109, 0.8273867595818817]

# set active view
SetActiveView(renderView1)

# change scalar bar placement
v1LUTColorBar.Position = [0.15729242521596326, 0.8195534600836075]

# Properties modified on v1LUTColorBar
v1LUTColorBar.AddRangeAnnotations = 1

# Properties modified on v1LUTColorBar
v1LUTColorBar.AddRangeAnnotations = 0

# Properties modified on v1LUTColorBar
v1LUTColorBar.AutomaticAnnotations = 1

# Properties modified on v1LUTColorBar
v1LUTColorBar.AutomaticAnnotations = 0

# Properties modified on v1LUTColorBar
v1LUTColorBar.AddRangeLabels = 0

# Properties modified on v1LUTColorBar
v1LUTColorBar.AddRangeLabels = 1

# Properties modified on v1LUTColorBar
v1LUTColorBar.AutomaticLabelFormat = 0

# change scalar bar placement
v1LUTColorBar.ScalarBarLength = 0.7836535303776733

# change scalar bar placement
v1LUTColorBar.Position = [0.07252231027343453, 0.8107729722787295]

# change scalar bar placement
v1LUTColorBar.ScalarBarLength = 0.7735960591133026

# change scalar bar placement
v1LUTColorBar.Position = [0.13717748268722757, 0.8156510210592173]

# Properties modified on v1LUTColorBar
v1LUTColorBar.AddRangeLabels = 0

# Properties modified on v1LUTColorBar
v1LUTColorBar.AutomaticLabelFormat = 1
v1LUTColorBar.AddRangeLabels = 1

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
renderView3.CameraPosition = [-1910410.9071315592, 494268767.27860534, 4467462039.9058485]
renderView3.CameraFocalPoint = [-1910410.9071315592, 494268767.27860534, 0.0]
renderView3.CameraParallelScale = 979085589.9049242

# get layout
layout1 = GetLayout()

# save animation
SaveAnimation('/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet2/fresh/movie/v1_v2_sch.png', layout1, SaveAllViews=1,
    ImageResolution=[1892, 1025],
    FrameRate=2,
    FrameWindow=[0, 200])

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
renderView3.CameraPosition = [-1910410.9071315592, 494268767.27860534, 4467462039.9058485]
renderView3.CameraFocalPoint = [-1910410.9071315592, 494268767.27860534, 0.0]
renderView3.CameraParallelScale = 979085589.9049242

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
