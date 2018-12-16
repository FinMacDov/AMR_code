# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find view
renderView2 = FindViewOrCreate('RenderView2', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView2.ViewSize = [643, 1025]

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [694, 1025]

# get active view
renderView3 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView3.ViewSize = [555, 1025]

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

# get layout
layout1 = GetLayout()

# save screenshot
SaveScreenshot('/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet2/fresh/Te_rho_tr_eq_0.png', layout1, SaveAllViews=1,
    ImageResolution=[1892, 1025])

# get active source.
clip1 = GetActiveSource()

# set active source
SetActiveSource(clip1)

# set active view
SetActiveView(renderView1)

# get display properties
clip1Display = GetDisplayProperties(clip1, view=renderView1)

# get color transfer function/color map for 'rho'
rhoLUT = GetColorTransferFunction('rho')

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=clip1)
annotateTimeFilter1.Format = 'Time: %g'
annotateTimeFilter1.Scale = 2.14683

# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

# trace defaults for the display properties.
annotateTimeFilter1Display.FontFile = ''

# find source
jetB50V40_0 = FindSource('jetB50V40_0*')

# update the view to ensure updated data information
renderView1.Update()

# Rescale transfer function
rhoLUT.RescaleTransferFunction(1.43337795246e-14, 1.87055992917e-09)

# get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho')

# Rescale transfer function
rhoPWF.RescaleTransferFunction(1.43337795246e-14, 1.87055992917e-09)

# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.WindowLocation = 'LowerLeftCorner'

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

# save animation
SaveAnimation('/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet2/fresh/movie/rho_Te_tr.png', layout1, SaveAllViews=1,
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
renderView3.CameraPosition = [0.0, 496179178.1857369, 4467462039.9058485]
renderView3.CameraFocalPoint = [0.0, 496179178.1857369, 0.0]
renderView3.CameraParallelScale = 979085589.9049242

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
