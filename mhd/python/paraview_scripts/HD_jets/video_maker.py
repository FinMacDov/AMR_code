# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
import glob
import os
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# colour_wheel = ['Black-Body Radiation', 'Cool to Warm', 'Cool to Warm (Extended)']
colour_wheel = ['Black, Orange and White', 'GREEN-WHITE_LINEAR', 'Grayscale']
col_black = [0.0, 0.0, 0.0]
col_white = [1.0, 1.0, 1.0]

def render_view_func():
    materialLibrary1 = GetMaterialLibrary()
    renderView = CreateView('RenderView')
    renderView.ViewSize = [1356, 713] # careful: set to a particular value (fixes 3rd panel)
    renderView.AxesGrid = 'GridAxes3DActor'
    renderView.StereoType = 0
    renderView.Background = [0.32, 0.34, 0.43]
    renderView.OSPRayMaterialLibrary = materialLibrary1
    
    # init the 'GridAxes3DActor' selected for 'AxesGrid'
    renderView.AxesGrid.Visibility = 1
    renderView.AxesGrid.YTitle = '     Y Axis'
    renderView.AxesGrid.XTitleFontFile = ''
    renderView.AxesGrid.XTitleBold = 1
    renderView.AxesGrid.XTitleFontSize = 15
    renderView.AxesGrid.YTitleFontFile = ''
    renderView.AxesGrid.YTitleBold = 1
    renderView.AxesGrid.YTitleFontSize = 15
    renderView.AxesGrid.ZTitleFontFile = ''
    renderView.AxesGrid.ZTitleBold = 1
    renderView.AxesGrid.ZTitleFontSize = 15
    renderView.AxesGrid.XLabelFontFile = ''
    renderView.AxesGrid.XLabelBold = 1
    renderView.AxesGrid.XLabelFontSize = 15
    renderView.AxesGrid.YLabelFontFile = ''
    renderView.AxesGrid.YLabelBold = 1
    renderView.AxesGrid.YLabelFontSize = 15
    renderView.AxesGrid.ZLabelFontFile = ''
    renderView.AxesGrid.ZLabelBold = 1
    renderView.AxesGrid.ZLabelFontSize = 15
    
    renderView.OrientationAxesVisibility = 0
    renderView.InteractionMode = '2D'
    
    renderView.CameraPosition = [62068695.867912374, 1006066501.1989748, 6607406493.1831665]
#    renderView.CameraPosition = [0.0, 1000000000.0, 4700000000.0]
    renderView.CameraFocalPoint = [0.0, 1000000000.0, 0.0]
    renderView.CameraParallelScale = 1405561883.6231568
    
    renderView.AxesGrid.XTitle = ''
    renderView.AxesGrid.YTitle = ''

    renderView.AxesGrid.YLabelColor = col_white
    renderView.AxesGrid.XLabelColor = col_white
    renderView.AxesGrid.GridColor = col_white
    # Properties modified on renderView.AxesGrid
#    XRANGE_MIN, XRANGE_MAX = -4e8, 4e8 
#    YRANGE_MIN, YRANGE_MAX  = 0, 2e9
#    renderView.AxesGrid.XAxisUseCustomLabels = 1
#    renderView.AxesGrid.XAxisLabels = np.linspace(XRANGE_MIN,XRANGE_MAX, 5)
#
#    renderView.AxesGrid.YAxisUseCustomLabels = 1
#    renderView.AxesGrid.YAxisLabels = np.linspace(YRANGE_MIN,YRANGE_MAX, 11)
    return renderView


def display_data(data):
    # displays the data in paraview
    dataDisplay = Show(data, render_view[0])
    # trace defaults for the display properties.
    dataDisplay.Representation = 'Surface'
    dataDisplay.ColorArrayName = [None, '']
    dataDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    dataDisplay.SelectOrientationVectors = 'None'
    dataDisplay.ScaleFactor = 899999948.8000001
    dataDisplay.SelectScaleArray = 'None'
    dataDisplay.GlyphType = 'Arrow'
    dataDisplay.GlyphTableIndexArray = 'None'
    dataDisplay.GaussianRadius = 44999997.44
    dataDisplay.SetScaleArray = [None, '']
    dataDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    dataDisplay.OpacityArray = [None, '']
    dataDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    dataDisplay.DataAxesGrid = 'GridAxesRepresentation'
    dataDisplay.SelectionCellLabelFontFile = ''
    dataDisplay.SelectionPointLabelFontFile = ''
    dataDisplay.PolarAxes = 'PolarAxesRepresentation'
    dataDisplay.ScalarOpacityUnitDistance = 167845592.2061122
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    dataDisplay.DataAxesGrid.XTitleFontFile = ''
    dataDisplay.DataAxesGrid.YTitleFontFile = ''
    dataDisplay.DataAxesGrid.ZTitleFontFile = ''
    dataDisplay.DataAxesGrid.XLabelFontFile = ''
    dataDisplay.DataAxesGrid.YLabelFontFile = ''
    dataDisplay.DataAxesGrid.ZLabelFontFile = ''
    return dataDisplay

def clip_display(clip,renderView):
    ## show data in view
    clipDisplay = Show(clip, renderView)
    # trace defaults for the display properties.
    clipDisplay.Representation = 'Surface'
    clipDisplay.ColorArrayName = [None, '']
    clipDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    clipDisplay.SelectOrientationVectors = 'None'
    clipDisplay.ScaleFactor = 200000000.0
    clipDisplay.SelectScaleArray = 'None'
    clipDisplay.GlyphType = 'Arrow'
    clipDisplay.GlyphTableIndexArray = 'None'
    clipDisplay.GaussianRadius = 10000000.0
    clipDisplay.SetScaleArray = [None, '']
    clipDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    clipDisplay.OpacityArray = [None, '']
    clipDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    clipDisplay.DataAxesGrid = 'GridAxesRepresentation'
    clipDisplay.SelectionCellLabelFontFile = ''
    clipDisplay.SelectionPointLabelFontFile = ''
    clipDisplay.PolarAxes = 'PolarAxesRepresentation'
    clipDisplay.ScalarOpacityUnitDistance = 42210192.22786043   
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    clipDisplay.DataAxesGrid.XTitleFontFile = ''
    clipDisplay.DataAxesGrid.YTitleFontFile = ''
    clipDisplay.DataAxesGrid.ZTitleFontFile = ''
    clipDisplay.DataAxesGrid.XLabelFontFile = ''
    clipDisplay.DataAxesGrid.YLabelFontFile = ''
    clipDisplay.DataAxesGrid.ZLabelFontFile = ''   
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    clipDisplay.PolarAxes.PolarAxisTitleFontFile = ''
    clipDisplay.PolarAxes.PolarAxisLabelFontFile = ''
    clipDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
    clipDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    return clipDisplay

def data_display(data, render_view):
    # set active source
    SetActiveSource(data)
    
    # show data in view
    data = Show(data, render_view)
    
    # trace defaults for the display properties.
    data.Representation = 'Surface'
    data.ColorArrayName = [None, '']
    data.OSPRayScaleArray = 'Te'
    data.OSPRayScaleFunction = 'PiecewiseFunction'
    data.SelectOrientationVectors = 'None'
    data.ScaleFactor = 0.35000000000000003
    data.SelectScaleArray = 'None'
    data.GlyphType = 'Arrow'
    data.GlyphTableIndexArray = 'None'
    data.GaussianRadius = 0.0175
    data.SetScaleArray = ['POINTS', 'Te']
    data.ScaleTransferFunction = 'PiecewiseFunction'
    data.OpacityArray = ['POINTS', 'Te']
    data.OpacityTransferFunction = 'PiecewiseFunction'
    data.DataAxesGrid = 'GridAxesRepresentation'
    data.SelectionCellLabelFontFile = ''
    data.SelectionPointLabelFontFile = ''
    data.PolarAxes = 'PolarAxesRepresentation'
    data.ScalarOpacityUnitDistance = 0.06443385322347264
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    data.ScaleTransferFunction.Points = [0.34215784072875977, 0.0, 0.5, 0.0, 2.905837297439575, 1.0, 0.5, 0.0]    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    data.OpacityTransferFunction.Points = [0.34215784072875977, 0.0, 0.5, 0.0, 2.905837297439575, 1.0, 0.5, 0.0]
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    data.DataAxesGrid.XTitleFontFile = ''
    data.DataAxesGrid.YTitleFontFile = ''
    data.DataAxesGrid.ZTitleFontFile = ''
    data.DataAxesGrid.XLabelFontFile = ''
    data.DataAxesGrid.YLabelFontFile = ''
    data.DataAxesGrid.ZLabelFontFile = ''
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    data.PolarAxes.PolarAxisTitleFontFile = ''
    data.PolarAxes.PolarAxisLabelFontFile = ''
    data.PolarAxes.LastRadialAxisTextFontFile = ''
    data.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    return data

master_dir = '/run/user/1000/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/v_jet/hd/2D'

## Input of changing vaibles.
#jet_angle = ['0.0','0.1','0.5','1.0','5.0','10.0','15.0','20.0','25.0','30.0']
jet_angle = ['0.1','0.5','1.0','5.0','10.0','15.0','20.0','25.0','30.0']
fps = 7
res = [1892, 1024]
rat = 1


for ii in range(len(jet_angle)):
    #creates folder to put movies in
    # save location
    data_path = master_dir+'/jet_a'+jet_angle[ii]
    file_sav_name = glob.glob1(data_path,'*.vtu')
    path = []
    for i in range(len(file_sav_name)):
        string = data_path+'/'+file_sav_name[i]
        path.append(string)
    
    data = XMLUnstructuredGridReader(FileName=path)
    
    # data names
    data.CellArrayStatus = ['Te', 'cs', 'e', 'm1', 'm2', 'p', 'rho', 'schrho', 'v1', 'v2']
    # first set of vids
    data_names = data.CellArrayStatus[5]    
    # get animation scene
    # Dont know what this does, something similar to below
    animationScene1 = GetAnimationScene()
    animationScene1.GoToLast()
    
    # update animation scene based on data timesteps
    # This dispays the nb of timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    render_view = render_view_func()
    
    # get layout
    layout1 = GetLayout()

    # place view in the layout
    layout1.AssignView(0, render_view)
    data_display(data, render_view)
    # reset view to fit data
    render_view.ResetCamera()
    #changing interaction mode based on data extents
    render_view.InteractionMode = '2D'    
    render_view.CameraPosition = [1.25, 0.75, 7.35627193011166]
    render_view.CameraFocalPoint = [1.25, 0.75, 0.0]
    render_view.CameraParallelScale = 1.3004188760781206

    # update the view to ensure updated data information
    render_view.Update()    #
    
   # get active source.
    xMLUnstructuredGridReader1 = GetActiveSource()
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1356, 713]
    
    # get display properties
    xMLUnstructuredGridReader1Display = GetDisplayProperties(xMLUnstructuredGridReader1, view=render_view)
    
    # set scalar coloring
    ColorBy(xMLUnstructuredGridReader1Display, ('POINTS', 'rho'))  

    # rescale color and/or opacity maps used to include current data range
    xMLUnstructuredGridReader1Display.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    xMLUnstructuredGridReader1Display.SetScalarBarVisibility(render_view, True)
    
    # get color transfer function/color map for 'rho'
    rhoLUT = GetColorTransferFunction('rho')
    
    # get color legend/bar for rhoLUT in view renderView1
    rhoLUTColorBar = GetScalarBar(rhoLUT, render_view)
    
    # change scalar bar placement
    rhoLUTColorBar.Orientation = 'Horizontal'
    rhoLUTColorBar.WindowLocation = 'AnyLocation'    
    # change scalar bar placement
    rhoLUTColorBar.Position = [0.3141740412979352, 0.06030855539971949]    
    # change scalar bar placement
    rhoLUTColorBar.ScalarBarLength = 0.46938053097345117    
    # change scalar bar placement
    rhoLUTColorBar.Position = [0.3053244837758113, 0.061711079943899017]
    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    rhoLUT.ApplyPreset('Linear Blue (8_31f)', True)
    # Properties modified on rhoLUTColorBar
    rhoLUTColorBar.LabelFontSize = 16    
    rhoLUTColorBar.TitleFontSize = 20

    # create a new 'Text'
    text1 = Text()
    text1.Text = file_sav_name[ii]
    # show data in view
    text1Display = Show(text1, render_view)
    # trace defaults for the display properties.
    text1Display.FontFile = ''
 
    
    sav_loc = master_dir+'/vids/jet_'+jet_angle[ii]
    if os.path.isdir(sav_loc) is False:
        os.makedirs(sav_loc)

    # save animation
    SaveAnimation(sav_loc+'\jet_a'+jet_angle[ii]+'.png', layout1, SaveAllViews=1, ImageResolution=[np.int(res[0]*rat), np.int(res[-1]*rat)], FrameRate=1, FrameWindow=[0, len(file_sav_name)])

    # destroy text1
    Delete(text1)
    del text1
    
    # set active source
    SetActiveSource(xMLUnstructuredGridReader1)
    
    # destroy xMLUnstructuredGridReader1
    Delete(xMLUnstructuredGridReader1)
    del xMLUnstructuredGridReader1
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    
    # destroy renderView1
    Delete(renderView1)
    del renderView1
    
    # get layout
    layout1 = GetLayoutByName("Layout #1")
    
    RemoveLayout(layout1)
    
    CreateLayout('Layout #1')
