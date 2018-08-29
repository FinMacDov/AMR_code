# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def render_view_func():
    renderView = CreateView('RenderView')
    renderView.ViewSize = [630, 1025] # careful: set to a particular value (fixes 3rd panel)
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
    
    renderView.CameraPosition = [0.0, 1000000000.0, 4700000000.0]
    renderView.CameraFocalPoint = [0.0, 1000000000.0, 0.0]
    renderView.CameraParallelScale = 1405561883.6231568

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

render_view = []

render_view.append(render_view_func())
# loads data
# create a new 'XML Unstructured Grid Reader'
data = XMLUnstructuredGridReader(FileName=['/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0000.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0001.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0002.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0003.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0004.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0005.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0006.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0007.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0008.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0009.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0010.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0011.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0012.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0013.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0014.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0015.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0016.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0017.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0018.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0019.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0020.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0021.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0022.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0023.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0024.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0025.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0026.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0027.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0028.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0029.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0030.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0031.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0032.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0033.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0034.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0035.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0036.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0037.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0038.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0039.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0040.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0041.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0042.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0043.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0044.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0045.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0046.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0047.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0048.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0049.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0050.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0051.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0052.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0053.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0054.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0055.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0056.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0057.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0058.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0059.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0060.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0061.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0062.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0063.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0064.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0065.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0066.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0067.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0068.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0069.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0070.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0071.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0072.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0073.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0074.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0075.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0076.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0077.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0078.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0079.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0080.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0081.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0082.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0083.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0084.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0085.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0086.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0087.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0088.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0089.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0090.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0091.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0092.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0093.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0094.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0095.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0096.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0097.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0098.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0099.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0100.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0101.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0102.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0103.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0104.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0105.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0106.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0107.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0108.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0109.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0110.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0111.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0112.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0113.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0114.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0115.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0116.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0117.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0118.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0119.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0120.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0121.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0122.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0123.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0124.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0125.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0126.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0127.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0128.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0129.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0130.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0131.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0132.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0133.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0134.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0135.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0136.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0137.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0138.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0139.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0140.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0141.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0142.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0143.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0144.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0145.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0146.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0147.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0148.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0149.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0150.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0151.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0152.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0153.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0154.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0155.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0156.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0157.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0158.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0159.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0160.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0161.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0162.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0163.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0164.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0165.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0166.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0167.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0168.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0169.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0170.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0171.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0172.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0173.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0174.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0175.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0176.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0177.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0178.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0179.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0180.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0181.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0182.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0183.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0184.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0185.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0186.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0187.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0188.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0189.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0190.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0191.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0192.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0193.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0194.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0195.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0196.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0197.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0198.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0199.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0200.vtu'])
# data names
data.CellArrayStatus = ['rho', 'v1', 'v2', 'v3', 'p', 'b1', 'b2', 'b3', 'trp1', 'Te', 'Alfv', 'divB', 'beta', 'schrho', 'cs']

# get animation scene
# Dont know what this does, something similar to below
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
# This dispays the nb of timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()
#displays data in paraview
dataDisplay = display_data(data)
# hide data in view
Hide(data, render_view[0])

#--------------------------------
# create a new 'Clip'
# !You will want to add file name for input!
clip1 = Clip(Input=jet_B30_V30_0)
clip1.ClipType = 'Plane'
clip1.Scalars = ['CELLS', 'Alfv']

# Clipping sdetting that we weant to change
clip1.ClipType.Origin = [0.0, 2000000000.0, 0.0]
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# Hides the plane outlines
Hide3DWidgets(proxy=clip1.ClipType)
# update the view to ensure updated data information
render_view[0].Update()

# makes sure three panels are the same size
size1 = 1.0/3.0
size2 = 0.5

# get the material library
materialLibrary1 = GetMaterialLibrary()

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(0, render_view[0])

# split cell
layout1.SplitHorizontal(0, size1)
layout1.SplitHorizontal(2, size2)

nb_panels = 3
for i in range(1,nb_panels):
    SetActiveView(None)
    render_view.append(render_view_func())
    layout1.AssignView(4+i, render_view[i])

data_names = [data.CellArrayStatus[0],data.CellArrayStatus[9],data.CellArrayStatus[2]] 
colour_wheel = ['Black-Body Radiation', 'Cool to Warm (Extended)', 'Cool to Warm (Extended)']
col_black = [0.0, 0.0, 0.0]
col_white = [1.0, 1.0, 1.0]
# find source
clip1 = FindSource('Clip1')
# set active source
SetActiveSource(clip1)
# SETING COLOR BAR PROPERTIES
colbar_pos_x = 0.85
colbar_pos_y = 0.125
scalar_bar_length = 0.75

clip_1display = []
rendLUT = []
rendLUTColorBar = []
for i in range(nb_panels):
    SetActiveView(render_view[i])
    clip_1display.append(clip_display(clip1,render_view[i]))    
    render_view[i].ResetCamera()
    # set scalar coloring
    ColorBy(clip_1display[i], ('CELLS', data_names[i]))
    # show color bar/color legend
    clip_1display[i].SetScalarBarVisibility(render_view[i], True)
    rendLUT.append(GetColorTransferFunction(data_names[i]))
    # where the color is change for vars
    rendLUT[i].ApplyPreset(colour_wheel[i], True) 
    
    rendLUTColorBar.append(GetScalarBar(rendLUT[i], render_view[i]))
    rendLUTColorBar[i].TitleColor = col_black
    rendLUTColorBar[i].LabelColor = col_black
    rendLUTColorBar[i].ScalarBarLength = scalar_bar_length
    rendLUTColorBar[i].Position = [colbar_pos_x, colbar_pos_y]

#
### Note: resize frame
##layout1.SetSplitFraction(0, 0.33333)


