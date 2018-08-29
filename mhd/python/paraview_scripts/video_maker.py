# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

render_view_func():
    renderView = CreateView('RenderView')
    renderView.ViewSize = [969, 415]
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
    return()


# loads data
# create a new 'XML Unstructured Grid Reader'
jet_B30_V30_0 = XMLUnstructuredGridReader(FileName=['/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0000.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0001.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0002.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0003.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0004.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0005.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0006.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0007.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0008.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0009.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0010.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0011.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0012.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0013.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0014.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0015.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0016.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0017.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0018.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0019.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0020.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0021.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0022.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0023.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0024.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0025.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0026.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0027.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0028.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0029.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0030.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0031.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0032.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0033.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0034.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0035.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0036.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0037.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0038.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0039.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0040.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0041.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0042.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0043.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0044.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0045.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0046.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0047.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0048.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0049.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0050.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0051.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0052.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0053.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0054.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0055.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0056.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0057.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0058.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0059.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0060.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0061.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0062.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0063.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0064.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0065.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0066.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0067.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0068.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0069.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0070.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0071.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0072.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0073.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0074.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0075.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0076.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0077.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0078.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0079.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0080.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0081.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0082.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0083.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0084.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0085.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0086.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0087.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0088.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0089.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0090.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0091.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0092.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0093.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0094.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0095.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0096.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0097.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0098.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0099.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0100.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0101.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0102.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0103.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0104.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0105.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0106.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0107.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0108.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0109.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0110.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0111.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0112.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0113.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0114.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0115.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0116.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0117.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0118.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0119.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0120.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0121.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0122.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0123.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0124.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0125.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0126.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0127.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0128.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0129.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0130.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0131.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0132.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0133.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0134.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0135.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0136.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0137.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0138.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0139.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0140.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0141.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0142.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0143.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0144.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0145.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0146.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0147.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0148.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0149.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0150.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0151.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0152.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0153.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0154.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0155.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0156.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0157.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0158.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0159.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0160.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0161.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0162.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0163.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0164.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0165.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0166.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0167.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0168.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0169.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0170.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0171.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0172.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0173.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0174.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0175.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0176.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0177.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0178.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0179.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0180.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0181.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0182.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0183.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0184.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0185.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0186.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0187.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0188.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0189.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0190.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0191.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0192.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0193.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0194.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0195.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0196.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0197.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0198.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0199.vtu', '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet/jet_B30_V30/jet_B30_V30_0200.vtu'])
# data names
jet_B30_V30_0.CellArrayStatus = ['rho', 'v1', 'v2', 'v3', 'p', 'b1', 'b2', 'b3', 'trp1', 'Te', 'Alfv', 'divB', 'beta', 'schrho', 'cs']

# get animation scene
# Dont know what this does, something similar to below
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
# This dispays the nb of timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [969, 415]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 0
renderView1.Background = [0.32, 0.34, 0.43]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.YTitle = '     Y Axis'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.XTitleBold = 1
renderView1.AxesGrid.XTitleFontSize = 15
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.YTitleBold = 1
renderView1.AxesGrid.YTitleFontSize = 15
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.ZTitleBold = 1
renderView1.AxesGrid.ZTitleFontSize = 15
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.XLabelBold = 1
renderView1.AxesGrid.XLabelFontSize = 15
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.YLabelBold = 1
renderView1.AxesGrid.YLabelFontSize = 15
renderView1.AxesGrid.ZLabelFontFile = ''
renderView1.AxesGrid.ZLabelBold = 1
renderView1.AxesGrid.ZLabelFontSize = 15
# uncomment following to set a specific view size
# renderView1.ViewSize = [1910, 1025]

# displays the data in paraview
jet_B30_V30_0Display = Show(jet_B30_V30_0, renderView1)

# trace defaults for the display properties.
jet_B30_V30_0Display.Representation = 'Surface'
jet_B30_V30_0Display.ColorArrayName = [None, '']
jet_B30_V30_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
jet_B30_V30_0Display.SelectOrientationVectors = 'None'
jet_B30_V30_0Display.ScaleFactor = 899999948.8000001
jet_B30_V30_0Display.SelectScaleArray = 'None'
jet_B30_V30_0Display.GlyphType = 'Arrow'
jet_B30_V30_0Display.GlyphTableIndexArray = 'None'
jet_B30_V30_0Display.GaussianRadius = 44999997.44
jet_B30_V30_0Display.SetScaleArray = [None, '']
jet_B30_V30_0Display.ScaleTransferFunction = 'PiecewiseFunction'
jet_B30_V30_0Display.OpacityArray = [None, '']
jet_B30_V30_0Display.OpacityTransferFunction = 'PiecewiseFunction'
jet_B30_V30_0Display.DataAxesGrid = 'GridAxesRepresentation'
jet_B30_V30_0Display.SelectionCellLabelFontFile = ''
jet_B30_V30_0Display.SelectionPointLabelFontFile = ''
jet_B30_V30_0Display.PolarAxes = 'PolarAxesRepresentation'
jet_B30_V30_0Display.ScalarOpacityUnitDistance = 167845592.2061122

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
jet_B30_V30_0Display.DataAxesGrid.XTitleFontFile = ''
jet_B30_V30_0Display.DataAxesGrid.YTitleFontFile = ''
jet_B30_V30_0Display.DataAxesGrid.ZTitleFontFile = ''
jet_B30_V30_0Display.DataAxesGrid.XLabelFontFile = ''
jet_B30_V30_0Display.DataAxesGrid.YLabelFontFile = ''
jet_B30_V30_0Display.DataAxesGrid.ZLabelFontFile = ''

# hide data in view
Hide(jet_B30_V30_0, renderView1)

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
renderView1.Update()

# makes sure three panels are the same size
size1 = 1.0/3.0
size2 = 0.5

# get the material library
materialLibrary1 = GetMaterialLibrary()

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(0, renderView1)

# split cell
layout1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# split cell
layout1.SplitHorizontal(2, 0.5)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [640, 480]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.StereoType = 0
renderView2.Background = [0.32, 0.34, 0.43]
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
layout1.AssignView(5, renderView2)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.ViewSize = [235, 415]
renderView3.AxesGrid = 'GridAxes3DActor'
renderView3.StereoType = 0
renderView3.Background = [0.32, 0.34, 0.43]
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

# set active view
SetActiveView(renderView1)

# splits the frames into three 
layout1.SetSplitFraction(0, size1)
layout1.SetSplitFraction(2, size2)

# find source
clip1 = FindSource('Clip1')

# set active source
SetActiveSource(clip1)

# show data in view
clip1Display = Show(clip1, renderView1)

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 200000000.0
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 10000000.0
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.SelectionCellLabelFontFile = ''
clip1Display.SelectionPointLabelFontFile = ''
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 42210192.22786043

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

# reset view to fit data
renderView1.ResetCamera()

# set active view
SetActiveView(renderView2)

# show data in view
clip1Display_1 = Show(clip1, renderView2)

# trace defaults for the display properties.
clip1Display_1.Representation = 'Surface'
clip1Display_1.ColorArrayName = [None, '']
clip1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display_1.SelectOrientationVectors = 'None'
clip1Display_1.ScaleFactor = 200000000.0
clip1Display_1.SelectScaleArray = 'None'
clip1Display_1.GlyphType = 'Arrow'
clip1Display_1.GlyphTableIndexArray = 'None'
clip1Display_1.GaussianRadius = 10000000.0
clip1Display_1.SetScaleArray = [None, '']
clip1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display_1.OpacityArray = [None, '']
clip1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display_1.DataAxesGrid = 'GridAxesRepresentation'
clip1Display_1.SelectionCellLabelFontFile = ''
clip1Display_1.SelectionPointLabelFontFile = ''
clip1Display_1.PolarAxes = 'PolarAxesRepresentation'
clip1Display_1.ScalarOpacityUnitDistance = 42210192.22786043

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

# set active view
SetActiveView(renderView3)

# show data in view
clip1Display_2 = Show(clip1, renderView3)

# trace defaults for the display properties.
clip1Display_2.Representation = 'Surface'
clip1Display_2.ColorArrayName = [None, '']
clip1Display_2.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display_2.SelectOrientationVectors = 'None'
clip1Display_2.ScaleFactor = 200000000.0
clip1Display_2.SelectScaleArray = 'None'
clip1Display_2.GlyphType = 'Arrow'
clip1Display_2.GlyphTableIndexArray = 'None'
clip1Display_2.GaussianRadius = 10000000.0
clip1Display_2.SetScaleArray = [None, '']
clip1Display_2.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display_2.OpacityArray = [None, '']
clip1Display_2.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display_2.DataAxesGrid = 'GridAxesRepresentation'
clip1Display_2.SelectionCellLabelFontFile = ''
clip1Display_2.SelectionPointLabelFontFile = ''
clip1Display_2.PolarAxes = 'PolarAxesRepresentation'
clip1Display_2.ScalarOpacityUnitDistance = 42210192.22786043

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

#set the camera views for each render
# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 1000000000.0, 4700000000.0]
renderView1.CameraFocalPoint = [0.0, 1000000000.0, 0.0]
renderView1.CameraParallelScale = 1405561883.6231568

# current camera placement for renderView2
renderView2.CameraPosition = [0.0, 1000000000.0, 4700000000.0]
renderView2.CameraFocalPoint = [0.0, 1000000000.0, 0.0]
renderView2.CameraParallelScale = 1405561883.6231568

# current camera placement for renderView3
renderView3.CameraPosition = [0.0, 1000000000.0, 4700000000.0]
renderView3.CameraFocalPoint = [0.0, 1000000000.0, 0.0]
renderView3.CameraParallelScale = 1405561883.6231568

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0
renderView2.OrientationAxesVisibility = 0
renderView3.OrientationAxesVisibility = 0

# set scalar coloring
ColorBy(clip1Display, ('CELLS', 'rho'))
ColorBy(clip1Display_1, ('CELLS', 'Te'))
ColorBy(clip1Display_2, ('CELLS', 'v2'))

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)
clip1Display_1.SetScalarBarVisibility(renderView2, True)
clip1Display_2.SetScalarBarVisibility(renderView3, True)

# SETING COLOR BAR PROPERTIES
colour_wheel = ['Black-Body', 'Cool to Warm', 'Cool to Warm (Extended)', 'Grayscale']
primvar = ['rho', 'Te', 'v2']
rend = []
LUTColorBar

test = [renderView1, ]

for i in range(len(primvar)):
     rend.append(GetColorTransferFunction(primvar[i]))
     LUTColorBar.append(GetScalarBar(rend[i], renderView1))     

# get color legend/bar for rhoLUT in view renderView1
rend1LUTColorBar = GetScalarBar(rend1LUT, renderView1)
rend2LUTColorBar = GetScalarBar(rend2LUT, renderView2)
rend3LUTColorBar = GetScalarBar(rend3LUT, renderView3)

col_black = [0.0, 0.0, 0.0]
col_white = [1.0, 1.0, 1.0]
# Properties modified on v2LUTColorBar
rend1LUTColorBar.TitleColor = col_black
rend2LUTColorBar.TitleColor = col_black
rend3LUTColorBar.TitleColor = col_black
# Properties modified on v2LUTColorBar
rend1LUTColorBar.LabelColor = col_black
rend2LUTColorBar.LabelColor = col_black
rend3LUTColorBar.LabelColor = col_black

# change scalar bar length
rend1LUTColorBar.ScalarBarLength = 0.75
rend2LUTColorBar.ScalarBarLength = 0.75
rend3LUTColorBar.ScalarBarLength = 0.75

# change scalar bar placement
colbar_pos_x = 0.85
colbar_pos_y = 0.125
rend1LUTColorBar.Position = [colbar_pos_x, colbar_pos_y]
rend2LUTColorBar.Position = [colbar_pos_x, colbar_pos_y]
rend3LUTColorBar.Position = [colbar_pos_x, colbar_pos_y]

## Note: resize frame
#layout1.SetSplitFraction(0, 0.33333)


