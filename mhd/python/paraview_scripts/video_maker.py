# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
import glob
import img2vid as i2v
import os
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

colour_wheel = ['Black-Body Radiation', 'Cool to Warm', 'Cool to Warm (Extended)']
col_black = [0.0, 0.0, 0.0]
col_white = [1.0, 1.0, 1.0]

def render_view_func():
    renderView = CreateView('RenderView')
    renderView.ViewSize = [630, 1025] # careful: set to a particular value (fixes 3rd panel)
    renderView.AxesGrid = 'GridAxes3DActor'
    renderView.StereoType = 0
    renderView.Background = [0.32, 0.34, 0.43]
#    materialLibrary1 = GetMaterialLibrary()
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

def movie_maker(mini_path, B, V, mini_sav_loc):
    # mini_Path path up till var changes
    # magnetic feild strength
    # Velocity
    # mini_sav_loc - path up till var changes
    B_str = str(np.int(B))
    V_str = str(np.int(V))
    file_sav_name = '/jet_B'+B_str+'_V'+V_str+'_'
    render_view = []

    render_view.append(render_view_func())
    # loads data
    # create a new 'XML Unstructured Grid Reader'
    myPath = mini_path+'/jet_B'+B_str+'_V'+V_str
    vtuCounter = len(glob.glob1(myPath,'*.vtu'))
    path = []
    for i in range(vtuCounter):
        string = myPath+file_sav_name
        number =  str(i).zfill(4)
        path.append(string + number + '.vtu')
    
    data = XMLUnstructuredGridReader(FileName=path)
    
    # data names
    data.CellArrayStatus = ['rho', 'v1', 'v2', 'v3', 'p', 'b1', 'b2', 'b3', 'trp1', 'Te', 'Alfv', 'divB', 'beta', 'schrho', 'cs']
    
    data_names = [data.CellArrayStatus[0],data.CellArrayStatus[9],data.CellArrayStatus[2]] 
    
    # get animation scene
    # Dont know what this does, something similar to below
    animationScene1 = GetAnimationScene()
    animationScene1.GoToLast()
    
    # update animation scene based on data timesteps
    # This dispays the nb of timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    
    #--------------------------------
    # create a new 'Clip'
    # !You will want to add file name for input!
    name = 'jet_B'+B_str+'_V'+V_str+'_0'
    clip1 = Clip(Input=data)
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
    
    # find source
    clip1 = FindSource('Clip1')
    # set active source
    SetActiveSource(clip1)
    # SETING COLOR BAR PROPERTIES
    colbar_pos_x = 0.85
    colbar_pos_y = 0.12
    scalar_bar_length = 0.75
    scalar_bar_thickness = 15
    
    clip_1display = []
    rendLUT = []
    rendLUTColorBar = []
    rendPWF = []
    # goes to last time step
    for i in range(nb_panels):
        SetActiveView(render_view[i])
        clip_1display.append(clip_display(clip1,render_view[i]))    
    #    render_view[i].ResetCamera()
        # set scalar coloring
        ColorBy(clip_1display[i], ('CELLS', data_names[i]))
        # show color bar/color legend
        clip_1display[i].SetScalarBarVisibility(render_view[i], True)
        rendLUT.append(GetColorTransferFunction(data_names[i]))
        # Makes parameter log scale
        if data_names[i] == 'rho':
            rendLUT[i].MapControlPointsToLogSpace()
            rendLUT[i].UseLogScale = 1
        rendPWF.append(GetOpacityTransferFunction(data_names[i]))
        # where the color is change for vars
        rendLUT[i].ApplyPreset(colour_wheel[i], True) 
        rendLUTColorBar.append(GetScalarBar(rendLUT[i], render_view[i]))
        # colour bar text col, pos and font size
        rendLUTColorBar[i].TitleColor = col_black
        rendLUTColorBar[i].LabelColor = col_black
        rendLUTColorBar[i].ScalarBarLength = scalar_bar_length
        rendLUTColorBar[i].ScalarBarThickness = scalar_bar_thickness
        rendLUTColorBar[i].Position = [colbar_pos_x, colbar_pos_y]
        rendLUTColorBar[i].TitleFontSize = 14
        rendLUTColorBar[i].LabelFontSize = 13
        rendLUTColorBar[i].HorizontalTitle = 1
        rendLUTColorBar[i].TitleJustification = 'Left'
        if data_names[i] == 'v2': # or 'v1' or 'b1' or 'b2':
            small_data = data.CellData.GetArray('v2').GetRange()
            v_range = np.max(abs(np.asarray(small_data)))
            rendLUT[i].RescaleTransferFunction(-v_range, v_range)
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    
    animationScene1.GoToFirst()
    
    # create a new 'Annotate Time Filter'
    annotateTimeFilter1 = AnnotateTimeFilter(Input=clip1)
    annotateTimeFilter1.Format = 'Time: %g'
    annotateTimeFilter1.Scale = 2.14683
    
    # show data in view
    annotateTimeFilter1Display = Show(annotateTimeFilter1, render_view[-1])
    # trace defaults for the display properties.
    annotateTimeFilter1Display.FontFile = ''
    # Properties modified on annotateTimeFilter1Display
    annotateTimeFilter1Display.WindowLocation = 'LowerRightCorner'
    
    # create a new 'Text'
    text1 = Text()
    # Properties modified on text1
    text1.Text = 'V = '+V_str+' km s-1, B = ' +B_str+' G'
    text1Display = Show(text1, render_view[0])
    text1Display.FontFile = ''
    
    res = [1892, 1024]
    rat = 1
    
    if os.path.isdir(sav_loc+'/jet_B'+B_str+'_V'+V_str) is False:
        os.makedirs(sav_loc+'/jet_B'+B_str+'_V'+V_str)
    # save animation
    SaveAnimation(sav_loc+'/jet_B'+B_str+'_V'+V_str+file_sav_name+'.png', layout1, SaveAllViews=1, ImageResolution=[np.int(res[0]*rat), np.int(res[-1]*rat)], FrameRate=2, FrameWindow=[0, vtuCounter-2])
    
    print('Finished')
    
    #for i in range(nb_panels):
    del render_view
    
    RemoveLayout(layout1)
    
    CreateLayout('Layout #1')
    
    # find source
    annotateTimeFilter1 = FindSource('AnnotateTimeFilter1')
    
    # destroy annotateTimeFilter1
    Delete(annotateTimeFilter1)
    del annotateTimeFilter1
    
    # find source
    clip1 = FindSource('Clip1')
    
    # destroy clip1
    Delete(clip1)
    del clip1
    
    # get active source.
    text1 = GetActiveSource()
    
    # destroy text1
    Delete(text1)
    del text1
    
    # find source
    xMLUnstructuredGridReader1 = FindSource('XMLUnstructuredGridReader1')
    
    # destroy xMLUnstructuredGridReader1
    Delete(xMLUnstructuredGridReader1)
    del xMLUnstructuredGridReader1
    
    # get animation scene
    animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
## not using
## Properties modified on v2LUTColorBar
#rendLUTColorBar[-1].UseCustomLabels = 1
#label_values = np.floor(np.linspace(-v_range,v_range, 13))
#
## This to round the number nicely for vis
#OFFSET = 4
#for j in range(len(label_values)):
#    digit_count = -len(str(label_values[j]))
#    print(digit_count, label_values[j])
#    label_values[j] = round(label_values[j], digit_count+OFFSET)
#    
#rendLUTColorBar[-1].CustomLabels = label_values


#rendPWF[-1].RescaleTransferFunction(-v_range,v_range)
#
### Note: resize frame
##layout1.SetSplitFraction(0, 0.33333)


mini_path = '/run/user/1001/gvfs/smb-share:server=uosfstore.shef.ac.uk,share=shared/mhd_jet1/User/smp16fm/sims/jet'
sav_loc = mini_path+'/vids'

## Input of changing vaibles.
B = ['30', '40', '50', '60', '70', '80']
V = ['30', '40', '50', '60']
fps = 7

for ii in range(len(B)):
    for jj in range(len(V)):
        movie_maker(mini_path=mini_path, B=B[ii], V=V[jj], mini_sav_loc=sav_loc)
        # not sure if that . should be int he name
#        file_sav_name = 'jet_B'+str(B[ii])+'_V'+str(V[jj])+'_'
#        i2v.image2video(filepath=sav_loc+'/jet_B'+str(B[ii])+'_V'+str(V[jj])+'.', prefix=file_sav_name, in_extension='png', 
#                        output_name=prefix+'_video', out_extension='avi', 
#                        fps=fps, n_loops=1, delete_images=True, 
#                        delete_old_videos=True, res=1080, overlay=False, cover_page=False)
#        

