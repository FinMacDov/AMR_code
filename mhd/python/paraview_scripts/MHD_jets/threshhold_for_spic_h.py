# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:45:50 2018

@author: fionnlagh
"""

# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
jetB30V30_0 = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [720, 841]

# get display properties
jetB30V30_0Display = GetDisplayProperties(jetB30V30_0, view=renderView1)

# set scalar coloring
ColorBy(jetB30V30_0Display, ('CELLS', 'rho'))

# rescale color and/or opacity maps used to include current data range
jetB30V30_0Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
jetB30V30_0Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'rho'
rhoLUT = GetColorTransferFunction('rho')

# create a new 'Threshold'
threshold1 = Threshold(Input=jetB30V30_0)
threshold1.Scalars = ['CELLS', 'Alfv']
threshold1.ThresholdRange = [226593.28125, 70748952.0]

# get animation scene
animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

# Properties modified on threshold1
threshold1.Scalars = ['CELLS', 'trp1']
threshold1.ThresholdRange = [2.0, 100.0]

# show data in view
threshold1Display = Show(threshold1, renderView1)

# get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['CELLS', 'rho']
threshold1Display.LookupTable = rhoLUT
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 18933332.8
threshold1Display.SelectScaleArray = 'None'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'None'
threshold1Display.GaussianRadius = 946666.64
threshold1Display.SetScaleArray = [None, '']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = [None, '']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.SelectionCellLabelFontFile = ''
threshold1Display.SelectionPointLabelFontFile = ''
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityFunction = rhoPWF
threshold1Display.ScalarOpacityUnitDistance = 14309320.26929398

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
threshold1Display.DataAxesGrid.XTitleFontFile = ''
threshold1Display.DataAxesGrid.YTitleFontFile = ''
threshold1Display.DataAxesGrid.ZTitleFontFile = ''
threshold1Display.DataAxesGrid.XLabelFontFile = ''
threshold1Display.DataAxesGrid.YLabelFontFile = ''
threshold1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
threshold1Display.PolarAxes.PolarAxisTitleFontFile = ''
threshold1Display.PolarAxes.PolarAxisLabelFontFile = ''
threshold1Display.PolarAxes.LastRadialAxisTextFontFile = ''
threshold1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(jetB30V30_0, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 4499999744.0, 20205472460.030487]
renderView1.CameraFocalPoint = [0.0, 4499999744.0, 0.0]
renderView1.CameraParallelScale = 5276974277.364563

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).