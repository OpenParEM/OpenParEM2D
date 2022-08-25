from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

EM2Dpvd = GetActiveSource() 
EM2Dpvd.CellArrays = ['attribute']
EM2Dpvd.PointArrays = ['Et_im', 'Et_re', 'Ez_im', 'Ez_re', 'Ht_im', 'Ht_re', 'Hz_im', 'Hz_re', 'Pz_re', 'Pz_im']

EM2DView = GetActiveViewOrCreate('RenderView')
Hide(EM2Dpvd, EM2DView)

#############################################################
# Et fields
#############################################################

et_reLUT = GetColorTransferFunction('Et_re')
et_reLUT.NumberOfTableValues = 64

Et_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Et real magnitude', Et_calculator)
Et_calculatorDisplay = Show(Et_calculator, EM2DView)
Et_calculatorDisplay.Representation = 'Surface'
Et_calculatorDisplay.ColorArrayName = ['POINTS', 'Et_re']
Et_calculatorDisplay.LookupTable = et_reLUT
Et_calculator.Function = 'Et_re_X*iHat+Et_re_Y*jHat'

Et_glyph = Glyph(Input=Et_calculator,GlyphType='2D Glyph')
RenameSource("Et real vector",Et_glyph)
Et_glyphDisplay = Show(Et_glyph, EM2DView)
ColorBy(Et_glyphDisplay, None)
Et_glyph.GlyphType = '2D Glyph'
Et_glyph.OrientationArray = ['POINTS', 'Result']
Et_glyph.ScaleArray = ['POINTS', 'Result']
#Et_glyph.ScaleFactor = 1e-9

#############################################################
# Ht fields
#############################################################

ht_reLUT = GetColorTransferFunction('Ht_re')
ht_reLUT.NumberOfTableValues = 64

Ht_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Ht real magnitude', Ht_calculator)
Ht_calculatorDisplay = Show(Ht_calculator, EM2DView)
Ht_calculatorDisplay.Representation = 'Surface'
Ht_calculatorDisplay.ColorArrayName = ['POINTS', 'Ht_re']
Ht_calculatorDisplay.LookupTable = ht_reLUT
Ht_calculator.Function = 'Ht_re_X*iHat+Ht_re_Y*jHat'

Ht_glyph = Glyph(Input=Ht_calculator,GlyphType='2D Glyph')
RenameSource("Ht real vector",Ht_glyph)
Ht_glyphDisplay = Show(Ht_glyph, EM2DView)
ColorBy(Ht_glyphDisplay, None)
Ht_glyph.GlyphType = '2D Glyph'
Ht_glyph.OrientationArray = ['POINTS', 'Result']
Ht_glyph.ScaleArray = ['POINTS', 'Result']

#############################################################
# Pz fields
#############################################################

pz_reLUT = GetColorTransferFunction('Pz_re')
pz_reLUT.NumberOfTableValues = 64

Pz_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Pz real magnitude', Pz_calculator)
Pz_calculatorDisplay = Show(Pz_calculator, EM2DView)
Pz_calculatorDisplay.Representation = 'Surface'
Pz_calculatorDisplay.ColorArrayName = ['POINTS', 'Pz_re']
Pz_calculatorDisplay.LookupTable = pz_reLUT
Pz_calculator.Function = ''

#############################################################
# window operations
#############################################################

# Et magnitude window

# Ht magnitude window - hide
htmagnitude = FindSource('Ht real magnitude')
SetActiveSource(htmagnitude)
Hide(htmagnitude,EM2DView)

# Pz magnitude window - hide
pzmagnitude = FindSource('Pz real magnitude')
SetActiveSource(pzmagnitude)
Hide(pzmagnitude,EM2DView)

# Et vector window - set color
etvector = FindSource('Et real vector')
SetActiveSource(etvector)
etvectorDisplay = GetDisplayProperties(etvector, view=EM2DView)
etvectorDisplay.AmbientColor = [0.16, 0.16, 0.16]
etvectorDisplay.DiffuseColor = [0.16, 0.16, 0.16]

#############################################################

EM2DView.Update()

