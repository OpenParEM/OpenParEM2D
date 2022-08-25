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

# Real

et_reLUT = GetColorTransferFunction('Et_re')
et_reLUT.NumberOfTableValues = 64

Et_re_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Et real magnitude', Et_re_calculator)
Et_re_calculatorDisplay = Show(Et_re_calculator, EM2DView)
Et_re_calculatorDisplay.Representation = 'Surface'
Et_re_calculatorDisplay.ColorArrayName = ['POINTS', 'Et_re']
Et_re_calculatorDisplay.LookupTable = et_reLUT
Et_re_calculator.Function = 'Et_re_X*iHat+Et_re_Y*jHat+Ez_re*kHat'

Et_re_glyph = Glyph(Input=Et_re_calculator,GlyphType='2D Glyph')
RenameSource("E real vector",Et_re_glyph)
Et_re_glyphDisplay = Show(Et_re_glyph, EM2DView)
ColorBy(Et_re_glyphDisplay, None)
Et_re_glyph.GlyphType = '2D Glyph'
Et_re_glyph.OrientationArray = ['POINTS', 'Result']
Et_re_glyph.ScaleArray = ['POINTS', 'Result']
Et_re_glyph.ScaleFactor = 1e-9

# Imaginary

et_imLUT = GetColorTransferFunction('Et_im')
et_imLUT.NumberOfTableValues = 64

Et_im_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Et imag magnitude', Et_im_calculator)
Et_im_calculatorDisplay = Show(Et_im_calculator, EM2DView)
Et_im_calculatorDisplay.Representation = 'Surface'
Et_im_calculatorDisplay.ColorArrayName = ['POINTS', 'Et_im']
Et_im_calculatorDisplay.LookupTable = et_imLUT
Et_im_calculator.Function = 'Et_im_X*iHat+Et_im_Y*jHat+Ez_im*kHat'

Et_im_glyph = Glyph(Input=Et_im_calculator,GlyphType='2D Glyph')
RenameSource("E imag vector",Et_im_glyph)
Et_im_glyphDisplay = Show(Et_im_glyph, EM2DView)
ColorBy(Et_im_glyphDisplay, None)
Et_im_glyph.GlyphType = '2D Glyph'
Et_im_glyph.OrientationArray = ['POINTS', 'Result']
Et_im_glyph.ScaleArray = ['POINTS', 'Result']
Et_im_glyph.ScaleFactor = 1e-9

#############################################################
# Ez fields
#############################################################

# Real

ez_reLUT = GetColorTransferFunction('Ez_re')
ez_reLUT.NumberOfTableValues = 64

Ez_re_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Ez real', Ez_re_calculator)
Ez_re_calculatorDisplay = Show(Ez_re_calculator, EM2DView)
Ez_re_calculatorDisplay.Representation = 'Surface'
Ez_re_calculatorDisplay.ColorArrayName = ['POINTS', 'Ez_re']
Ez_re_calculatorDisplay.LookupTable = ez_reLUT
Ez_re_calculator.Function = 'Ez_re'

# Imaginary

ez_imLUT = GetColorTransferFunction('Ez_im')
ez_imLUT.NumberOfTableValues = 64

Ez_im_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Ez imag', Ez_im_calculator)
Ez_im_calculatorDisplay = Show(Ez_im_calculator, EM2DView)
Ez_im_calculatorDisplay.Representation = 'Surface'
Ez_im_calculatorDisplay.ColorArrayName = ['POINTS', 'Ez_im']
Ez_im_calculatorDisplay.LookupTable = ez_imLUT
Ez_im_calculator.Function = 'Ez_im'

#############################################################
# Ht fields
#############################################################

# Real

ht_reLUT = GetColorTransferFunction('Ht_re')
ht_reLUT.NumberOfTableValues = 64

Ht_re_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Ht real magnitude', Ht_re_calculator)
Ht_re_calculatorDisplay = Show(Ht_re_calculator, EM2DView)
Ht_re_calculatorDisplay.Representation = 'Surface'
Ht_re_calculatorDisplay.ColorArrayName = ['POINTS', 'Ht_re']
Ht_re_calculatorDisplay.LookupTable = ht_reLUT
Ht_re_calculator.Function = 'Ht_re_X*iHat+Ht_re_Y*jHat+Hz_re*kHat'

Ht_re_glyph = Glyph(Input=Ht_re_calculator,GlyphType='2D Glyph')
RenameSource("H real vector",Ht_re_glyph)
Ht_re_glyphDisplay = Show(Ht_re_glyph, EM2DView)
ColorBy(Ht_re_glyphDisplay, None)
Ht_re_glyph.GlyphType = '2D Glyph'
Ht_re_glyph.OrientationArray = ['POINTS', 'Result']
Ht_re_glyph.ScaleArray = ['POINTS', 'Result']

# Imaginary

ht_imLUT = GetColorTransferFunction('Ht_im')
ht_imLUT.NumberOfTableValues = 64

Ht_im_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Ht imag magnitude', Ht_im_calculator)
Ht_im_calculatorDisplay = Show(Ht_im_calculator, EM2DView)
Ht_im_calculatorDisplay.Representation = 'Surface'
Ht_im_calculatorDisplay.ColorArrayName = ['POINTS', 'Ht_im']
Ht_im_calculatorDisplay.LookupTable = ht_imLUT
Ht_im_calculator.Function = 'Ht_im_X*iHat+Ht_im_Y*jHat+Hz_re*kHat'

Ht_im_glyph = Glyph(Input=Ht_im_calculator,GlyphType='2D Glyph')
RenameSource("H imag vector",Ht_im_glyph)
Ht_im_glyphDisplay = Show(Ht_im_glyph, EM2DView)
ColorBy(Ht_im_glyphDisplay, None)
Ht_im_glyph.GlyphType = '2D Glyph'
Ht_im_glyph.OrientationArray = ['POINTS', 'Result']
Ht_im_glyph.ScaleArray = ['POINTS', 'Result']

#############################################################
# Hz fields
#############################################################

# Real

hz_reLUT = GetColorTransferFunction('Hz_re')
hz_reLUT.NumberOfTableValues = 64

Hz_re_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Hz real', Hz_re_calculator)
Hz_re_calculatorDisplay = Show(Hz_re_calculator, EM2DView)
Hz_re_calculatorDisplay.Representation = 'Surface'
Hz_re_calculatorDisplay.ColorArrayName = ['POINTS', 'Hz_re']
Hz_re_calculatorDisplay.LookupTable = hz_reLUT
Hz_re_calculator.Function = 'Hz_re'

# Imaginary

hz_imLUT = GetColorTransferFunction('Hz_im')
hz_imLUT.NumberOfTableValues = 64

Hz_im_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Hz imag', Hz_im_calculator)
Hz_im_calculatorDisplay = Show(Hz_im_calculator, EM2DView)
Hz_im_calculatorDisplay.Representation = 'Surface'
Hz_im_calculatorDisplay.ColorArrayName = ['POINTS', 'Hz_im']
Hz_im_calculatorDisplay.LookupTable = hz_imLUT
Hz_im_calculator.Function = 'Hz_im'

#############################################################
# Pz fields
#############################################################

# Real

pz_reLUT = GetColorTransferFunction('Pz_re')
pz_reLUT.NumberOfTableValues = 64

Pz_re_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Pz real magnitude', Pz_re_calculator)
Pz_re_calculatorDisplay = Show(Pz_re_calculator, EM2DView)
Pz_re_calculatorDisplay.Representation = 'Surface'
Pz_re_calculatorDisplay.ColorArrayName = ['POINTS', 'Pz_re']
Pz_re_calculatorDisplay.LookupTable = pz_reLUT
Pz_re_calculator.Function = ''

# Imaginary

pz_imLUT = GetColorTransferFunction('Pz_im')
pz_imLUT.NumberOfTableValues = 64

Pz_im_calculator = Calculator(Input=EM2Dpvd)
RenameSource('Pz imag magnitude', Pz_im_calculator)
Pz_im_calculatorDisplay = Show(Pz_im_calculator, EM2DView)
Pz_im_calculatorDisplay.Representation = 'Surface'
Pz_im_calculatorDisplay.ColorArrayName = ['POINTS', 'Pz_im']
Pz_im_calculatorDisplay.LookupTable = pz_imLUT
Pz_im_calculator.Function = ''

#############################################################
# window operations
#############################################################

# Et magnitude window - hide

et_im_magnitude = FindSource('Et imag magnitude')
SetActiveSource(et_im_magnitude)
Hide(et_im_magnitude,EM2DView)

# E vector window - hide

e_im_vector = FindSource('E imag vector')
SetActiveSource(e_im_vector)
Hide(e_im_vector,EM2DView)

# Ez - hide

ez_re = FindSource('Ez real')
SetActiveSource(ez_re)
Hide(ez_re,EM2DView)

ez_im = FindSource('Ez imag')
SetActiveSource(ez_im)
Hide(ez_im,EM2DView)

# Ht magnitude window - hide

ht_re_magnitude = FindSource('Ht real magnitude')
SetActiveSource(ht_re_magnitude)
Hide(ht_re_magnitude,EM2DView)

ht_im_magnitude = FindSource('Ht imag magnitude')
SetActiveSource(ht_im_magnitude)
Hide(ht_im_magnitude,EM2DView)

# Hz - hide

hz_re = FindSource('Hz real')
SetActiveSource(hz_re)
Hide(hz_re,EM2DView)

hz_im = FindSource('Hz imag')
SetActiveSource(hz_im)
Hide(hz_im,EM2DView)

# H vector window - hide

h_re_vector = FindSource('H real vector')
SetActiveSource(h_re_vector)
Hide(h_re_vector,EM2DView)

h_im_vector = FindSource('H imag vector')
SetActiveSource(h_im_vector)
Hide(h_im_vector,EM2DView)

# Pz magnitude window - hide

pz_re_magnitude = FindSource('Pz real magnitude')
SetActiveSource(pz_re_magnitude)
Hide(pz_re_magnitude,EM2DView)

pz_im_magnitude = FindSource('Pz imag magnitude')
SetActiveSource(pz_im_magnitude)
Hide(pz_im_magnitude,EM2DView)

# E vector window - set color

e_re_vector = FindSource('E real vector')
SetActiveSource(e_re_vector)
e_re_vectorDisplay = GetDisplayProperties(e_re_vector, view=EM2DView)
e_re_vectorDisplay.AmbientColor = [0.16, 0.16, 0.16]
e_re_vectorDisplay.DiffuseColor = [0.16, 0.16, 0.16]

e_im_vector = FindSource('E imag vector')
SetActiveSource(e_im_vector)
e_im_vectorDisplay = GetDisplayProperties(e_im_vector, view=EM2DView)
e_im_vectorDisplay.AmbientColor = [0.16, 0.16, 0.16]
e_im_vectorDisplay.DiffuseColor = [0.16, 0.16, 0.16]

#############################################################

EM2DView.Update()

