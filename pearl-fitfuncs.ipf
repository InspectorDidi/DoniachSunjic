#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma IgorVersion = 6.2
#pragma ModuleName = PearlFitFuncs
#pragma version = 1.01
#include "mm-physconst", version >= 1.05

// various fit functions for photoelectron spectroscopy

// $Id$
// author: matthias.muntwiler@psi.ch
// Copyright (c) 2013-14 Paul Scherrer Institut

// Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//   http://www.apache.org/licenses/LICENSE-2.0

//------------------------------------------------------------------------------
// Doniach-Sunjic fit functions
//------------------------------------------------------------------------------

threadsafe function DoniachSunjic(x, amp, pos, sing, fwhm)
	// Doniach-Sunjic line shape
	// [S. Doniach, M. Sunjic, J. Phys. C 3 (1970) 285]
	variable x	// independent variable
	variable amp // amplitude
	variable pos // position
	variable sing // singularity index (0 <= sing < 1)
	variable fwhm // width

	variable nom, denom
	nom = cos(pi * sing / 2 + (1 - sing) * atan((x - pos) / fwhm * 2))
	denom = ((x - pos)^2 + fwhm^2 / 4)^((1 - sing) / 2)
	
	return amp * nom / denom * fwhm / 2
end

threadsafe function ds1_bg(w, x): FitFunc
	// Doniach-Sunjic fit function
	// 0 <= sing < 1
	wave w	// coefficients - see below
	variable x	// independent variable

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = DoniachSunjic(x, amp, pos, sing, fwhm) + bg
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = bg
	//CurveFitDialog/ w[1] = amp
	//CurveFitDialog/ w[2] = pos
	//CurveFitDialog/ w[3] = sing
	//CurveFitDialog/ w[4] = FWHM

	return DoniachSunjic(x, w[1], w[2], w[3], w[4]) + w[0]	
end

threadsafe function ds2_bg(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = w_0+( w_1*cos(pi*w_3/2+(1-w_3)*atan((x-w_2)/w_4)))/(((x-w_2)^2+w_4^2)^((1-w_3)/2))   +(w_5*cos(pi*w_7/2+(1-w_7)*atan((x-(w_6))/w_8)))/(((x-w_6)^2+w_8^2)^((1-w_7)/2))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 9
	//CurveFitDialog/ w[0] = bg
	//CurveFitDialog/ w[1] = amp1
	//CurveFitDialog/ w[2] = pos1
	//CurveFitDialog/ w[3] = sing1
	//CurveFitDialog/ w[4] = wid1
	//CurveFitDialog/ w[5] = amp2
	//CurveFitDialog/ w[6] = pos2
	//CurveFitDialog/ w[7] = sing2
	//CurveFitDialog/ w[8] = wid2

	variable ds1 = DoniachSunjic(x, w[1], w[2], w[3], w[4])
	variable ds2 = DoniachSunjic(x, w[5], w[6], w[7], w[8])
	
	return w[0] + ds1 + ds2
End

Function ds4_bg(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = w_0+( w_1*cos(pi*w_3/2+(1-w_3)*atan((x-w_2)/w_4)))/(((x-w_2)^2+w_4^2)^((1-w_3)/2))   +(w_5*cos(pi*w_7/2+(1-w_7)*atan((x-(w_6))/w_8)))/(((x-w_6)^2+w_8^2)^((1-w_7)/2))   +( w_9*cos(pi*w_11/2+(1-w_11)*atan((x-w_10)/w_12)))/(((x-w_10)^2+w_12^2)^((1-w_11)/2))      +( w_13*cos(pi*w_15/2+(1-w_15)*atan((x-w_14)/w_16)))/(((x-w_14)^2+w_16^2)^((1-w_15)/2))   
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 17
	//CurveFitDialog/ w[0] = w_0
	//CurveFitDialog/ w[1] = w_11
	//CurveFitDialog/ w[2] = w_12
	//CurveFitDialog/ w[3] = w_13
	//CurveFitDialog/ w[4] = w_14
	//CurveFitDialog/ w[5] = w_21
	//CurveFitDialog/ w[6] = w_22
	//CurveFitDialog/ w[7] = w_23
	//CurveFitDialog/ w[8] = w_24
	//CurveFitDialog/ w[9] = w_31
	//CurveFitDialog/ w[10] = w_32
	//CurveFitDialog/ w[11] = w_33
	//CurveFitDialog/ w[12] = w_34
	//CurveFitDialog/ w[13] = w_41
	//CurveFitDialog/ w[14] = w_42
	//CurveFitDialog/ w[15] = w_43
	//CurveFitDialog/ w[16] = w_44
	Variable ds1, ds2, ds3, ds4
	ds1=( w[1]*cos(pi*w[3]/2+(1-w[3])*atan((x-w[2])/w[4])))/(((x-w[2])^2+w[4]^2)^((1-w[3])/2))
	ds2=( w[5]*cos(pi*w[7]/2+(1-w[7])*atan((x-w[6])/w[8])))/(((x-w[6])^2+w[8]^2)^((1-w[7])/2))
	ds3=( w[9]*cos(pi*w[11]/2+(1-w[11])*atan((x-w[10])/w[12])))/(((x-w[10])^2+w[12]^2)^((1-w[11])/2))
	ds4=( w[13]*cos(pi*w[15]/2+(1-w[15])*atan((x-w[14])/w[16])))/(((x-w[14])^2+w[16]^2)^((1-w[15])/2))

	
	return w[0]+ds1+ds2+ds3+ds4
	

End

Function ds6_bg(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ 
	//CurveFitDialog/ Variable g, ds1, ds2, ds3, ds4, ds5, ds6
	//CurveFitDialog/ ds1=( w_11*cos(pi*w_13/2+(1-w_13)*atan((x-w_12)/w_14)))/(((x-w_12)^2+w_14^2)^((1-w_13)/2))
	//CurveFitDialog/ ds2=( w_21*cos(pi*w_23/2+(1-w_23)*atan((x-w_22)/w_24)))/(((x-w_22)^2+w_24^2)^((1-w_23)/2))
	//CurveFitDialog/ ds3=( w_31*cos(pi*w_33/2+(1-w_33)*atan((x-w_32)/w_34)))/(((x-w_32)^2+w_34^2)^((1-w_33)/2))
	//CurveFitDialog/ ds4=( w_41*cos(pi*w_43/2+(1-w_43)*atan((x-w_42)/w_44)))/(((x-w_42)^2+w_44^2)^((1-w_43)/2))
	//CurveFitDialog/ ds5=( w_51*cos(pi*w_53/2+(1-w_53)*atan((x-w_52)/w_54)))/(((x-w_52)^2+w_54^2)^((1-w_53)/2))
	//CurveFitDialog/ ds6=( w_61*cos(pi*w_63/2+(1-w_63)*atan((x-w_62)/w_64)))/(((x-w_62)^2+w_64^2)^((1-w_63)/2))
	//CurveFitDialog/ 
	//CurveFitDialog/ f(x) =w_0+ds1+ds2+ds3+ds4+ds5+ds6
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 25
	//CurveFitDialog/ w[0] = w_0
	//CurveFitDialog/ w[1] = w_11
	//CurveFitDialog/ w[2] = w_12
	//CurveFitDialog/ w[3] = w_13
	//CurveFitDialog/ w[4] = w_14
	//CurveFitDialog/ w[5] = w_21
	//CurveFitDialog/ w[6] = w_22
	//CurveFitDialog/ w[7] = w_23
	//CurveFitDialog/ w[8] = w_24
	//CurveFitDialog/ w[9] = w_31
	//CurveFitDialog/ w[10] = w_32
	//CurveFitDialog/ w[11] = w_33
	//CurveFitDialog/ w[12] = w_34
	//CurveFitDialog/ w[13] = w_41
	//CurveFitDialog/ w[14] = w_42
	//CurveFitDialog/ w[15] = w_43
	//CurveFitDialog/ w[16] = w_44
	//CurveFitDialog/ w[17] = w_51
	//CurveFitDialog/ w[18] = w_52
	//CurveFitDialog/ w[19] = w_53
	//CurveFitDialog/ w[20] = w_54
	//CurveFitDialog/ w[21] = w_61
	//CurveFitDialog/ w[22] = w_62
	//CurveFitDialog/ w[23] = w_63
	//CurveFitDialog/ w[24] = w_64

	
	Variable ds1, ds2, ds3, ds4, ds5, ds6
	ds1=( w[1]*cos(pi*w[3]/2+(1-w[3])*atan((x-w[2])/w[4])))/(((x-w[2])^2+w[4]^2)^((1-w[3])/2))
	ds2=( w[5]*cos(pi*w[7]/2+(1-w[7])*atan((x-w[6])/w[8])))/(((x-w[6])^2+w[8]^2)^((1-w[7])/2))
	ds3=( w[9]*cos(pi*w[11]/2+(1-w[11])*atan((x-w[10])/w[12])))/(((x-w[10])^2+w[12]^2)^((1-w[11])/2))
	ds4=( w[13]*cos(pi*w[15]/2+(1-w[15])*atan((x-w[14])/w[16])))/(((x-w[14])^2+w[16]^2)^((1-w[15])/2))
	ds5=( w[17]*cos(pi*w[19]/2+(1-w[19])*atan((x-w[18])/w[20])))/(((x-w[18])^2+w[20]^2)^((1-w[19])/2))
	ds6=( w[21]*cos(pi*w[23]/2+(1-w[23])*atan((x-w[22])/w[24])))/(((x-w[22])^2+w[24]^2)^((1-w[23])/2))
	
	return w[0]+ds1+ds2+ds3+ds4+ds5+ds6
	
End

structure DoniachSunjicStruct
	// data structure for DoniachSunjicBroadS structural function fit

	// waves populated by the FuncFit operation	
	wave pw
	wave yw
	wave xw

	// convolution parameters to be set upon creation of the structure
	variable precision
	variable oversampling

	// auxiliary fields used internally by DoniachSunjicBroadS
	// do not touch these
	wave xdw
	wave model
	wave broadening
	wave convolution
EndStructure

//------------------------------------------------------------------------------
threadsafe function DoniachSunjicBroadS(s) : FitFunc
//------------------------------------------------------------------------------
	// Two-peak (bulk + surface) Doniach-Sunjic line shape with Gaussian broadening (convolution).
	// Hold parameter 5 at 0 to fit just one peak.
	
	// Structural fit function for efficient fitting in procedures.
	// Calculating the convolution requires auxiliary waves and additional, non-fitting parameters.
	// To eliminate the time-consuming overhead of creating and killing the auxiliary waves,
	// these waves are held in the fitting structure.
	
	// Caution: The function on its own is thread-safe.
	// However, since FuncFit uses the same structure in all threads, the fitting cannot run in parallel.
	// Set /NTHR=1.
	
	// See also Fit_DoniachSunjicBroad (example), DoniachSunjicBroad (conventional fit function)
	Struct DoniachSunjicStruct &s

	// pw[0] = bulk amplitude
	// pw[1] = bulk position
	// pw[2] = Lorentzian FWHM
	// pw[3] = Donjach-Sunjic singularity index (0..1)
	// pw[4] = surface shift
	// pw[5] = surface/bulk ratio
	// pw[6] = Gaussian FWHM
	// pw[7] = constant background
	// pw[8] = linear background
	
	wave xw = s.xw
	wave yw = s.yw
	wave pw = s.pw
	
	variable precision = s.precision
	variable oversampling = s.oversampling

	if (WaveExists(s.xdw))
		wave xdw = s.xdw
		wave model = s.model
		wave broadening = s.broadening
		wave convolution = s.convolution
	else
		make /n=0 /free xdw, model, broadening, convolution
		redimension /d xdw, model, broadening, convolution
		wave fs.xdw = xdw
		wave fs.model = model
		wave fs.broadening = broadening
		wave fs.convolution = convolution
	endif

	// calculate wave spacing based on minimum spacing of desired x points
	differentiate /p xw /d=xdw
	xdw = abs(xdw)
	variable xd = wavemin(xdw) / oversampling
	
	// calculate broadening wave size based on width and precision variable
	variable x0b = pw[6] * precision
	variable nb = 2 * floor(x0b / xd) + 1
	
	// calculate model size based on desired range for yw
	variable x0m = max(abs(wavemax(xw) - pw[1]), abs(wavemin(xw) - pw[1])) + x0b
	variable nm = 2 * floor(x0m / xd) + 1
	nb = min(nb, nm * 10)	// limit wave size to avoid runtime errors for unphysically large parameter
	
	// create and calculate initial waves, normalize exponential
	redimension /n=(nb) broadening
	redimension /n=(nm) model
	setscale/i x -x0b, x0b, "", broadening
	setscale/i x -x0m, x0m, "", model
	
	broadening = exp( - (x / pw[6])^2 * 4 * ln(2))
	variable nrm = area(broadening)
	broadening /= nrm
	model = DoniachSunjic(x, 1, 0, pw[3], pw[2]) // bulk
	model += DoniachSunjic(x, pw[5], pw[4], pw[3], pw[2]) // surface
	
	// calculate the convolution
	Convolve /a broadening, model
	variable scale = pw[0] / wavemax(model)
	model *= scale
	
	// prepare output
	nm = numpnts(model)
	x0m = xd * (nm - 1) / 2
	setscale/i x -x0m, x0m, "", model
	
	yw = model(xw[p] - pw[1]) + pw[7] + pw[8] * xw[p]
	yw = numtype(yw) ? 0 : yw
end

//------------------------------------------------------------------------------
function DoniachSunjicBroad(pw, yw, xw) : FitFunc
//------------------------------------------------------------------------------
	// Two-peak (bulk + surface) Doniach-Sunjic line shape with Gaussian broadening (convolution).
	// Hold parameter 5 at 0 to fit just one peak.
	// Conventional fit function for use with the curve-fitting dialog.
	// Compared to DoniachSunjicBroadS this function incurs extra overhead
	// because auxiliary waves are created and killed between function calls.
	// See also DoniachSunjicBroadS (optimized structural fit function)
	Wave pw
	Wave yw
	Wave xw

	// pw[0] = bulk amplitude
	// pw[1] = bulk position
	// pw[2] = Lorentzian FWHM
	// pw[3] = Donjach-Sunjic singularity index (0..1)
	// pw[4] = surface shift
	// pw[5] = surface/bulk ratio
	// pw[6] = Gaussian FWHM
	// pw[7] = constant background
	// pw[8] = linear background
	
	// set up data structure
	struct DoniachSunjicStruct fs
	fs.precision = 5
	fs.oversampling = 4
	
	wave fs.pw = pw
	wave fs.xw = xw
	wave fs.yw = yw
	
	// create temporary calculation waves in a global folder
	dfref df = root:packages:pearl_fitfuncs:doniach_sunjic
	if (DataFolderRefStatus(df) == 0)
		newdatafolder root:packages:pearl_fitfuncs:doniach_sunjic
		dfref df = root:packages:pearl_fitfuncs:doniach_sunjic
	endif
	
	wave /z /sdfr=df fs.xdw = xdw
	wave /z /sdfr=df fs.model = model
	wave /z /sdfr=df fs.broadening = broadening
	wave /z /sdfr=df fs.convolution = convolution

	if (WaveExists(fs.xdw) == 0)
		dfref savedf = GetDataFolderDFR()
		setdatafolder df
		make /n=0 /d xdw, model, broadening, convolution
		wave fs.xdw = xdw
		wave fs.model = model
		wave fs.broadening = broadening
		wave fs.convolution = convolution
		setdatafolder savedf
	endif

	// calculate
	DoniachSunjicBroadS(fs)
	
	yw = fs.yw
end

//------------------------------------------------------------------------------
function Calc_DoniachSunjicBroad(pw, yw)
//------------------------------------------------------------------------------
	// Calculate the DoniachSunjicBroadS line shape
	Wave pw // coefficient wave
	Wave yw // output wave, correct x-scaling required on input
	
	struct DoniachSunjicStruct fs
	fs.precision = 5
	fs.oversampling = 4
	
	duplicate /free pw, fs.pw
	duplicate /free yw, fs.xw
	fs.xw = x
	duplicate /free yw, fs.yw
	
	DoniachSunjicBroadS(fs)
	
	yw = fs.yw
end

//------------------------------------------------------------------------------
Function Fit_DoniachSunjicBroad(pw, yw, xw, ww)
//------------------------------------------------------------------------------
	// Fit the DoniachSunjicBroadS line shape.
	// The function applies constraints which assume that the energy scale is in eV.
	// Returns chi^2.
	wave pw // coefficient wave- pre-load it with initial guess
	wave yw
	wave /z xw
	wave /z ww // weights (standard deviation)
	
	struct DoniachSunjicStruct fs
	fs.precision = 5
	fs.oversampling = 4
	
	duplicate /free pw, fs.pw
	if (WaveExists(xw))
		duplicate /free xw, fs.xw
	else
		duplicate /free yw, fs.xw
		fs.xw = x
	endif
	duplicate /free yw, fs.yw
	
	variable v_chisq = nan
	variable V_FitMaxIters = 100
	make /n=1 /t /free constraints = {"K0 >= 0", "K2 > 0", "K2 < 10", "K3 >= 0", "K3 < 1", "K4 >= -10", "K4 <= 10", "K5 >= 0", "K5 <= 1", "K6 >= 0", "K6 < 10"}
	// note: only single thread allowed
	FuncFit /NTHR=1 DoniachSunjicBroadS, pw, yw /X=xw /D /STRC=fs /C=constraints /NWOK /I=1 /W=ww
	
	return v_chisq
End

//------------------------------------------------------------------------------
// peak-specific fit functions
//------------------------------------------------------------------------------

function Au4f(w, x): fitfunc
	// fit function for a nitrogen 1s-pi* absorption spectrum
	// modelled as multiple Voigt shapes on a constant background
	// similar to the Igor VoigtFit function
	// but with a constant gaussian width (instrumental broadening) for all peaks
	// gaussian and lorentzian widths are specified as FWHM
	wave w // parameters
		// w[0] constant background
		// w[1] linear background
		// w[2] global gaussian FWHM
		// w[3 + 0 + (n-1) * 3] peak n area
		// w[3 + 1 + (n-1) * 3] peak n position
		// w[3 + 2 + (n-1) * 3] peak n lorentzian FWHM
		// length of wave defines number of peaks
		
		// for compatibility with older code the linear background term can be omitted.
		// if the number of parameters divides by 3, the linear background term is added,
		// otherwise only the constant background.
	variable x
	
	variable np = 15
	variable ip, ip0
	
	variable bg = w[0]
	variable v = bg
	if (mod(np, 3) == 0)
		v += w[1] * x
		ip0 = 3
	else
		ip0 = 2
	endif

	variable vc1, vc2, vc3, vc4
	vc2 = 2 * sqrt(ln(2)) / w[ip0-1]
	for (ip = ip0; ip < np; ip += 3)
		vc1 = w[ip] / sqrt(pi) * vc2
		vc3 = w[ip+1]
		vc4 = vc2 * w[ip+2] / 2
		v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)
	endfor

	return v
	
end

function Au4f_2p2(w, x): fitfunc
	// Au 4f 5/2 and 7/2 2-component Voigt fit with a common gaussian width
	// gaussian and lorentzian widths are specified as FWHM
	wave w // parameters
		// w[0] constant background
		// w[1] linear background
		// w[2] global gaussian FWHM
		// w[3] 5/2 bulk area
		// w[4] 5/2 bulk position
		// w[5] 5/2 lorentzian FWHM
		// w[6] 7/2 bulk area
		// w[7] 7/2 bulk position
		// w[8] 7/2 lorentzian FWHM
		// w[9] surface/bulk area ratio
		// w[10] surface core level shift
	variable x

	variable bg = w[0] + w[1] * x
	variable v = bg

	variable vc1 // amplitude
	variable vc2 // width
	variable vc3 // position
	variable vc4 // shape
	vc2 = 2 * sqrt(ln(2)) / w[2]

	// 5/2 bulk
	vc1 = w[3] / sqrt(pi) * vc2
	vc3 = w[4]
	vc4 = vc2 * w[5] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 5/2 surface
	vc1 = w[3] / sqrt(pi) * vc2 * w[9]
	vc3 = w[4] + w[10]
	vc4 = vc2 * w[5] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 7/2 bulk
	vc1 = w[6] / sqrt(pi) * vc2
	vc3 = w[7]
	vc4 = vc2 * w[8] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 7/2 surface
	vc1 = w[6] / sqrt(pi) * vc2 * w[9]
	vc3 = w[7] + w[10]
	vc4 = vc2 * w[8] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	return v
	
end

function ShowComponents_Au4f_2p2(coef_wave, fit_wave)
	wave coef_wave
	wave fit_wave

	duplicate /free coef_wave, coef1, coef2
	coef1[9] = 0
	coef2[3] *= coef_wave[9]
	coef2[4] += coef_wave[10]
	coef2[6] *= coef_wave[9]
	coef2[7] += coef_wave[10]
	coef2[9] = 0
	
	string s_fit_wave = NameOfWave(fit_wave)
	string s_fit_p1 = s_fit_wave + "_p1"
	string s_fit_p2 = s_fit_wave + "_p2"
	duplicate /o fit_wave, $(s_fit_p1) /wave=fit_p1
	duplicate /o fit_wave, $(s_fit_p2) /wave=fit_p2
	
	fit_p1 = Au4f_2p2(coef1, x)
	fit_p2 = Au4f_2p2(coef2, x)

	string traces = TraceNameList("", ";", 1)
	if ((WhichListItem(s_fit_wave, traces, ";") >= 0) && (WhichListItem(s_fit_p1, traces, ";") < 0))
		appendtograph fit_p1, fit_p2
		ModifyGraph lstyle($s_fit_p1)=2
		ModifyGraph lstyle($s_fit_p2)=2
		ModifyGraph rgb($s_fit_p1)=(0,0,65280)
		ModifyGraph rgb($s_fit_p2)=(0,0,65280)
	endif
end

function Au4f_2p3(w, x): fitfunc
	// Au 4f 5/2 and 7/2 3-component Voigt fit with a common gaussian width
	// gaussian and lorentzian widths are specified as FWHM
	wave w // parameters
		// w[0] constant background
		// w[1] linear background
		// w[2] global gaussian FWHM
		// w[3] 5/2 bulk area
		// w[4] 5/2 bulk position
		// w[5] 5/2 lorentzian FWHM
		// w[6] 7/2 bulk area
		// w[7] 7/2 bulk position
		// w[8] 7/2 lorentzian FWHM
		// w[9] surface/bulk area ratio
		// w[10] surface core level shift
		// w[11] 2nd layer/bulk area ratio
		// w[12] 2nd layer core level shift
	variable x

	variable bg = w[0] + w[1] * x
	variable v = bg

	variable vc1 // amplitude
	variable vc2 // width
	variable vc3 // position
	variable vc4 // shape
	vc2 = 2 * sqrt(ln(2)) / w[2]

	// 5/2 bulk
	vc1 = w[3] / sqrt(pi) * vc2
	vc3 = w[4]
	vc4 = vc2 * w[5] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 5/2 surface
	vc1 = w[3] / sqrt(pi) * vc2 * w[9]
	vc3 = w[4] + w[10]
	vc4 = vc2 * w[5] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 5/2 2nd layer
	vc1 = w[3] / sqrt(pi) * vc2 * w[11]
	vc3 = w[4] + w[12]
	vc4 = vc2 * w[5] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 7/2 bulk
	vc1 = w[6] / sqrt(pi) * vc2
	vc3 = w[7]
	vc4 = vc2 * w[8] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 7/2 surface
	vc1 = w[6] / sqrt(pi) * vc2 * w[9]
	vc3 = w[7] + w[10]
	vc4 = vc2 * w[8] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	// 7/2 2nd layer
	vc1 = w[6] / sqrt(pi) * vc2 * w[11]
	vc3 = w[7] + w[12]
	vc4 = vc2 * w[8] / 2
	v += vc1 * VoigtFunc(vc2 * (x - vc3), vc4)

	return v
	
end

function ShowComponents_Au4f_2p3(coef_wave, fit_wave)
	wave coef_wave
	wave fit_wave

	duplicate /free coef_wave, coef1, coef2, coef3
	coef1[9] = 0
	coef1[11] = 0
	
	coef2[3] *= coef_wave[9]
	coef2[4] += coef_wave[10]
	coef2[6] *= coef_wave[9]
	coef2[7] += coef_wave[10]
	coef2[9] = 0
	coef2[11] = 0
	
	coef3[3] *= coef_wave[11]
	coef3[4] += coef_wave[12]
	coef3[6] *= coef_wave[11]
	coef3[7] += coef_wave[12]
	coef3[9] = 0
	coef3[11] = 0
	
	string s_fit_wave = NameOfWave(fit_wave)
	string s_fit_p1 = s_fit_wave + "_p1"
	string s_fit_p2 = s_fit_wave + "_p2"
	string s_fit_p3 = s_fit_wave + "_p3"
	duplicate /o fit_wave, $(s_fit_p1) /wave=fit_p1
	duplicate /o fit_wave, $(s_fit_p2) /wave=fit_p2
	duplicate /o fit_wave, $(s_fit_p3) /wave=fit_p3
	
	fit_p1 = Au4f_2p2(coef1, x)
	fit_p2 = Au4f_2p2(coef2, x)
	fit_p3 = Au4f_2p2(coef3, x)

	string traces = TraceNameList("", ";", 1)
	if ((WhichListItem(s_fit_wave, traces, ";") >= 0) && (WhichListItem(s_fit_p1, traces, ";") < 0))
		appendtograph fit_p1, fit_p2, fit_p3
		ModifyGraph lstyle($s_fit_p1)=2
		ModifyGraph lstyle($s_fit_p2)=2
		ModifyGraph lstyle($s_fit_p3)=2
		ModifyGraph rgb($s_fit_p1)=(0,0,65280)
		ModifyGraph rgb($s_fit_p2)=(0,0,65280)
		ModifyGraph rgb($s_fit_p3)=(0,0,65280)
	endif
end

/// convolution of Fermi-Dirac distribution and a Gaussian.
///
/// @arg pw[0] = constant background
/// @arg pw[1] = linear background
/// @arg pw[2] = amplitude
/// @arg pw[3] = Fermi level in eV
/// @arg pw[4] = temperature in K
/// @arg pw[5] = gaussian width = FWHM / 1.66511
///
function FermiGaussConv(pw, yw, xw) : FitFunc
	WAVE pw, yw, xw
	
	// half width of temporary gaussian wave is pw[5] multiplied by this factor (may be fractional)
	variable precision_g = 5
	variable oversampling = 4
	
	// calculate wave spacing based on minimum spacing of desired x points
	duplicate /free xw, xdw
	differentiate /p xw /d=xdw
	xdw = abs(xdw)
	variable xd = wavemin(xdw) / oversampling
	
	// calculate gausswave size based on pw[5] and precision variable
	variable x0g = abs(pw[5]) * precision_g
	variable ng = 2 * floor(x0g / xd) + 1
	
	// calculate fermiwave size based on desired range for yw
	variable emax = wavemax(xw)
	variable emin = wavemin(xw)
	variable x0f = max(abs(emax - pw[3]), abs(emin - pw[3])) + x0g
	variable ne = 2 * floor(x0f / xd) + 1
	
	// create and calculate initial waves, normalize exponential
	make /d /n=(ng) /free gausswave
	make /d /n=(ne) /free fermiwave
	setscale/i x -x0g, x0g, "", gausswave
	setscale/i x -x0f, x0f, "", fermiwave
	
	gausswave = exp( - (x / pw[5] )^2 )
	fermiwave = 1 / (exp( x / (kBoltzmann * pw[4])) + 1.0 )

	// calculate the convolution
	duplicate /free fermiwave, resultwave
	Convolve /a gausswave, resultwave
	variable rmax = wavemax(resultwave)
	resultwave /= rmax
	
	// prepare output
	ng = numpnts(resultwave)
	x0g = xd * (ng - 1) / 2
	setscale/i x -x0g, x0g, "", resultwave
	
	yw = pw[2] * resultwave(xw[p] - pw[3]) + pw[0] + pw[1] * xw[p]
end
