#pragma IndependentModule= Localizer
#pragma TextEncoding = "UTF-8"		// For details execute DisplayHelpTopic "The TextEncoding Pragma"
#pragma rtGlobals=3

// Copyright 2008-2015 Peter Dedecker.
//
// This file is part of Localizer.
//
// Localizer is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Localizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Localizer.  If not, see <http://www.gnu.org/licenses/>.


#pragma IgorVersion = 6.30
// Igor-related constants
constant kMaxWaveName = 31
constant kUserAbort = 57

// constants for the dimensions of the interface
constant kViewerGraph_EdgeMargin = 10
constant kViewerControls_Width = 477
constant kViewerControls_Height = 40
constant kViewerControls_Margin = 10
constant kAnalysisControls_Width = 284
constant kAnalysisControls_Height = 315
constant kAnalysisControls_Margin = 10
constant kHistogram_Margin = 6
constant kLUTControls_Height = 110
constant kTabHeight = 20

// camera types (for data loading)
constant CAMERA_TYPE_WINSPEC = 0
constant CAMERA_TYPE_ANDOR = 1
constant CAMERA_TYPE_HAMAMATSU = 2
constant CAMERA_TYPE_TIFF = 3
constant CAMERA_TYPE_PDE = 4	// a custom, very simple image format.
constant CAMERA_TYPE_ZEISS = 5		// Zeiss .lsm files. Currently unused since these are really TIFF files
constant CAMERA_TYPE_IGOR_WAVE = 6	// load data from an Igor wave
constant CAMERA_TYPE_MULTIFILE_TIFF = 8	// 7 is reserved for matlab matrices

// threshold methods
constant THRESHOLD_METHOD_GLRT = 0
constant THRESHOLD_METHOD_ISODATA = 1
constant THRESHOLD_METHOD_TRIANGLE = 2
constant THRESHOLD_METHOD_DIRECT = 3
constant THRESHOLD_METHOD_SMOOTHSIGMA = 4

// preprocessing
constant  PREPROCESSOR_NONE = 0
constant  PREPROCESSOR_3X3MEDIAN = 1
constant  PREPROCESSOR_5X5MEDIAN = 2
constant  PREPROCESSOR_1X1GAUSSIAN = 3
constant  PREPROCESSOR_2X2GAUSSIAN = 4
constant  PREPROCESSOR_3X3MEAN = 5
constant  PREPROCESSOR_5X5MEAN = 6

// postprocessing
constant POSTPROCESSOR_NONE = 0
constant POSTPROCESSOR_REMOVE_ISOLATED = 1

// localization methods
constant LOCALIZATION_GAUSS_FITTING = 0
constant LOCALIZATION_GAUSS_FITTING_FIX = 1
constant LOCALIZATION_MULTIPLICATION = 2
constant LOCALIZATION_CENTROID = 3
constant LOCALIZATION_ZEISSPALM = 4
constant LOCALIZATION_ELLIPSOIDAL2DGAUSS = 5
constant LOCALIZATION_MLEWG = 6
constant LOCALIZATION_ELLIPSGAUSS_ASTIG = 7
constant LOCALIZATION_ASTIG_3D = 8

// particle finding methods
constant PARTICLEFINDER_ADJACENT4 = 0
constant PARTICLEFINDER_ADJACENT8 = 1
constant PARTICLEFINDER_RADIUS = 2

// particle verification methods
constant PARTICLEVERIFIER_NONE = 0
constant PARTICLEVERIFIER_2DGAUSS = 1
constant PARTICLEVERIFIER_ELLIPS_SYMM = 2
constant PARTICLEVERIFIER_OVERLAP = 3
constant PARTICLEVERIFIER_ELLIPS_ASTIG = 4
constant PARTICLEVERIFIER_FIXEDWIDTH = 5

// CCD processing methods
constant PROCESSCCDIMAGE_SUBAVG = 0
constant PROCESSCCDIMAGE_DIFFIMAGE = 1
constant PROCESSCCDIMAGE_CONVERTFORMAT = 2
constant PROCESSCCDIMAGE_CROP = 3
constant PROCESSCCDIMAGE_CONVERTPHOTONS = 4

// CCD processing output formats
constant IMAGE_OUTPUT_TYPE_MOV = -1
constant IMAGE_OUTPUT_TYPE_TIFF = 0
constant IMAGE_OUTPUT_TYPE_COMPR_TIFF = 1
constant IMAGE_OUTPUT_TYPE_IGOR = 2
constant IMAGE_OUTPUT_TYPE_PDE = 3
constant IMAGE_OUTPUT_TYPE_MULTITIFF = 4

// CCD images analysis methods
constant ANALYZECCD_SUMMEDINTENSITYTRACE = 0
constant ANALYZECCD_AVERAGETRACE = 1
constant ANALYZECCD_AVERAGE_IMAGE = 2
constant ANALYZECCD_VARIANCE_IMAGE = 3
constant ANALYZECCD_BLEACHING_ANALYSIS = 4

// PALM bitmap error estimators
constant PALMBITMAP_DEV_SAME = 0
constant PALMBITMAP_DEV_FITUNCERTAINTY = 1
constant PALMBITMAP_DEV_GAUSSIANMASK = 2

// structure that makes it more easy to pass data to e.g. the fitting routines
structure FitData
	Variable directThresholdLevel
	Variable cameraType
	variable background		// currently unused
	variable currentImage
	variable minDistanceBetweenParticles	// currently unused
	variable PFA
	variable smoothSigmaFactor
	variable PSFWidth
	variable xSize	// currently unused
	variable ySize	// currently unused
	variable numberOfImages
	variable firstFrameToAnalyze
	variable lastFrameToAnalyze
	variable thresholdMethod
	variable preprocessing
	variable postprocessing
	variable particlefinder
	variable localizationMethod
	string CCDFilePath
	
	// use or don't use a particular particle verification algorithm
	// to not use it set the corresponding value to zero
	// to use it set it to the corresponding constant defined above
	variable PVerOverlap
	variable PVerSymm
	variable PVerEllipseSymm
	variable PVerEllipseAstig
endstructure

// structure used to specify a progress reporter function to the operations
Constant kLocalizerProgStructVersion = 1
Structure LocalizerProgStruct
	uint32 version
	FUNCREF ProgressFunctionPrototype func
EndStructure

// included for ease of writing and reading .pde files if required
// (provided for convenience; not currently used in this Igor code)
Structure PDEFormatHeader
	uint32 magic			// set to 27
	uint32 version			// set to 1
	uint32 nImages
	uint32 xSize
	uint32 ySize
	uint32 storageType		// see immediately below
EndStructure
// use the following constants to specify the storage type
// (this is an excerpt from the XOP source code)
//const int STORAGE_TYPE_INT4 = 0;
//const int STORAGE_TYPE_UINT4 = 1;
//const int STORAGE_TYPE_INT8 = 2;
//const int STORAGE_TYPE_UINT8 = 3;
//const int STORAGE_TYPE_INT16 = 4;
//const int STORAGE_TYPE_UINT16 = 5;
//const int STORAGE_TYPE_INT32 = 6;
//const int STORAGE_TYPE_UINT32 = 7;
//const int STORAGE_TYPE_INT64 = 8;
//const int STORAGE_TYPE_UINT64 = 9;
//const int STORAGE_TYPE_FP32 = 10;
//const int STORAGE_TYPE_FP64 = 11;

Function WritePositionsToFile(pos, refNum)
	wave pos
	variable refNum
	
	// write the metadata
	fprintf refnum, GetPositionsSaveTextHeader(pos)
	
	// write the actual positions
	variable nPos = DimSize(pos, 0), nColumns = DimSize(pos, 1), i, j
	for (i = 0; i < nPos; i+=1)
		for (j = 0; j < nColumns; j+=1)
			if (j < (nColumns - 1))
				fprintf refNum, "%.16g\t", pos[i][j]
			else
				fprintf refNum, "%.16g\n", pos[i][j]
			endif
		endfor
	endfor
End

Function WriteTracksToFile(W_Tracks, refNum)
	wave /WAVE W_Tracks
	variable refNum
	
	variable nTracks = DimSize(W_Tracks, 0)
	variable i
	fprintf refNum, "Particle tracks assembled using Localizer\n"
	fprintf refNum, "Contains %d tracks\n", nTracks
	for (i = 0; i < nTracks; i+=1)
		fprintf refNum, "TRACK %d FOLLOWS\n", i
		wave thisPos = W_Tracks[i]
		WritePositionsToFile(thisPos, refNum)
		fprintf refNum, "END TRACK %d\n", i
	endfor
End

Function /S GetPositionsSaveTextHeader(pos)
	wave pos
	
	// construct a header that will identify the filetype and columns
	string header, columnLabels
	variable localizationMethod = NumberByKey("LOCALIZATION METHOD", note(pos))
	if (NumType(localizationMethod) == 2)
		localizationMethod = LOCALIZATION_GAUSS_FITTING
	endif
	
	switch (localizationMethod)
		case LOCALIZATION_GAUSS_FITTING:
			header = "Localized positions using symmetric 2D Gauss fitting\n"
			columnLabels = "First frame\tIntegrated intensity\tFitted PSF standard deviation\tX position (pixel)\tY position (pixel)\tBackground\tIntensity deviation\tPSF width deviation\tX position deviation\tY position deviation\tBackground deviation\tNumber of frames where this emitter is present\n"
			break
		case LOCALIZATION_GAUSS_FITTING_FIX:
			header = "Localized positions using symmetric 2D Gauss fitting with fixed PSF width\n"
			columnLabels = "First frame\tIntegrated intensity\tX position (pixel)\tY position (pixel)\tBackground\tIntensity deviation\tX position deviation\tY position deviation\tBackground deviation\tNumber of frames where this emitter is present\n"
			break
		case LOCALIZATION_MULTIPLICATION:
			header = "Localized positions using iterative multiplication\n"
			columnLabels = "First frame\tUsed PSF standard deviation\tX position (pixel)\tY position (pixel)\tNumber of frames where this emitter is present\n"
			break
		case LOCALIZATION_CENTROID:
			header = "Localized positions using centroid calculation\n"
			columnLabels = "First frame\tX position (pixel)\tY position (pixel)\tNumber of frames where this emitter is present\n"
			break
		case LOCALIZATION_ZEISSPALM:
			header = "Positions localized with the Zeiss PALM system\n"
			columnLabels = "First frame\tIntegrated number of photons\tX position (pixel)\tY position (pixel)\tLocalization precision (pixel)\tNumber of frames where this emitter is present\n"
			break
		case LOCALIZATION_ELLIPSOIDAL2DGAUSS:
		case LOCALIZATION_ELLIPSGAUSS_ASTIG:
			header = "Localized positions using ellipsoidal 2D Gauss fitting\n"
			columnLabels = "First frame\tIntegrated intensity\tFitted PSF standard deviation along x (pixel)\tFitted PSF standard deviation along y (pixel)\tX position (pixel)\tY position (pixel)\tCorrelation between x and y\t"
			columnLabels += "Background\tIntensity deviation\tPSF width deviation along x (pixel)\tPSF width deviation along y (pixel)\tX position deviation (pixel)\tY position deviation (pixel)\tCorrelation deviation\tBackground deviation\tNumber of frames where this emitter is present\n"
			break
		case LOCALIZATION_MLEWG:
			header = "Localized positions using MLEwG localization\n"
			columnLabels = "First frame\tIntegrated intensity\tFitted PSF standard deviation\tX position (pixel)\tY position (pixel)\tBackground\tPosition deviation\t\tNumber of frames where this emitter is present\n"
			break
		default:
			Abort "Unknown localization method"
	EndSwitch
	
	string fullHeader = header
	header += "IGOR WAVENOTE FOLLOWS\n"
	header += note(pos) + "\n"
	header += "DATA FOLLOWS\n"
	header += columnLabels
	return header
End

Function LoadPositionsFromTextFile()
	// load PALM positions from a text file
	// currently only the Zeiss output format is supported
	
	string inputFilePaths
	variable refNum
	
	// allow the user to select a series of text files
	Open /D/R/F="Text Files:.txt; All Files:.*" /M="Select one or more text files" /MULT=1 refNum
	inputFilePaths = S_fileName
	if (strlen(inputFilePaths) == 0)	// cancelled
		return 0
	endif
	
	// for each of the files, check if it is a Zeiss file or a file created by this software
	string currentFilePath, singleLine
	string ZeissPALMFilepaths = ""
	string thisSoftwareFiles = ""
	variable i
	
	for (i = 0; i < ItemsInList(inputFilePaths, "\r"); i+=1)
		currentFilePath = StringFromList(i, inputFilePaths, "\r")
		Open /Z=1 /R refNum as currentFilePath
		if (V_flag != 0)
			// an error occurred opening this file
			Printf "An error occurred opening the file at %s, skipping\r", currentFilePath
			continue
		endif
		
		// check if it is a valid Zeiss PALM file by looking for the unusual [mm] units specified
		FReadLine refNum, singleLine
		Close refNum
		if ((StringMatch(singleLine, "*[mm]*") == 1) || (StringMatch(singleLine, "*[nm]*") == 1))
			// probably a Zeiss PALM file
			ZeissPALMFilepaths += currentFilePath + "\r"
			continue
		endif
		
		// check if it is a file created by this software
		if (StringMatch(singleLine, "*Localized positions using*") == 1)
			// probably created by this software
			thisSoftwareFiles += currentFilePath + "\r"
			continue
		endif
		
		// unknown or unsupported file type, ignore
		Printf "The file at %s does not seem to be a file containing PALM positions or is not supported, skipping\r", currentFilePath
		i -= 1
		inputFilePaths = RemoveFromList(inputFilePaths, currentFilePath, "\r")
	endfor
	
	// load the Zeiss files
	if (ItemsInList(ZeissPALMFilepaths, "\r") > 0)
		LoadZeissPALMPositions(ZeissPALMFilepaths)
	endif
	
	// load files created with this software
	if (ItemsInList(thisSoftwareFiles, "\r") > 0)
		LoadXOPPositionsFiles(thisSoftwareFiles)
	endif
End	
			
Function LoadZeissPALMPositions(ZeissPALMFilepaths)
	// load PALM positions from the text output format used by the Zeiss software
	// Zeiss works in units of nanometers, but since we work in pixels we need to convert
	
	// the zeiss positions have their own custom positions format with 6 columns
	// signaled by LOCALIZATION METHOD:LOCALIZATION_ZEISSPALM
	
	// the files to load should be specified as a list separated by carriage returns (\r)
	string ZeissPALMFilepaths
	
	Assert(strlen(ZeissPALMFilepaths) > 0)
	
	variable nDataFiles, pixelSize = 100, nPositions, i, j, k, nEntries, offset, nFrames
	variable xSize = 256, ySize = 256
	string currentMacintoshFilePath, waveOutputName, promptString
	string waveNote, singleLine
	variable refNum, toNmConversionFactor
	
	NewDataFolder /O root:'Localized Positions'
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	
	nDataFiles = ItemsInList(ZeissPALMFilepaths, "\r")
	
	for (i = 0; i < nDataFiles; i+=1)
		currentMacintoshFilePath = ParseFilePath(5, StringFromList(i, ZeissPALMFilepaths, "\r"), ":", 0, 0)
		waveOutputName = ParseFilePath(3, currentMacintoshFilePath, ":", 0, 0)
		
		sprintf promptString, "Loading positions from %s", currentMacintoshFilePath
		Prompt pixelSize, "Optical CCD pixel size (in nm)"
		Prompt xSize, "X-size (in pixels)"
		Prompt ySize, "Y-size (in pixels)"
		DoPrompt promptString, pixelSize, xSize, ySize
		if (V_flag == 1)	// cancel
			continue
		endif
		
		// read the first line of the file to check the units that it uses (nm or mm)
		Open /R refNum as currentMacintoshFilePath
		FReadLine refNum, singleLine
		Close refNum
		
		if (StringMatch(singleLine, "*[mm]*") == 1)
			toNmConversionFactor = 1e6
		elseif (StringMatch(singleLine, "*[nm]*") == 1)
			toNmConversionFactor = 1
		else
			toNmConversionFactor = 1
		endif
		
		LoadWave /Q/A/G/O/B="N='_skip_';N=ZeissFirstFrame,T=96;N=ZeissNFrames,T=96;N=ZeissFramesMissing;N=ZeissXPosition,T=4; N=ZeissYPosition,T=4; N=ZeissPrecision,T=4;N=ZeissPhotons,T=96;" currentMacintoshFilePath
		
		wave ZeissFirstFrame, ZeissNFrames, ZeissFramesMissing, ZeissXPosition, ZeissYPosition, ZeissPrecision, ZeissPhotons
		nEntries = DimSize(ZeissFirstFrame, 0)
		
		// The Zeiss output format does not sort the positions according to the ascending frame number in which a positions first appears
		// some routines here assume that the positions are ordered
		// do the sorting now
		Sort ZeissFirstFrame, ZeissFirstFrame, ZeissNFrames, ZeissFramesMissing, ZeissXPosition, ZeissYPosition, ZeissPrecision, ZeissPhotons
		
		nPositions = 0
		for (j = 0; j < nEntries; j += 1)
			if ((ZeissFramesMissing[j] != 0) || (ZeissPhotons[j] == 0))
							// the 'ZeissFramesMissing' column is unusual, so ignore it for now
			continue		// occasionally the Zeiss file contains strange data positions containing zero photons
			endif
			nPositions += 1
		endfor
		
		Make /D/O/N=(nPositions,6) root:'Localized Positions':$waveOutputName
		wave outputWave = root:'Localized Positions':$waveOutputName
		FastOP outputWave = 0
		
		offset = 0
		for (j = 0; j < nEntries; j+=1)
			if ((ZeissFramesMissing[j] != 0) || (ZeissPhotons[j] == 0))
				continue
			endif
			
			outputWave[offset][0] = ZeissFirstFrame[j]
			outputWave[offset][1] = ZeissPhotons[j]
			outputWave[offset][2] = ZeissXPosition[j] / pixelSize * toNmConversionFactor	// positions are reported in mm by the Zeiss software
			outputWave[offset][3] = ZeissYPosition[j] / pixelSize * toNmConversionFactor
			outputWave[offset][4] = ZeissPrecision[j] / pixelSize	// precision is reported in nm
			outputWave[offset][5] = ZeissNFrames[j]
			
			offset += 1
		endfor
		
		// append a note with as much additional information as possible
		waveNote = "LOCALIZATION METHOD:" + num2str(LOCALIZATION_ZEISSPALM) + ";"
		waveNote += "ORIGINAL FILE PATH:" + currentMacintoshFilePath + ";"
		waveNote += "X SIZE:" + num2str(xSize) + ";"
		waveNote += "Y SIZE:" + num2str(ySize) + ";"
		waveNote += "NUMBER OF IMAGES:" + num2str(nFrames) + ";"
		waveNote += "X PIXEL SIZE:" + num2str(pixelSize) + ";"
		waveNote += "Y PIXEL SIZE:" + num2str(pixelSize) + ";"
		
		Note /K outputWave	// remove a previous note (should not occur)
		Note /NOCR outputWave, waveNote
	endfor
		
		
	KillWaves /Z ZeissFirstFrame, ZeissNFrames, ZeissFramesMissing, ZeissXPosition, ZeissYPosition, ZeissPrecision, ZeissPhotons
End

Function LoadXOPPositionsFiles(PALMFilePaths)
	// load PALM positions from the text output format used by the Zeiss software
	// Zeiss works in units of nanometers, but since we work in pixels we need to convert
	
	// the zeiss positions have their own custom positions format with 6 columns
	// signaled by LOCALIZATION METHOD:LOCALIZATION_ZEISSPALM
	
	// the files to load should be specified as a list separated by carriage returns (\r)
	string PALMFilePaths
	
	Assert(strlen(PALMFilePaths) > 0)
	
	variable refNum, nDataFiles, pixelSize = 100, nPositions, i, j, k, offset, nFrames
	variable firstLineWithData
	variable localizationMethod
	string currentMacintoshFilePath, waveOutputName, promptString, singleLine, waveNote
	DFREF savDF = GetDataFolderDFR()
	
	NewDataFolder /O root:'Localized Positions'
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	
	nDataFiles = ItemsInList(PALMFilePaths, "\r")
	
	for (i = 0; i < nDataFiles; i+=1)
		currentMacintoshFilePath = ParseFilePath(5, StringFromList(i, PALMFilePaths, "\r"), ":", 0, 0)
		
		waveOutputName = CleanupName(ParseFilePath(3, currentMacintoshFilePath, ":", 0, 0), 0)
		if (strlen(waveOutputName) > 27)
			// loadwave easily complains of a name that's too long, so guard for that
			waveOutputName = waveOutputName[0, 26]
		endif
		
		waveNote = ""
		
		// get some information from the header
		Open /R refNum as currentMacintoshFilePath
		FReadLine refNum, singleLine
		
		// determine the type of positions from the first line
		if (StringMatch(singleLine, "*symmetric 2D Gauss fitting with fixed PSF width*") == 1)
			localizationMethod = LOCALIZATION_GAUSS_FITTING_FIX
		elseif (StringMatch(singleLine, "*symmetric 2D Gauss fitting*") == 1)
			localizationMethod = LOCALIZATION_GAUSS_FITTING
		elseif (StringMatch(singleLine, "*ellipsoidal 2D Gauss fitting*") == 1)
			localizationMethod = LOCALIZATION_ELLIPSOIDAL2DGAUSS
		elseif (StringMatch(singleLine, "*centroid calculation*") == 1)
			localizationMethod = LOCALIZATION_CENTROID
		elseif (StringMatch(singleLine, "*iterative multiplication*") == 1)
			localizationMethod = LOCALIZATION_MULTIPLICATION
		elseif (StringMatch(singleLine, "*MLEwG localization*") == 1)
			localizationMethod = LOCALIZATION_MLEWG
		elseif (StringMatch(singleLine, "*ellipsoidal 2D Gauss fitting for astigmatism imaging*") == 1)
			localizationMethod = LOCALIZATION_ELLIPSGAUSS_ASTIG
		else
			Printf "The file at %s does not appear to be a valid positions file, skipping...\r", currentMacintoshFilePath
			Close refNum
			continue
		endif
		
		waveNote = "LOCALIZATION METHOD:" + num2str(localizationMethod) + ";"
		
		// parse the waveNote from the file, but stop when "DATA FOLLOWS" is encountered
		for (j = 1; ; j+=1)
			FReadLine refNum, singleLine
			if (strlen(singleLine) == 1)	// empty line
				continue
			endif
			if (stringmatch(singleLine, "DATA FOLLOWS*") == 1)
				firstLineWithData = j + 2
				break
			endif
			waveNote += singleLine + ";"
		endfor
		
		Close refNum
		
		// get rid of possible extra carriage returns produced by Igor
		waveNote = ReplaceString("\r", waveNote, "")
		waveNote += "POSITIONS FILE PATH:" + currentMacintoshFilePath + ";"
		
		SetDataFolder root:'Localized Positions'
		LoadWave /Q /A=$waveOutputName /G/ K=1 /M /L={firstLineWithData, firstLineWithData, 0, 0, 0} currentMacintoshFilePath
		SetDataFolder savDF
		
		if (ItemsInList(S_waveNames) != 1)
			Abort "An error occurred loading the data at " + currentMacintoshFilePath
		endif
		
		wave loadedWave = root:'Localized Positions':$StringFromList(0, S_waveNames)
		Note /K loadedWave, waveNote
	endfor
End

Function /WAVE GetPositionsWithinLimits(positions, xStart, xEnd, yStart, yEnd)
	wave positions
	variable xStart, xEnd, yStart, yEnd
	
	// xStart, xEnd, yStart, yEnd are all specified in pixels
	// returns a free wave
	
	variable nPos = DimSize(positions, 0)
	variable xCol, yCol, zCol
	variable nPositionsWithinLimits = 0, i, offset
	
	// get the indices of the columns containing the x and y coordinates
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1))
		Abort "The positions passed to ExtractSubsetOfPositions_menu() do not appear to contain any (x,y) information"
	endif
	
	
	for (i = 0; i < nPos; i+= 1)
		if ((positions[i][xCol] >= xStart) && (positions[i][yCol] >= yStart) && (positions[i][xCol] <= xEnd) && (positions[i][yCol] <= yEnd))
			nPositionsWithinLimits += 1
		endif
	endfor
	
	Make /D/FREE/N=(nPositionsWithinLimits, DimSize(positions, 1)) M_SelectedPositions
	
	offset = 0
	for (i = 0; i < nPos; i+= 1)
		if ((positions[i][xCol] >= xStart) && (positions[i][yCol] >= yStart) && (positions[i][xCol] <= xEnd) && (positions[i][yCol] <= yEnd))
			M_SelectedPositions[offset][] = positions[i][q]
			offset += 1
		endif
	endfor
	
	Note /K M_SelectedPositions, note(positions)
	
	return M_SelectedPositions
End

Function /WAVE GetTracksWithinLimits(tracksWave, modeStr, xStart, xEnd, yStart, yEnd)
	wave /WAVE tracksWave
	string modeStr	// "allinbox" or "centroidinbox"
	variable xStart, xEnd, yStart, yEnd
	
	variable nTracks = DimSize(tracksWave, 0)
	Make /FREE/N=(nTracks) /WAVE tracksWithinLimits
	
	variable allInBox = 0
	StrSwitch (modeStr)
		case "allinbox":
			allInBox = 1
			break
		case "centroidinbox":
			allInBox = 0
			break
		default:
			Abort "Unknown mode passed to GetTracksWithinLimits()"
			break
	EndSwitch
	
	variable i, xCol, yCol, zCol, nPosInTrack
	variable offset = 0
	for (i = 0; i < nTracks; i+=1)
		wave thisTrack = tracksWave[i]
		nPosInTrack = DimSize(thisTrack, 0)
		GetColumnsForEmitterPositions(thisTrack, xCol, yCol, zCol)
		if (xCol == -1)
			Abort
		endif
		
		MatrixOP /FREE W_xPos = col(thisTrack, xCol)
		MatrixOP /FREE W_yPos = col(thisTrack, yCol)
		if (allInBox)
			if ((WaveMin(W_xPos) < xStart) || (WaveMin(W_yPos) < yStart) || (WaveMax(W_xPos) > xEnd) || (WaveMax(W_yPos) > yEnd))
				continue
			endif
		else
			variable meanX = Mean(W_xPos)
			variable meanY = Mean(W_yPos)
			if ((meanX < xStart) || (meanY < yStart) || (meanX > xEnd) || (meanY > yEnd))
				continue
			endif
		endif
		
		tracksWithinLimits[offset] = thisTrack
		offset += 1
	endfor
	
	Redimension /N=(offset) tracksWithinLimits
	Note /K tracksWithinLimits, note(tracksWave)
	return tracksWithinLimits
End

Function IsPositionsWave(w)
	wave w
	
	// Require that positions waves have a "LOCALIZATION METHOD" entry in the wave note
	string waveNote = Note(w)
	if (StringMatch(waveNote, "*LOCALIZATION METHOD*"))
		return 1
	else
		return 0
	endif
End

Function IsTrackingWave(w)
	wave w
	
	// Require that tracking waves have "TRACKING MAX SHIFT" entry in the wave note
	string waveNote = Note(w)
	if (StringMatch(waveNote, "*TRACKING MAX SHIFT*"))
		return 1
	else
		return 0
	endif
End

Function IsRegistrationMap(w)
	wave w
	
	// Require that tracking waves have "TYPE:RegistrationMap" entry in the wave note
	string waveNote = Note(w)
	if (StringMatch(waveNote, "*TYPE:RegistrationMap*"))
		return 1
	else
		return 0
	endif
End

Function IsAstigmatism3DCalibrationMap(w)
	wave w
	
	// Require that tracking waves have "TYPE:RegistrationMap" entry in the wave note
	string waveNote = Note(w)
	if (StringMatch(waveNote, "*TYPE:Astigmatic3DCalibrationMap*"))
		return 1
	else
		return 0
	endif
End

Function IsSOFICombinationWeightsWave(w)
	wave w
	
	string waveNote = Note(w)
	if (StringMatch(waveNote, "*TYPE:SOFICombinationWeights*"))
		return 1
	else
		return 0
	endif
End

ThreadSafe Function /WAVE ExtractPositionsInFrame(PositionsWave, frameNumber, [positionsHaveBeenConsolidated])
	wave PositionsWave
	variable frameNumber
	variable positionsHaveBeenConsolidated
	
	// given a list of fitted positions, extract all positions that are present in frame frameNumber
	// do not return just the positions wave, but rather return a wave containing 2 wave references
	// the first reference is the positions wave, the second wave contains the indices of the found positions in the
	// original positions wave. If positionsHaveBeenConsolidated is provided then this routine does not attempt to
	// determine this parameter.
	
	variable checkForConsolidation
	if (ParamIsDefault(positionsHaveBeenConsolidated))
		checkForConsolidation = 1
	else
		checkForConsolidation = 0
	endif
	
	variable nPos = DimSize(PositionsWave, 0)
	variable start = 0, stop = nPos - 1, middle
	variable i
	variable nExtractedPos
	
	Make /O/N=2 /WAVE /FREE M_ExtractedPositions
	
	variable nFramesPresentCol
	GetColumnForNFramesPresent(PositionsWave, nFramesPresentCol)
	
	if (checkForConsolidation)
		// check if the positions have been consolidated before
		// if they haven't then use a binary search, if they have then perform a full loop
		ImageTransform /G=(nFramesPresentCol) sumCol, PositionsWave
		positionsHaveBeenConsolidated = V_value != DimSize(PositionsWave, 0)
	endif
	
	if (positionsHaveBeenConsolidated == 0)
		ImageTransform /G=0 getCol, PositionsWave
		wave W_ExtractedCol
		middle = BinarySearch(W_ExtractedCol, frameNumber)
		if ((middle < 0) || (PositionsWave[middle][0] != frameNumber))	// not found
			Make /D/FREE/O/N=(0, DimSize(PositionsWave, 1)) M_Positions
			Make /I/U/FREE/O/N=0 W_PositionIndices
		else
			for (i = middle; (i >= 0) && (PositionsWave[i][0] == frameNumber); i -= 1)
			endfor
			start = i + 1
			for (i = middle; (i < nPos) && (PositionsWave[i][0] == frameNumber); i += 1)
			endfor
			stop = i - 1
			
			nExtractedPos = stop - start + 1
			
			Make /D/FREE/O/N=(nExtractedPos, DimSize(PositionsWave, 1)) M_Positions
			Make /FREE/O/N=(nExtractedPos)/I/U W_PositionIndices
			
			M_Positions = PositionsWave[p + start][q]
			W_PositionIndices = p + start
		endif
	else		// positions have been consolidated
			// run over all frames up to the current frame to see if any of those emitters is still active
		Make /D/FREE/O/N=(0, DimSize(PositionsWave, 1)) M_Positions
		Make /I/U/FREE/O/N=0 W_PositionIndices
		variable offset = 0
		for (i = 0; (i < nPos) && (PositionsWave[i][0] <= frameNumber); i+=1)
			if ((PositionsWave[i][0] <= frameNumber) && ((PositionsWave[i][0] + PositionsWave[i][nFramesPresentCol]) > frameNumber))
				Redimension /N=(DimSize(M_Positions, 0) + 1, DimSize(M_Positions, 1)) M_Positions
				Redimension /N=(DimSize(W_PositionIndices, 0) + 1) W_PositionIndices
				M_Positions[offset][] = PositionsWave[i][q]
				W_PositionIndices[offset] = i
				offset += 1
			endif
		endfor
	endif
	
	M_ExtractedPositions[0] = M_Positions
	M_ExtractedPositions[1] = W_PositionIndices
	Note /NOCR /K M_ExtractedPositions[0], note(PositionsWave)
	
	KillWaves /Z W_ExtractedCol
	
	return M_ExtractedPositions
End

Function /WAVE ExtractTracksInFrame(tracksWave, frameNumber)
	wave /WAVE tracksWave
	variable frameNumber
	
	variable nTracks = DimSize(tracksWave, 0)
	Make /N=(nTracks) /WAVE /FREE W_ExtractedTracks
	Make /N=(nTracks) /D /FREE W_ExtractedTrackIndices
	variable nPositionsInTrack, firstOccurrence, lastOccurrence
	variable i, offset = 0
	for (i = 0; i < nTracks; i+=1)
		wave currentTrack = tracksWave[i]
		nPositionsInTrack = DimSize(currentTrack, 0)
		firstOccurrence = currentTrack[0][0]
		lastOccurrence = currentTrack[nPositionsInTrack - 1][0]
		
		if (firstOccurrence > frameNumber)
			break
		endif
		
		if ((firstOccurrence <= frameNumber) && (lastOccurrence >= frameNumber))
			Duplicate /FREE currentTrack, currentTrackCopy
			W_ExtractedTracks[offset] = currentTrackCopy
			W_ExtractedTrackIndices[offset] = i
			offset += 1
		endif
	endfor
	
	Redimension /N=(offset) W_ExtractedTracks, W_ExtractedTrackIndices
	Make /N=2 /WAVE /FREE W_ExtractedResults = {W_ExtractedTracks, W_ExtractedTrackIndices}
	return W_ExtractedResults
End

Function CountPositionsInTracks(tracksWave)
	wave /wave tracksWave
	
	variable nTracks = DimSize(tracksWave, 0)
	variable i, nPositions = 0
	for (i = 0; i < nTracks; i+=1)
		wave thisTrack = tracksWave[i]
		nPositions += DimSize(thisTrack, 0)
	endfor
	
	return nPositions
End

Function /WAVE FitMLEGaussToPositions(M_Positions)
	wave M_Positions
	
	variable nCoordinates = DimSize(M_Positions, 0)
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(M_Positions, xCol, yCol, zCol)
	
	MatrixOP /FREE W_xCol = col(M_Positions, xCol)
	WaveStats /Q W_xCol
	variable xx = V_avg
	variable sigmaX = V_sDev
	MatrixOP /FREE W_yCol = col(M_Positions, yCol)
	WaveStats /Q W_yCol
	variable yy = V_avg
	variable sigmaY = V_sDev
	variable corr = StatsCorrelation(W_xCol, W_yCol)
	
	Make /N=5/D/FREE W_GaussParams = {xx, yy, sigmaX, sigmaY, corr}
	
	return W_GaussParams
End

Function MLEGaussFitFunc(xx, yy, W_GaussParams)
	variable xx, yy
	wave W_GaussParams
	
	variable x0 = W_GaussParams[0]
	variable y0 = W_GaussParams[1]
	variable sX = W_GaussParams[2]
	variable sY = W_GaussParams[3]
	variable corr = W_GaussParams[4]
	
	variable amplitude = 1 / (2 * pi * sX * sY * (1 - corr^2))
	variable expScale = 1 / (2 * (1 - corr^2))
	variable result = (amplitude * exp(- expScale * ((xx - x0)^2 / sX^2 + (yy - y0)^2 / sY^2 - 2 * corr * (xx - x0) * (yy - y0) / (sX * sY))))
	
	return result
End

Function /WAVE SimulateRandomClusteringBounds(doLFunction, generateRandomPositions, positions, nSimulations, nPointsToSample, calculationRange, nBins, xSize, ySize, [positions2])
	variable doLFunction				// non-zero for L function, otherwise pairwise correlation
	variable generateRandomPositions	// non-zero for random points, otherwise randomly mixes positions from positions and positions2. Ignored for non-bivariate.
	wave positions	// only used as a template, not modified by this function
	variable  nSimulations // number of simulations to peform
	variable nPointsToSample	// number of points to sample from each positions wave
	variable calculationRange, nBins // passed to the /RNGE flag of the operation
	variable xSize, ySize // the dimensions of the region to simulate over
	wave /Z positions2	// if present then do bivariate simulations
	
	// the percentiles will be returned as a wave containing two wave references
	// elem 0 is the low percentile, elem 1 the high
	
	variable isBivariate = (!ParamIsDefault(positions2) && WaveExists(positions2))
	variable nPoints1 = DimSize(positions, 0)
	variable nPoints2 = DimSize(positions2, 0)
	variable nPointsSampledFrom1 = min(((nPointsToSample <= 0) ? nPoints1 : nPointsToSample), nPoints1)
	variable nPointsSampledFrom2 = min(round(nPointsSampledFrom1 * nPoints2 / nPoints1), nPoints2)
	
	DFREF tempDF = NewFreeDataFolder()
	DFREF savDF = GetDataFolderDFR()
	
	variable xCol1, yCol1, zCol1, xCol2, yCol2, zCol2
	GetColumnsForEmitterPositions(positions, xCol1, yCol1, zCol1)
	if (isBivariate)
		GetColumnsForEmitterPositions(positions2, xCol2, yCol2, zCol2)
	endif
	
	Duplicate /FREE positions, M_RandomSample
	Redimension /N=(nPointsSampledFrom1, -1), M_RandomSample
	if (isBivariate)
		Duplicate /FREE positions2, M_RandomSample2
		Redimension /N=(nPointsSampledFrom2, -1), M_RandomSample2
	endif
	if (isBivariate && !generateRandomPositions)
		Make /FREE /D/N=(nPoints1 + nPoints2) W_Random, W_ShuffledIndices
	endif
	Make /N=(nSimulations, nBins) /D /FREE M_SimulationResults
	Make /O/N=(nBins) /FREE/D W_LowPercentile, W_HighPercentile
	Make /N=2 /FREE /WAVE W_CombinedResults
	
	variable i, posIndex
	for (i = 0; i < nSimulations; i+=1)
		LocalizerProgressWindow("MonteCarlo", i, nSimulations)
		if (generateRandomPositions || !isBivariate)
			// RandomSample contains fully random coordinates
			M_RandomSample[][xCol1] = enoise(xSize / 2)
			M_RandomSample[][yCol1] = enoise(ySize / 2)
			if (isBivariate)
				M_RandomSample2[][xCol2] = enoise(xSize / 2)
				M_RandomSample2[][yCol2] = enoise(ySize / 2)
			endif
		else
			// RandomSample contains randomly-reassigned positions from positions and positions2
			W_Random = gnoise(1)
			WaveTransform index W_ShuffledIndices
			Sort W_Random, W_ShuffledIndices
			variable offset = 0, posIndexToTake
			for (posIndex = 0; posIndex < nPointsSampledFrom1; posIndex += 1)
				posIndexToTake = W_ShuffledIndices[offset]
				if (posIndexToTake < nPoints1)
					M_RandomSample[posIndex][] = positions[posIndexToTake][q]
				else
					M_RandomSample[posIndex][] = positions2[posIndexToTake - nPoints1][q]
				endif
				offset += 1
			endfor
			for (posIndex = 0; posIndex < nPointsSampledFrom2; posIndex += 1)
				posIndexToTake = W_ShuffledIndices[offset]
				if (posIndexToTake < nPoints1)
					M_RandomSample2[posIndex][] = positions[posIndexToTake][q]
				else
					M_RandomSample2[posIndex][] = positions2[posIndexToTake - nPoints1][q]
				endif
				offset += 1
			endfor
		endif
		
		SetDataFolder tempDF
		if (doLFunction)
			if (!isBivariate)
				RipleyLFunctionClustering /RNGE={calculationRange, nBins} M_RandomSample
			else
				RipleyLFunctionClustering /RNGE={calculationRange, nBins} M_RandomSample, M_RandomSample2
			endif
			wave W_Result = W_LFunction
			W_Result -= x
		else
			if (!isBivariate)
				PairwiseCorrelationClustering /RNGE={calculationRange, nBins} M_RandomSample
			else
				PairwiseCorrelationClustering /RNGE={calculationRange, nBins} M_RandomSample, M_RandomSample2
			endif
			wave W_Result = W_PairwiseCorrelation
		endif
		SetDataFolder savDF
		M_SimulationResults[i][] = W_Result[q]
	endfor
	
	SetScale /P x, DimOffset(W_Result, 0), DimDelta(W_Result, 0), W_LowPercentile, W_HighPercentile
	
	for (i = 0; i < nBins; i+=1)
		// get the percentiles of each bin
		MatrixOP /FREE W_ExtractedCol = col(M_SimulationResults, i)
		Sort W_ExtractedCol, W_ExtractedCol
		W_LowPercentile[i] = W_ExtractedCol(0.025 * (nSimulations - 1) * DimDelta(W_ExtractedCol, 0) + DimOffset(W_ExtractedCol, 0))
		W_HighPercentile[i] = W_ExtractedCol(0.975 * (nSimulations - 1) * DimDelta(W_ExtractedCol, 0) + DimOffset(W_ExtractedCol, 0))
	endfor
	
	W_CombinedResults[0] = W_LowPercentile
	W_CombinedResults[1] = W_HighPercentile
	
	return W_CombinedResults
End

Function /WAVE CalculateMSDs(tracksWave, maxTimeLag, minTrackLength, maxTrackLength)
	wave /wave tracksWave
	variable maxTimeLag, minTrackLength, maxTrackLength
	
	// returns 2 column wave - first column is MSD, second column is error on the MSD
	
	variable nTracks = DimSize(tracksWave, 0)
	variable nLags = maxTimeLag + 1	// +1 because we also include a point at time lag zero, which is always equal to zero
	
	Make /FREE/N=(nLags) /D W_MSD = 0, W_MSDError = 0, W_nEstimates = 0, W_nContributingTracks = 0, W_thisTrackContributions = 0
	Make /FREE/N=(nLags) /D W_Delta = 0, W_Mean = 0, W_M2 = 0	// temp variables for calculating online variance (http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm)
	
	SetScale /P x, 0, 1, W_MSD, W_MSDError, W_nEstimates, W_nContributingTracks
	
	variable i, lagTime, nPosInTrack, xCol, yCol, zCol
	variable first, second, timeLagBetweenThesePos, squaredDistance
	for (i = 0; i < nTracks; i+=1)
		wave thisTrack = tracksWave[i]
		W_thisTrackContributions = 0
		nPosInTrack = DimSize(thisTrack, 0)
		getColumnsForEmitterPositions(thisTrack, xCol, yCol, zCol)
		
		if ((minTrackLength > 0) && (nPosInTrack < minTrackLength))
			continue
		endif
		if ((maxTrackLength > 0) && (nPosInTrack > maxTrackLength))
			continue
		endif
		
		for (first = 0; first < nPosInTrack - 1; first += 1)
			for (second = first + 1; second < nPosInTrack; second += 1)
				timeLagBetweenThesePos = thisTrack[second][0] - thisTrack[first][0]
				if (timeLagBetweenThesePos > maxTimeLag)
					break
				endif
				
				squaredDistance = (thisTrack[second][xCol] - thisTrack[first][xCol])^2 + (thisTrack[second][yCol] - thisTrack[first][yCol])^2
				W_MSD[timeLagBetweenThesePos] += squaredDistance
				W_nEstimates[timeLagBetweenThesePos] += 1
				W_thisTrackContributions[timeLagBetweenThesePos] = 1
				
				// running variance calculation
				W_Delta[timeLagBetweenThesePos] = squaredDistance - W_Mean[timeLagBetweenThesePos]
				W_Mean[timeLagBetweenThesePos] += W_Delta[timeLagBetweenThesePos] / W_nEstimates[timeLagBetweenThesePos]
				W_M2[timeLagBetweenThesePos] += W_Delta[timeLagBetweenThesePos] * (squaredDistance - W_Mean[timeLagBetweenThesePos])
			endfor
		endfor
		W_nContributingTracks += W_thisTrackContributions
	endfor
	
	W_MSD = W_MSD[p] / W_nEstimates[p]
	W_MSDError = sqrt(W_M2[p] / (W_nEstimates[p] - 1)) / sqrt(W_nEstimates[p])
	MatrixOP /FREE W_MSDError = ReplaceNaNs(W_MSDError, 0)
	W_MSD[0] = 0
	W_MSDError[0] = 0
	W_nContributingTracks[0] = NaN
	
	Make /FREE /N=(nLags, 3) /D M_MSDResult
	SetScale /P x, 0, 1, M_MSDResult
	M_MSDResult[][0] = W_MSD[p]
	M_MSDResult[][1] = W_MSDError[p]
	M_MSDResult[][2] = W_nContributingTracks[p]
	
	return M_MSDResult
End

Function /WAVE EstimatePerTrackDisplacements(tracksWave)
	wave /WAVE tracksWave
	
	variable nTracks = DimSize(tracksWave, 0)
	Make /FREE/N=(nTracks)/D W_SqDisplacement
	
	variable trackIndex, xCol, yCol, zCol, nPosInTrack, posIndex, sqDisplacement, stepDuration
	for (trackIndex = 0; trackIndex < nTracks; trackIndex += 1)
		wave thisTrack = tracksWave[trackIndex]
		nPosInTrack = DimSize(thisTrack, 0)
		if (nPosInTrack < 2)
			W_SqDisplacement[trackIndex] = NaN
		endif
		getColumnsForEmitterPositions(thisTrack, xCol, yCol, zCol)
		
		sqDisplacement = 0
		for (posIndex = 1; posIndex < nPosInTrack; posIndex += 1)
			stepDuration = thisTrack[posIndex][0] - thisTrack[posIndex - 1][0]
			sqDisplacement += ((thisTrack[posIndex - 1][xCol] - thisTrack[posIndex][xCol])^2 + (thisTrack[posIndex - 1][yCol] - thisTrack[posIndex][yCol])^2) / stepDuration
		endfor
		sqDisplacement /= nPosInTrack - 1
		W_SqDisplacement[trackIndex] = sqDisplacement
	endfor
	
	WaveTransform /O sqrt W_SqDisplacement
	
	return W_SqDisplacement
End

Function GetMinAndMaxDisplacements(tracksWave, minDisplacement, maxDisplacement)
	wave /WAVE tracksWave
	variable &minDisplacement, &maxDisplacement
	
	minDisplacement = inf
	maxDisplacement = -1
	variable nTracks = DimSize(tracksWave, 0)
	
	variable trackIndex, xCol, yCol, zCol, nPosInTrack, posIndex, sqDisplacement, stepDuration
	for (trackIndex = 0; trackIndex < nTracks; trackIndex += 1)
		wave thisTrack = tracksWave[trackIndex]
		nPosInTrack = DimSize(thisTrack, 0)
		if (nPosInTrack < 2)
			continue
		endif
		getColumnsForEmitterPositions(thisTrack, xCol, yCol, zCol)
		
		for (posIndex = 1; posIndex < nPosInTrack; posIndex += 1)
			stepDuration = thisTrack[posIndex][0] - thisTrack[posIndex - 1][0]
			sqDisplacement = ((thisTrack[posIndex - 1][xCol] - thisTrack[posIndex][xCol])^2 + (thisTrack[posIndex - 1][yCol] - thisTrack[posIndex][yCol])^2) / stepDuration
			minDisplacement = Min(minDisplacement, sqDisplacement)
			maxDisplacement = Max(maxDisplacement, sqDisplacement)
		endfor
	endfor
	
	minDisplacement = sqrt(minDisplacement)
	maxDisplacement = sqrt(maxDisplacement)
End

Function /WAVE CalculateCDFs(tracksWave, maxTimeLag, maxSqDisplacement, nCDFPoints, minTrackLength, maxTrackLength)
	wave /wave tracksWave
	variable maxTimeLag, maxSqDisplacement, nCDFPoints, minTrackLength, maxTrackLength
	
	// return is 2D wave: rows contain different lag times and columns contain CDFs for each lag time
	
	variable nTracks = DimSize(tracksWave, 0)
	variable nLags = maxTimeLag
	variable sqDisplacementDelta = maxSqDisplacement / (nCDFPoints - 1)
	
	Make /FREE/N=(nLags, nCDFPoints) /D M_CDF = 0
	variable nDisplacements = 0
	SetScale /P x, 1, 1, M_CDF
	SetScale /P y, 0, sqDisplacementDelta, M_CDF
	Make /FREE /N=(nLags)/D W_nEstimates = 0
	
	variable i, j, allowableDisplacement, lagTime, nPosInTrack, xCol, yCol, zCol
	variable first, second, timeLagBetweenThesePos, squaredDistance
	for (i = 0; i < nTracks; i+=1)
		wave thisTrack = tracksWave[i]
		nPosInTrack = DimSize(thisTrack, 0)
		getColumnsForEmitterPositions(thisTrack, xCol, yCol, zCol)
		
		if ((minTrackLength > 0) && (nPosInTrack < minTrackLength))
			continue
		endif
		if ((maxTrackLength > 0) && (nPosInTrack > maxTrackLength))
			continue
		endif
		
		for (first = 0; first < nPosInTrack - 1; first += 1)
			for (second = first + 1; second < nPosInTrack; second += 1)
				timeLagBetweenThesePos = thisTrack[second][0] - thisTrack[first][0]
				if (timeLagBetweenThesePos > maxTimeLag)
					break
				endif
				
				W_nEstimates[timeLagBetweenThesePos - 1] += 1
				squaredDistance = (thisTrack[second][xCol] - thisTrack[first][xCol])^2 + (thisTrack[second][yCol] - thisTrack[first][yCol])^2
				if (squaredDistance <= maxSqDisplacement)
					M_CDF[timeLagBetweenThesePos - 1][ceil(squaredDistance / sqDisplacementDelta), ] += 1
				endif
			endfor
		endfor
	endfor
	
	M_CDF /= W_nEstimates[p]
	
	return M_CDF
End

Function /WAVE DisplacementHistograms(tracksWave, histogramRange, nHistogramBins, maxTimeLag, minTrackLength, maxTrackLength)
	wave /wave tracksWave
	variable histogramRange, nHistogramBins, maxTimeLag, minTrackLength, maxTrackLength
	
	variable nTracks = DimSize(tracksWave, 0)
	
	variable nLags = maxTimeLag
	variable binWidth = histogramRange / nHistogramBins
	variable sqHistogramRange = histogramRange^2
	Make /N=(nLags, nHistogramBins) /D/FREE M_DisplacementHistograms = 0
	SetScale /P x, 1, 1, M_DisplacementHistograms
	SetScale /P y, 0, binWidth, M_DisplacementHistograms
	
	variable i, j, lagTime, nPosInTrack, xCol, yCol, zCol
	variable first, second, timeLagBetweenThesePos, squaredDistance
	for (i = 0; i < nTracks; i+=1)
		wave thisTrack = tracksWave[i]
		nPosInTrack = DimSize(thisTrack, 0)
		getColumnsForEmitterPositions(thisTrack, xCol, yCol, zCol)
		
		if ((minTrackLength > 0) && (nPosInTrack < minTrackLength))
			continue
		endif
		if ((maxTrackLength > 0) && (nPosInTrack > maxTrackLength))
			continue
		endif
		
		for (first = 0; first < nPosInTrack - 1; first += 1)
			for (second = first + 1; second < nPosInTrack; second += 1)
				timeLagBetweenThesePos = thisTrack[second][0] - thisTrack[first][0]
				if (timeLagBetweenThesePos > maxTimeLag)
					break
				endif
				
				squaredDistance = (thisTrack[second][xCol] - thisTrack[first][xCol])^2 + (thisTrack[second][yCol] - thisTrack[first][yCol])^2
				if (squaredDistance < sqHistogramRange)
					M_DisplacementHistograms[timeLagBetweenThesePos - 1][floor(sqrt(squaredDistance) / binWidth)] += 1
				endif
			endfor
		endfor
	endfor
	
	return M_DisplacementHistograms
End






































































//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
// D. Dibble created functions-------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------

Function ExportScatterData(positions, xSize, ySize, pixelSize, FilteringType, imageFileLocation)
	wave positions
	variable xSize, ySize, pixelSize, FilteringType
	string imageFileLocation
	
	Assert(WaveExists(positions))
	Assert((xSize > 0) && (ySize > 0))
	
	variable nPositions, nImages, pixelSizeUM = pixelSize / 1000
	

	// create blank lists to store the X, Y, and "color" aka position the event occurs at

	Make /D/O Scatter_X, Scatter_Y,  Filter_X, Filter_Y
	Make /I/O Scatter_Color, Filter_Color
	




	// reload the image for manipulation outside of the localizer XOP 

	// if you can't reload don't process anything else
	if(numtype(strlen(imageFileLocation)) != 2)

		//-------------------------------------------------------------------------
		// Open the image file for processing
		//-------------------------------------------------------------------------
		ImageLoad/T=tiff/S=0/C=-1/Q imageFileLocation
		String/G ChannelOneReference = StringFromList(0,s_waveNames)

		


		//-----------------------------------------------------------------------------------------------------------
		// for the non-track data, pull out the points, for track data pull the points out of the individual tracks
		//-----------------------------------------------------------------------------------------------------------
		if(FilteringType != 2)

			nPositions = DimSize(positions, 0)
			nImages = positions[nPositions - 1][0] - positions[0][0] + 1

			redimension/N=(nPositions)  Scatter_X, Scatter_Y, Scatter_Color

			// get the indices of the columns containing the x and y coordinates
			variable xCol, yCol, zCol
			getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
			if ((xCol == -1) || (yCol == -1))
				Abort "The positions passed to ExportScatterData() do not appear to contain any (x,y) information"
			endif

			// get the data from main localizer utility
			Scatter_X = positions[p][xCol]
			Scatter_Y = positions[p][yCol]
			Scatter_Color = positions[p][0]
		elseif(FilteringType == 2)
			ExtractTracksToScatter(positions, Scatter_X, Scatter_Y)
		endif	
		//---------------------------------------------------------------------------------------------------------	
		

		//--------------------------------------------------------------------------------------------------------------------
		// filter the data for duplicate points and points on the edge of the image is that filtering is requested
		// print a status update on the command line
		// if tracks are requested, still filter edge points
		//--------------------------------------------------------------------------------------------------------------------

		if(FilteringType == 1)

			Print "There are: " + num2str(nPositions) + " unfiltered points localized" 
			ScatterDataFilterImageEdge(Scatter_X, Scatter_Y, Scatter_Color)
			Print "There are: " + num2str(DimSize(Scatter_X, 0)) + " unfiltered points localized that are not on the edge of the image"
			ScatterDataFilter(Scatter_X, Scatter_Y, Scatter_Color, Filter_X, Filter_Y, Filter_Color)	
			Print "There are: " + num2str(DimSize(Filter_X, 0)) + " filtered points localized (removing duplicated integer coordinates)" 

		elseif(FilteringType == 2)

			// create a color wave (doesn't exist for tracks), needed for interface
			Redimension/N=(numpnts(Scatter_X)) Scatter_Color 
			Scatter_color = 1

			Print "There are: " + num2str(Numpnts(Scatter_X)) + " unfiltered points localized" 
			ScatterDataFilterImageEdge(Scatter_X, Scatter_Y, Scatter_Color)
			Print "There are: " + num2str(DimSize(Scatter_X, 0)) + " unfiltered points localized that are not on the edge of the image"
			ScatterDataFilter(Scatter_X, Scatter_Y, Scatter_Color, Filter_X, Filter_Y, Filter_Color)	
			Print "There are: " + num2str(DimSize(Filter_X, 0)) + " filtered points localized (removing duplicated integer coordinates)"


			// duplicate Scatter into filter waves, needed for interface
			//Redimension/N=(numpnts(Scatter_X)) Filter_X, Filter_Y 
			//Filter_X = floor(Scatter_X)
			//Filter_Y = floor(Scatter_Y)

		endif
		//--------------------------------------------------------------------------------------------------------------------	
		

		//----------------------------------------------------------------------------------------------------------------		
		// set up a place to store all of the trace data and/or clean it up if it already exists (maybe not)
		//----------------------------------------------------------------------------------------------------------------		
		if(DataFolderExists("root:SingleParticleTraceData"))
			SetDataFolder root:SingleParticleTraceData
			KillWaves/A/Z
			SetDataFolder root:
		else
			NewDataFolder root:SingleParticleTraceData
			SetDataFolder root:
		endif

		// set up a place to store all of the curve fit data and/or clean it up if it already exists (maybe not)

		if(DataFolderExists("root:SingleParticleTraceData:Statistics"))
			SetDataFolder root:SingleParticleTraceData:Statistics
			KillWaves/A/Z
			SetDataFolder root:
		else
			NewDataFolder root:SingleParticleTraceData:Statistics
			SetDataFolder root:
		endif



	//----------------------------------------------------------------------------
	// create waves associated with the statistics, empty for now to be added to 
	// when evauluated  
	//----------------------------------------------------------------------------

	String savedDataFolder = GetDataFolder(1)	// Save
	SetDataFolder root:SingleParticleTraceData:Statistics

		String/G StatWaveNames = "Stat_X;Stat_Y;Stat_CurveFitBegin;Stat_CurveFitEnd;Stat_EventStart;Stat_EventEnd;Stat_EventHeight;Stat_MeanHeight"

		// make as text waves to directly output to file when done
		variable/G StatWaveSize = 1
		Make /O/T/N=(StatWaveSize) Stat_X
		Make /O/T/N=(StatWaveSize) Stat_Y
		Make /O/T/N=(StatWaveSize) Stat_CurveFitBegin
		Make /O/T/N=(StatWaveSize) Stat_CurveFitEnd
		Make /O/T/N=(StatWaveSize) Stat_EventStart
		Make /O/T/N=(StatWaveSize) Stat_EventEnd
		Make /O/T/N=(StatWaveSize) Stat_EventHeight
		Make /O/T/N=(StatWaveSize) Stat_MeanHeight

		// raw number data for later conversion or processing stored in a matrix
		variable ParametersToStore = 8
		Make /O/N=(ParametersToStore) Stat_TempMatrix
		// implement this if we start doing data processing in app
		//Make /O/N=(StatWaveSize,ParametersToStore) Stat_RawMatrix

		

		// put in the header for each
		Stat_X[0] = "X"
		Stat_Y[0] = "Y"
		Stat_CurveFitBegin[0] = "Curve Fit Start"
		Stat_CurveFitEnd[0] = "Curve Fit End"
		Stat_EventStart[0] = "Event Start"
		Stat_EventEnd[0] = "Event End"
		Stat_EventHeight[0] = "Event Height"
		Stat_MeanHeight[0] = "Event Mean Height"

	 

	SetDataFolder savedDataFolder




		// set up an arrays to contain the names of all the data waves and their colors for later viewing in a plot
		Make /T/O/N=(DimSize(Filter_X, 0)) WaveNameArray
		// to contain names of baseline-adjusted waves
		Make /T/O/N=(DimSize(Filter_X, 0)) BaselineWaveNameArray


		Make /T/O/N=(DimSize(Filter_X, 0)) WaveColorArray
		// set up arrays to contain the X and Y coordinates that corresponding to each data wave
		Make /O/N=(DimSize(Filter_X, 0)) WaveNameArray_X
		Make /O/N=(DimSize(Filter_X, 0)) WaveNameArray_Y

		// user decision to keep or throw away a wave, >=1 keep, <= -1 throw away, otherwise undecided
		Make /I/O/N=(DimSize(Filter_X, 0)) KeepWaveArray
		KeepWaveArray = 0
		// information on if a wave has been plotted or not. only plot if required for efficiencies sake
		Make /I/O/N=(DimSize(Filter_X, 0)) WavePlottedArray
		WavePlottedArray = 0


	

		// display the plots and interface for manipulation
		PlotWaveSelection()


		// clean up time
		// killWaves Scatter_X, Scatter_Y, Scatter_Color, Filter_X, Filter_Y, Filter_Color


	// can't load print error and quit	
	else
		Print "Error loading image file"
		
	endif

End







//---------------------------------------------------------------------------------------------------------
// pulls out all the x,y pairs out of the tracks and place them in seperate x,y waves for later editing
//---------------------------------------------------------------------------------------------------------
Function ExtractTracksToScatter(positions, Scatter_X, Scatter_Y)
	Wave /WAVE positions
	wave Scatter_X, Scatter_Y


	//variable nTracks = DimSize(positions, 0)
	variable nTracks = NumPnts(positions)
	variable NumberOfPoints = 0
	variable TrackItor, PointItor
	variable ScatterIndexStart, xCol, yCol, zCol
	
	// itorate over the tracks
	for (TrackItor = 0; TrackItor < nTracks; TrackItor++)
		
		wave CurrentPosition = positions[TrackItor]

		GetColumnsForEmitterPositions(CurrentPosition, xCol, yCol, zCol)
		
		// get the number of points and increase the size of the scatter waves 
		ScatterIndexStart = NumberOfPoints
		NumberOfPoints = NumberOfPoints + DimSize(CurrentPosition,0)

		// change wave sizes
		Redimension/N=(NumberOfPoints) Scatter_X, Scatter_Y  

		// itorate over the individual points in the tracks storing them
		for(PointItor = 0; PointItor < DimSize(CurrentPosition,0); PointItor++)
			Scatter_X[ScatterIndexStart + PointItor] = CurrentPosition[PointItor][xCol]
			Scatter_Y[ScatterIndexStart + PointItor] = CurrentPosition[PointItor][yCol]
		endfor

	endfor

End

//-------------------------------------------------------------------
// Function to remove points that are on the very edge of the image
// all points must be at least one pixel away from the edge for
// signal averaging
//-------------------------------------------------------------------


Function ScatterDataFilterImageEdge(Scatter_X, Scatter_Y, Scatter_Color)
	Wave Scatter_X, Scatter_Y, Scatter_Color
	

	//------------------------------------------------------------------------
	// get information about the image dimensions for removing border points
	//------------------------------------------------------------------------

	svar ChannelOneReference=root:ChannelOneReference
	// ChannelOneReference is in root directory, switch to get the reference for the wave then back to current directory
	String savedDataFolder = GetDataFolder(1)	// Save
	SetDataFolder root:	
	Wave ChannelOne=$ChannelOneReference
	SetDataFolder savedDataFolder
	// get information about the image
	variable image_width = dimsize(ChannelOne, 0)
	variable image_height = dimsize(ChannelOne, 1)

	variable itor, X_Point, Y_Point
	// itorate over all the points checking for edge cases
	for(itor = dimsize(Scatter_X,0)-1; itor >=0; itor--)
		//----------------------------------------------------------------------------------
		// use integers for the comparison, indexed at 0 so the floor value is appropriate
		//----------------------------------------------------------------------------------
		X_Point = floor(Scatter_X[itor])
		Y_Point = floor(Scatter_Y[itor])

		if((X_Point > 0)&&(Y_Point > 0)&&(X_Point < (image_width - 1))&&(Y_Point < (image_height - 1)))
			// do nothing, the point is within range
		else
			//delete the point from the list of points and frames 
			deletepoints/M=0 itor,1,Scatter_X
			deletepoints/M=0 itor,1,Scatter_Y
			deletepoints/M=0 itor,1,Scatter_Color
			

		endif
	endfor

End







// function filters the data for duplicated points 
Function ScatterDataFilter(Scatter_X, Scatter_Y, Scatter_Color, Filter_X, Filter_Y, Filter_Color)
	Wave Scatter_X, Scatter_Y, Scatter_Color, Filter_X, Filter_Y, Filter_Color
	

	variable itor_X, itor_Y, itor_List, itor_Check, filter_Size, itor
	int duplicated = 0


	// enter the first point into the filtered data
	Filter_Size = 1
	Redimension/N=(filter_Size) Filter_X
	Redimension/N=(filter_Size) Filter_Y
	Redimension/N=(filter_Size) Filter_Color

	Filter_X[0] =  Scatter_X[0]
	Filter_Y[0] =  Scatter_Y[0]
	Filter_Color[0] =  Scatter_Color[0]

	// create a temporary integer version of the X and Y data, floating point is tricky on the comparison and some duplicates are missed
	// use the floor value, pixels are indexed from 0 and you can't round up to a non-existent coordinate

	Make /I/O/N=(dimsize(Scatter_X,0)) Integer_X
	Make /I/O/N=(dimsize(Scatter_Y,0)) Integer_Y
	for(itor = 0; itor < dimsize(Scatter_X,0); itor++)
		Integer_X[itor] = floor(Scatter_X[itor])
		Integer_Y[itor] = floor(Scatter_Y[itor])
	endfor	


	// start by going through the unfiltered list item by item 
	for(itor_List = 1; itor_List < dimsize(Scatter_X,0); itor_List++)

	// check if the previous items do not match, if they don't
	// then include the new item in the growing filtered list

		for(itor_Check = 0; itor_Check < itor_List; itor_Check++)
			if((Integer_X[itor_List]==Integer_X[itor_Check])&&(Integer_Y[itor_List]==Integer_Y[itor_Check]))
				duplicated = 1				
			endif
		endfor

		// if not duplicated , add to the list, if it is reset the check switch

		if(duplicated == 0)
			filter_Size++
			Redimension/N=(filter_Size) Filter_X
			Redimension/N=(filter_Size) Filter_Y
			Redimension/N=(filter_Size) Filter_Color

			// old stored the integer values
			//Filter_X[filter_Size - 1] = Integer_X[itor_List]
			//Filter_Y[filter_Size - 1] = Integer_Y[itor_List]

			// new stores the floating point values
			Filter_X[filter_Size - 1] = Scatter_X[itor_List]
			Filter_Y[filter_Size - 1] = Scatter_Y[itor_List]


			Filter_Color[filter_Size - 1] =  Scatter_Color[itor_List]
		else
			duplicated = 0	
		endif	

	endfor

	// clean up
	KillWaves Integer_X
	killWaves Integer_Y

End


Function RetrieveIntensityPlot(Filter_X, Filter_Y, Filter_Frame, PointIndex)
	Wave Filter_X, Filter_Y, Filter_Frame
	variable PointIndex

	svar ChannelOneReference=root:ChannelOneReference
	// ChannelOneReference is in root directory, switch to get the reference for the wave then back to current directory
	String savedDataFolder = GetDataFolder(1)	// Save
	SetDataFolder root:	
	Wave ChannelOne=$ChannelOneReference
	SetDataFolder savedDataFolder


	// get information about the image
	variable image_width = dimsize(ChannelOne, 0)
	variable image_height = dimsize(ChannelOne, 1)


	// make sure that the point is at least 1 pixel from the edge before including, this may be duplicate of what filtering does above
	if((floor(Filter_X[PointIndex]) > 0)&&(floor(Filter_Y[PointIndex]) > 0)&&(floor(Filter_X[PointIndex]) < (image_width - 1))&&(floor(Filter_Y[PointIndex]) < (image_height - 1)))

		Make/O/N=(dimsize(ChannelOne, 2)) average
		average = 0
		Make/O/N=(dimsize(ChannelOne, 2), 3) RGB_Value


		variable itor_X, itor_Y 
		// sum all of the waves in a three by three grid
		for(itor_X = -1; itor_X < 2; itor_X++)
			for(itor_Y = -1; itor_Y < 2; itor_Y++)
				ImageTransform/beam={(floor(Filter_X[PointIndex]) + itor_X), (floor(Filter_Y[PointIndex]) + itor_Y)} getbeam ChannelOne
				wave W_Beam
				average = average + W_Beam
			endfor
		endfor
		average = average / 9.0		


		// fill all with zero be default
		if((Filter_Frame[PointIndex]-1) < -0.9)
				RGB_Value[0,0][0] = 0
		else		
				RGB_Value[0,Filter_Frame[PointIndex]-1][0] = 0
		endif		
		// fill all remaining with red based on the value of filter color
		RGB_Value[Filter_Frame[PointIndex],][0] = 65535
			
		RGB_Value[][1] = 0
		RGB_Value[][2] = 0

		StoreIndividualIntensityWave(Filter_X[PointIndex], Filter_Y[PointIndex], average, RGB_Value, PointIndex)
		average = 0


	endif


End


// stores the information about the wave *and* creates a baseline adjusted wave for later processing if needed
Function StoreIndividualIntensityWave(point_X, point_Y, average, RGB_Value, PointIndex)
	variable point_X, point_Y
	wave average, RGB_Value
	variable PointIndex	

	wave/T WaveList = root:WaveNameArray
	wave/T BaseWaveList = root:BaselineWaveNameArray
	wave WaveList_X = root:WaveNameArray_X
	wave WaveList_Y = root:WaveNameArray_Y
	wave/T WaveColorList = root:WaveColorArray
	wave/I KeepWaveList = root:KeepWaveArray
	wave/I WavePlottedList = root:WavePlottedArray
	
	//------------------------------
	// generate the new wave names
	//------------------------------
	String IntensityData_Name = "X" + num2str(floor(point_X)) + "Y" + num2str(floor(point_Y))
	String BaselineIntensityData_Name = "X" + num2str(floor(point_X)) + "Y" + num2str(floor(point_Y)) + "B"	
	String Frame_Data_Name = "X" + num2str(floor(point_X)) + "Y" + num2str(floor(point_Y)) + "C"
	//---------------------------------------------------------
	// store the new wave names and some basic info about them
	//---------------------------------------------------------
	WaveList[PointIndex] = IntensityData_Name
	BaseWaveList[PointIndex] = BaselineIntensityData_Name
	WaveList_X[PointIndex] = point_x
	WaveList_Y[PointIndex] = point_y
	WaveColorList[PointIndex] = Frame_Data_Name
	KeepWaveList[PointIndex] = 0
	WavePlottedList[PointIndex] = 1

	//-------------------------------------------------------------------------------------------------------
	// change data folder to the location where the recently generated waves are to be stored and save them
	// including generating a baseline adjusted wave
	//-------------------------------------------------------------------------------------------------------

	SetDataFolder root:SingleParticleTraceData

		duplicate average $IntensityData_Name
		duplicate RGB_Value $Frame_Data_Name

		//------------------------------------------------------------------------------------------------------------------------------------
		// this chunk generates a baseline adjusted wave and saves it
		// be sure to check out:
		// Friedrichs, M. S., A model-free algorithm for the removal of baseline artifacts. Journal of Biomolecular NMR 1995, 5 (2), 147–153.
		//------------------------------------------------------------------------------------------------------------------------------------
		// get the baseline 
		duplicate /O average $BaselineIntensityData_Name
		smooth /M=0 (floor(numPnts($BaselineIntensityData_Name)/3)), $BaselineIntensityData_Name
		//subtract it
		Wave BaseL = $BaselineIntensityData_Name
		Wave Intens = $IntensityData_Name	
		BaseL = Intens - BaseL

	SetDataFolder root:	


End



//-------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
// section associated with the plotting waves, selecting them, and curve fitting if needed
//-------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------

Function PlotWaveSelection()

	
	variable /G PointIndex = 0
	variable /G PointIndexOld = 0
	variable /G TotalPoints
	variable /G CurveFitRepeat = 10
	wave/T WaveNameList = root:WaveNameArray
	wave/T BaseWaveNameList = root:BaselineWaveNameArray
	wave/T WaveColorList = root:WaveColorArray
	// integer flag to indicate filter status
	wave/I KeepWaveList = root:KeepWaveArray
	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color

	TotalPoints = NumPnts(X_Wave) - 1

	// create the very first plot 
	RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, PointIndex)

	SetDataFolder root:SingleParticleTraceData
	Variable/G MinCurveFitFrame = 0
	Variable/G MaxCurveFitFrame = numpnts($WaveNameList[0]) - 1
	// variable defining the number of trials of simulated annealing to loop over
	Variable/G NumberTrialFits = 10
	// wave to store the curve defining the fit 
	Make /O/N=(numpnts($WaveNameList[0])) CurveFitWave
	CurveFitWave = 0


	//-----------------------------------------------------
	// create a panel to set up the interface in
	// as well as to listen to certain keystroke events 
	//-----------------------------------------------------
	NewPanel /W=(100,100,775,600)/N=PointChooserPanel
	SetDrawLayer UserBack

	// listen for arrow keystrokes to press as well as buttons
	SetWindow PointChooserPanel, hook(MyHook)=KeyboardWindowHook

	//----------------------------------------------------------------
	//set up the chart to include a data wave and a curve fit wave
	//----------------------------------------------------------------

	Display /W=(20,20,520,430)/N=PointChooser/HOST=# $WaveNameList[0] $BaseWaveNameList[0] CurveFitWave

	//----------------------------------------------------------------------------------------------------------------------------
	// section to set the initial range and information about the axis, ignoring the baseline wave and curve fit (not calculated)
	//----------------------------------------------------------------------------------------------------------------------------

	// zoom out by 10 percent on the vertical data to avoid hitting the top and bottom for ease of viewing 
	variable range = wavemax ($WaveNameList[PointIndex]) - wavemin ($WaveNameList[PointIndex])
	variable delta = range * 0.1 / 2.0	
	SetAxis /W=PointChooserPanel#PointChooser left,(wavemin($WaveNameList[PointIndex]) - delta),(wavemax($WaveNameList[PointIndex]) + delta)	
	// set mirror axis with a thickness of 3, tick inside, major tick length 7, major tick thick 3, no minor tick, no axis standoff, no labels on top and right axis 
	ModifyGraph /W=PointChooserPanel#PointChooser mirror=1,axThick=3,tick=2,btLen=7,btThick=3,minor=0,standoff=0
	//set font of axis at bold 18 
	ModifyGraph /W=PointChooserPanel#PointChooser fstyle=1,fsize=18,font="arial"
	Label bottom "Frame"
	Label left "Intensity"

	//----------------------------------------------------------------------------------------------------------
	// section to change the background color of the chart to indicate the selection state of the current wave
	//----------------------------------------------------------------------------------------------------------

	// paint background white (undecided), extra care here because we aren't dealing with integers
	if(abs(KeepWaveList[PointIndex]) < 0.5)
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,65535,65535)
	// paint background redish (throwin it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[0]) > 0.5)&&(KeepWaveList[0] < 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,40863,40863)
	// paint background greenish (keepin it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[PointIndex]) > 0.5)&&(KeepWaveList[0] > 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(40863,65535,40863)
	endif
	
	//--------------------------------------------------------------------
	//section to modify the traces including thickness, type, color, etc
	//--------------------------------------------------------------------

	// set the line thickness of the data wave,color the data wave to indicate the point that localizer chose acording to WaveColorList
	ModifyGraph /W=PointChooserPanel#PointChooser lSize($WaveNameList[PointIndex])=3,zColor($WaveNameList[PointIndex])={$WaveColorList[PointIndex],*,*,directRGB}
	// similar settings for the baseline wave, but make it all re
	ModifyGraph /W=PointChooserPanel#PointChooser lSize($BaseWaveNameList[PointIndex])=3,rgb($BaseWaveNameList[PointIndex])=(65535,0,0)
	// for the curve fit, set the initial line thickness (3), type (dashed), and color(blue)
	ModifyGraph /W=PointChooserPanel#PointChooser rgb(CurveFitWave)=(0,0,65535),lstyle(CurveFitWave)=3, lSize(CurveFitWave)=3

	//---------------------------------------------------------------------
	// section to navigate through the plots using buttons
	//---------------------------------------------------------------------
	Button Last,size={120,20},pos={540,10},title="Previous",proc=PreviousPlot
	Button Next,size={120,20},pos={540,40},title="Next",proc=NextPlot
	SetVariable CurrentTrace,pos={540.00,70.00},size={120.00,18.00},title="\\F'Arial'\\f01Current Trace",limits={0,(numpnts(X_Wave)) - 1,0},value=root:PointIndex,proc=ChangeCurrentTrace
	//CheckBox RandomPick,pos={547.00,100.00},size={106.00,16.00},title="Random Plot Pick",value= 0

	ValDisplay TotalTraces,pos={545.00,95.00},size={115.00,18.00},title="\\F'Arial'\\f01Total Traces:",limits={0,0,0},barmisc={0,1000},value=#"TotalPoints"

	//accepting or rejecting plots
	Button Accept,size={120,20},pos={540,160},title="Accept",proc=AcceptPlot
	Button Reject,size={120,20},pos={540,190},title="Reject",proc=RejectPlot 
	
	//exporting selected traces
	Button ExportAccepted,size={120,20},pos={540,240},title="Export Accepted",proc=ExportAccept
	Button ExportRejected,size={120,20},pos={540,270},title="Export Rejected",proc=ExportReject
	Button ExportUndecided,size={120,20},pos={540,300},title="Export Undecided",proc=ExportUndecided

	// exporting baseline adjusted plots and curve fit statistics
	Button ExportBaselineAdjusted,size={120,40},pos={540,360},title="Export\rBaseline-Adjusted",proc=ExportBaselined

	// export the curve fit data
	Button ExportCurveFitData,size={120,20},pos={540,410},title="Export Curve Fits",proc=ExportCurveFits

	// leaving the window and on to the next task
	Button ReturnMain,size={120,20},pos={540,470},title="Return",proc=ReturnMain

	//-----------------------------------------------------------
	// area to do the curve fitting of the selected trace
	//-----------------------------------------------------------
	GroupBox TraceFitting,pos={20.00,440.00},size={500.00,50.00},title="Curve Fitting"
	SetVariable SetMinX,pos={25,460.00},size={100.00,18.00},limits={0,(numpnts($WaveNameList[0])) - 1,0},value=MinCurveFitFrame,title="X Minimum",proc=CurveFitSetMinX
	SetVariable SetMaxX,pos={130,460.00},size={100.00,18.00},limits={0,(numpnts($WaveNameList[0])) - 1,0},value=MaxCurveFitFrame,title="X Maximum",proc=CurveFitSetMaxX
	SetVariable NumberRepeat,pos={234.00,461.00},size={150.00,19.00},title="\\F'Arial'\\f01Curve Fit Repeat:",limits={1,inf,1},value=CurveFitRepeat

	Button ApplyCurveFit,size={62,20},pos={388,461},title="Apply Fit",proc=CurveFitPerform,disable=2
	Button SaveCurveFit,size={62,20},pos={453,461},title="Accept Fit",proc=CurveFitInclude,disable=2

	SetDataFolder root:	

End


Function CurveFitSetMinX(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	nvar MinX = root:SingleParticleTraceData:MinCurveFitFrame
	nvar MaxX = root:SingleParticleTraceData:MaxCurveFitFrame

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update

			if(sva.dval < MaxX)
				MinX = sva.dval
			else
				MinX = 0
			endif

			break
	endswitch

	return 0
End

Function CurveFitSetMaxX(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	nvar MinX = root:SingleParticleTraceData:MinCurveFitFrame
	nvar MaxX = root:SingleParticleTraceData:MaxCurveFitFrame
	wave/T WaveNameList = root:WaveNameArray

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update

			if(sva.dval > MinX)
				MaxX = sva.dval
			else
				MaxX = numpnts($("root:SingleParticleTraceData:"+WaveNameList[0]))
			endif

			break
	endswitch

	return 0
End




//--------------------------------------------------------------------------------------------------------------------------------
// Functions to fit a square wave pulse to the raw data 
//--------------------------------------------------------------------------------------------------------------------------------
Function CurveFitPerform(CtrlName) : ButtonControl
	String CtrlName

	nvar MinimumX = root:SingleParticleTraceData:MinCurveFitFrame
	nvar MaximumX = root:SingleParticleTraceData:MaxCurveFitFrame
	nvar NumberRuns = root:CurveFitRepeat
	nvar Point = root:PointIndex
	Wave/T BaseWaveNameList = root:BaselineWaveNameArray
	Wave CurveFitWave = root:SingleParticleTraceData:CurveFitWave
	Wave BaseOriginalWave = $("root:SingleParticleTraceData:" + BaseWaveNameList[Point])
	Wave TempRawData = root:SingleParticleTraceData:Statistics:Stat_TempMatrix
	// needed to store as information about the wave
	nvar Point = root:PointIndex
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y




	// reset some values between runs
	RescaleBaselinePlot()
	Redimension /N=(numpnts(BaseOriginalWave)) CurveFitWave
	SetScale/I x,MinimumX,MaximumX,CurveFitWave
	
	//---------------------------------------------------------------------------------------------------------------------
	// optimization using simulated annealing operation
	//---------------------------------------------------------------------------------------------------------------------

	// the answer shows quite a bit of variation per run
	// lets try to run this function a user defined number of times, and select the best fit as the answer
	// store values in a multidimensional wave



	variable NumberParameters = 3
	//optimization input for the curve fit routine
	Make /O DummyWave
	Make /O/N=(3,2) ParameterRangeSA
	Make /O/N=3 ParameterSolution

	// optimization output
	Make /O/N=(NumberRuns,NumberParameters) RunOutPut 
	Make /O/N=(NumberRuns) RunMinOutput
	variable itor
	for(itor = 0; itor < NumberRuns; itor++)

		// set up the guesses for the fit, they are overwritten each run 
		// so redo and restart
		// first parameter is the location of the rise, should be between the the range set by the user
		ParameterRangeSA[0][0] = MinimumX
		ParameterRangeSA[0][1] = MaximumX
		// second parameter is the location of the fall, should be between the the range set by the user
		ParameterRangeSA[1][0] = MinimumX
		ParameterRangeSA[1][1] = MaximumX
		// Third parameter is the height of the peak relative too the baseline. must be a positive number
		// since we are talking about a baseline-adjusted peak here, no more than 10% of the max value seems a reasonable
		// bracket 
		ParameterRangeSA[2][0] = 0
		ParameterRangeSA[2][1] = Wavemax(BaseOriginalWave, MinimumX,MaximumX) + 0.1 * Wavemax(BaseOriginalWave, MinimumX,MaximumX)
	
		// guesses for the value of the squarewave
		// first guess is the bracket values 
		ParameterSolution[0] = MinimumX
		ParameterSolution[1] = MaximumX	
		// go with the maximum for a first guess when baseline adjusted
		ParameterSolution[2] =  Wavemax(BaseOriginalWave, MinimumX,MaximumX) 

		Optimize /M={3,0} /Q /TSA={0,0.5} /XSA=ParameterRangeSA /X=ParameterSolution BaselinePeakGeometryCalculation,DummyWave
		RunMinOutput[itor] = V_min	
		RunOutput[itor][0] = ParameterSolution[0]
		RunOutput[itor][1] = ParameterSolution[1]
		RunOutput[itor][2] = ParameterSolution[2]

		CurveFitWave = ParameterSolution[2] * (HeavisideApproximation(x,ParameterSolution[0]) - HeavisideApproximation(x,ParameterSolution[1]))
		ModifyGraph /W=PointChooserPanel#PointChooser lstyle(CurveFitWave)=0,rgb(CurveFitWave)=(0,0,0)
		doupdate /W=PointChooserPanel#PointChooser
		print "Fit Run: " + num2str(itor) + " Minimum: " + num2str(V_min)
	endfor	


	
	Button SaveCurveFit,disable=0

	WaveStats /Q RunMinOutput 
	variable BestFitIndex = V_minloc

	Print "Best Fit Run: "+num2str(BestFitIndex)

	variable PeakHeight, PeakStart, PeakEnd

	PeakStart = RunOutput[BestFitIndex][0]
	PeakEnd = RunOutput[BestFitIndex][1]
	PeakHeight = RunOutput[BestFitIndex][2]

	// store the data temporarily in case the user wants to accept the curve fit
	TempRawData[0] = X_Wave[Point]
	TempRawData[1] = Y_Wave[Point]
	TempRawData[2] = MinimumX
	TempRawData[3] = MaximumX
	TempRawData[4] = PeakStart
	TempRawData[5] = PeakEnd
	TempRawData[6] = PeakHeight
	TempRawData[7] = mean(BaseOriginalWave,floor(PeaKStart),floor(PeakEnd))
	
	// print out some info in the command window about what is happening
	print "Event Length: "+num2str(PeakEnd - PeakStart)+" Event Height: " + num2str(PeakHeight)

	PeakHeight = mean(BaseOriginalWave,floor(PeaKStart),floor(PeakEnd))
	print "Height During Event Corrected to Mean Value: " +  num2str(PeakHeight) + " Residual Number: " + num2str(BaselinePeakGeometryCalculation(DummyWave,PeakStart,PeakEnd,PeakHeight))

	//--------------------------------------------------
	// adjust the graph to display the new curve fit	
	//--------------------------------------------------	
	SetScale/I x,MinimumX,MaximumX,CurveFitWave
	CurveFitWave = PeakHeight * (HeavisideApproximation(x,PeakStart) - HeavisideApproximation(x,PeakEnd))
	ModifyGraph /W=PointChooserPanel#PointChooser lstyle(CurveFitWave)=0,rgb(CurveFitWave)=(0,0,65535)

	
	Button SaveCurveFit,disable=0
End











Function CurveFitInclude(CtrlName) : ButtonControl
	String CtrlName

	// wave needed to store number as text for output 
	String FormattedNumber


	// wave containing the current data
 	Wave TempRawData = root:SingleParticleTraceData:Statistics:Stat_TempMatrix
 	 

	svar WaveNameData = root:SingleParticleTraceData:Statistics:StatWaveNames
	nvar WaveSizeData = root:SingleParticleTraceData:Statistics:StatWaveSize
	Wave/T XData = root:SingleParticleTraceData:Statistics:Stat_X
	Wave/T YData = root:SingleParticleTraceData:Statistics:Stat_Y
	Wave/T FitBeginData = root:SingleParticleTraceData:Statistics:Stat_CurveFitBegin
	Wave/T FitEndData = root:SingleParticleTraceData:Statistics:Stat_CurveFitEnd
	Wave/T EventStartData = root:SingleParticleTraceData:Statistics:Stat_EventStart
	Wave/T EventEndData = root:SingleParticleTraceData:Statistics:Stat_EventEnd
	Wave/T EventHeightData = root:SingleParticleTraceData:Statistics:Stat_EventHeight
	Wave/T MeanHeightData = root:SingleParticleTraceData:Statistics:Stat_MeanHeight

	WaveSizeData+=1

	Redimension/N=(WaveSizeData) XData
	Redimension/N=(WaveSizeData) YData
	Redimension/N=(WaveSizeData) FitBeginData
	Redimension/N=(WaveSizeData) FitEndData
	Redimension/N=(WaveSizeData) EventStartData
	Redimension/N=(WaveSizeData) EventEndData
	Redimension/N=(WaveSizeData) EventHeightData
	Redimension/N=(WaveSizeData) MeanHeightData


	sprintf FormattedNumber, "%.10f", TempRawData[0]
	XData[WaveSizeData - 1] = FormattedNumber

	sprintf FormattedNumber, "%.10f", TempRawData[1]
	YData[WaveSizeData - 1] = FormattedNumber

	sprintf FormattedNumber, "%.10f", TempRawData[2]
	FitBeginData[WaveSizeData - 1] = FormattedNumber

	sprintf FormattedNumber, "%.10f", TempRawData[3]
	FitEndData[WaveSizeData - 1] = FormattedNumber

	sprintf FormattedNumber, "%.10f", TempRawData[4]
	EventStartData[WaveSizeData - 1] = FormattedNumber

	sprintf FormattedNumber, "%.10f", TempRawData[5]
	EventEndData[WaveSizeData - 1] = FormattedNumber

	sprintf FormattedNumber, "%.10f", TempRawData[6]
	EventHeightData[WaveSizeData - 1] = FormattedNumber

	sprintf FormattedNumber, "%.10f", TempRawData[7]
	MeanHeightData[WaveSizeData - 1] = FormattedNumber

	Button SaveCurveFit,disable=2	


End






//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
// special functions for either simulated annealing or curve fitting
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
// fit function for a square wave to the current wave under examination
// assumption of this version is that we are working with the baseline adjusted plot
// given this, additional baseline variable is not included in the calculation  
// subtract trial wave, integrate the residue, return the absolute value

Function BaselinePeakGeometryCalculation(DummyWave,SquareRise,SquareFall,SquareHeight)
	Wave DummyWave
	Variable SquareRise,SquareFall,SquareHeight

	nvar Point = root:PointIndex
	nvar MinimumX = root:SingleParticleTraceData:MinCurveFitFrame
	nvar MaximumX = root:SingleParticleTraceData:MaxCurveFitFrame
	Wave/T BaseWaveNameList = root:BaselineWaveNameArray
	Wave BaseOriginalWave = $("root:SingleParticleTraceData:" + BaseWaveNameList[Point])

	// first, if the SquareRise value is larger than or equal to the SquareFall, the 
	// curve fit is examining a nonsense solution, return a large postive number 
	// immediately to remove this possibility 
	if(SquareRise >= SquareFall)
		return 100000
	endif

	Variable WaveSize = numpnts($("root:SingleParticleTraceData:" + BaseWaveNameList[Point]))

	Wave Square = SquareWavePulse(WaveSize,SquareRise,SquareFall,0,SquareHeight)

	Make /O/N=(WaveSize) ToIntegrate 

	ToIntegrate = BaseOriginalWave - Square


	// crop the wave to fit the area of interest if the maximum is one less than the number of points in the wave

	// crop the end of the wave if the MaximumX is less than the maximum possible
	if(MaximumX < (numpnts(ToIntegrate) - 1))
		redimension /N=(MaximumX + 1) ToIntegrate
	endif	

	// crop the beginning if the minimum is less than the size of the wave
	// and also more the zero

	if((MinimumX < (numpnts(ToIntegrate) - 1))&&(MinimumX > 0))
		DeletePoints 0, MinimumX, ToIntegrate
	endif

	Integrate ToIntegrate /D=ReturnValue 


	// if the wave matches the experimental data, subtracting the two should lead to two things
	// 1. the integral should not vary much in time (the maximum minus the minimum should be low)
	// 2. the overall integral should be as close to zero as possible

	// condition 1
	Variable Condition1 = abs(wavemax(ReturnValue) - wavemin(ReturnValue))
	// condition 2
	Variable Condition2 = abs(ReturnValue(numPnts(ReturnValue) - 1))

	return Condition1 + Condition2
	

End



// function wrapper to use the standard curve futting routine
Function SquareWavePulseFit(w,x) : FitFunc
WAVE w
Variable x

// 0 is rise, 1 is fall, 2 is baseline height, 3 is peak height relative to baseline

return w[2] + w[3] * (HeavisideApproximation(x,w[0]) - HeavisideApproximation(x,w[1]))

End


//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
// functions that define the square wave fit
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
//generate a square wave pulse wave with fit parameters to subtract from real data
Function/WAVE SquareWavePulse(NumberPoints,SquareRise,SquareFall,BaselineHeight,SquareHeight)
	variable NumberPoints,SquareRise,SquareFall,BaselineHeight,SquareHeight

	String savedDataFolder = GetDataFolder(1)	// Save
	SetDataFolder root:SingleParticleTraceData:

	Make /O/N=(NumberPoints) SquareWave

	variable itor
	for(itor = 0; itor < NumberPoints; itor++)
		SquareWave[itor] = BaselineHeight + SquareHeight * (HeavisideApproximation(itor,SquareRise) - HeavisideApproximation(itor,SquareFall))
	endfor

	SetDataFolder savedDataFolder

	return SquareWave

End

//smooth approximation of a heaviside function
Function HeavisideApproximation(X_Point,Offset)
	variable X_Point,Offset

	return (1 / (1 + exp(-1000*(X_Point - Offset)))) 

End


//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------






















//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
// section for the right hand buttons in the single trace plots 
// first navigation and the like as well as some helper functions
//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------



Function RescalePlot()

	nvar Point = root:PointIndex
	Wave/T WaveNameList = root:WaveNameArray
	Wave OriginalWave = $("root:SingleParticleTraceData:" + WaveNameList[Point])
	Wave CurveFitWave = root:SingleParticleTraceData:CurveFitWave


	variable MaxRange, MinRange, range, delta


	// if the curve fit is basicall set to zero (the curve has no variance at all, just fit the current trace)
	if((wavemax(CurveFitWave) - wavemin(CurveFitWave)) < .0000001)

		// zoom out by 5 percent on the vertical data to avoid hitting the top and bottom for ease of viewing
		MaxRange =  wavemax(OriginalWave)
		MinRange = wavemin(OriginalWave)
		range = MaxRange - MinRange
		delta = range * 0.05 	
		SetAxis /W=PointChooserPanel#PointChooser left,(MinRange - delta),(MaxRange + delta)	
	else
		MaxRange = max(wavemax(OriginalWave),wavemax(CurveFitWave))
		MinRange = min(wavemin(OriginalWave),wavemin(CurveFitWave))
		range = MaxRange - MinRange
		delta = range * 0.05 	
		SetAxis /W=PointChooserPanel#PointChooser left,(MinRange - delta),(MaxRange + delta)
	Endif

End

Function RescaleBaselinePlot()

	nvar Point = root:PointIndex
	Wave/T BaseWaveNameList = root:BaselineWaveNameArray
	Wave BaseOriginalWave = $("root:SingleParticleTraceData:" + BaseWaveNameList[Point])
	Wave CurveFitWave = root:SingleParticleTraceData:CurveFitWave


	variable MaxRange, MinRange, range, delta


	// if the curve fit is basicall set to zero (the curve has no variance at all, just fit the current trace)
	if((wavemax(CurveFitWave) - wavemin(CurveFitWave)) < .0000001)

		// zoom out by 5 percent on the vertical data to avoid hitting the top and bottom for ease of viewing
		MaxRange =  wavemax(BaseOriginalWave)
		MinRange = wavemin(BaseOriginalWave)
		range = MaxRange - MinRange
		delta = range * 0.05 	
		SetAxis /W=PointChooserPanel#PointChooser left,(MinRange - delta),(MaxRange + delta)	
	else
		MaxRange = max(wavemax(BaseOriginalWave),wavemax(CurveFitWave))
		MinRange = min(wavemin(BaseOriginalWave),wavemin(CurveFitWave))
		range = MaxRange - MinRange
		delta = range * 0.05 	
		SetAxis /W=PointChooserPanel#PointChooser left,(MinRange - delta),(MaxRange + delta)
	Endif

End



Function ReturnMain(CtrlName) : ButtonControl
	String CtrlName


	// check if wave filtering has been applied and that there is at least one accepted wave before making a new subtraction filter menu item available
	if(WaveExists(root:KeepWaveArray))
		wave/I KeepWaveList = root:KeepWaveArray
		variable itor
		variable Accepted = 0
		for(itor = 0; itor < dimsize(KeepWaveList,0) - 1; itor++)
			if(KeepWaveList[itor] > 0.5)
				Accepted = 1
				break
			endif
		endfor
		if(Accepted > 0.5)
			String/G root:ExportedStr = "Subtract Background"
			BuildMenu "Localizer"
		endif	
	endif

	// kill all the waves in the single particle trace data folder, user should have saved them
	// and they are slowing everything down anyways
	SetDataFolder root:SingleParticleTraceData
	KillWaves/A/Z
	SetDataFolder root:

	SetDataFolder root:SingleParticleTraceData:Statistics
	KillWaves/A/Z
	SetDataFolder root:

	KillWindow/Z PointChooserPanel
End


Function PreviousPlot(CtrlName) : ButtonControl
	String CtrlName
	RunPreviousPlot()
End

Function NextPlot(CtrlName) : ButtonControl
	String CtrlName
	RunNextPlot()
End

Function AcceptPlot(CtrlName) : ButtonControl
	String CtrlName
	RunAcceptPlot()
End

Function RejectPlot(CtrlName) : ButtonControl
	String CtrlName
	RunRejectPlot()
End

Function ChangeCurrentTrace(sva) : SetVariableControl
    STRUCT WMSetVariableAction &sva
    if (sva.eventCode == 1 || sva.eventCode == 2 || sva.eventCode == 3)
    	variable NewIndex = sva.dval
    	RunChangePlot(NewIndex)
    endif
End

Function RunAcceptPlot()
	NVAR Index = root:PointIndex
	wave/I KeepWaveList = root:KeepWaveArray

	KeepWaveList[Index] = 1
	ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(40863,65535,40863)
	Button ApplyCurveFit,disable=0
	Button SaveCurveFit,disable=2
End

Function RunRejectPlot()
	NVAR Index = root:PointIndex
	wave/I KeepWaveList = root:KeepWaveArray

	KeepWaveList[Index] = -1
	ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,40863,40863)
	Button ApplyCurveFit,disable=2
	Button SaveCurveFit,disable=2
	//RunNextPlot()
End


Function RunNextPlot()
	
	NVAR Index = root:PointIndex
	NVAR IndexOld = root:PointIndexOld
	Wave/T WaveList = root:WaveNameArray
	Wave/T BaseWaveList = root:BaselineWaveNameArray
	wave/T WaveColorList = root:WaveColorArray
	wave/I KeepWaveList = root:KeepWaveArray
	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color
	Wave CurveFitWave = root:SingleParticleTraceData:CurveFitWave
	
	IndexOld = Index

	// loop back to beginning if at the end
	if((Index >= (dimsize(KeepWaveList,0) - 1))||((Index +1) > (dimsize(KeepWaveList,0) - 1))) // be careful, the index is a floating point number
		Index = 0
	else
		Index+=1
	Endif

	if(WavePlottedList[Index] == 0)	
		RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, Index)
	endif
		
	SetDataFolder root:SingleParticleTraceData



	// zoom the vertical axis properly
	// zoom out by 5 percent on the vertical data tp avoid hitting the top and bottom for ease of viewing 
	variable range = wavemax($WaveList[Index]) - wavemin ($WaveList[Index])
	variable delta = range * 0.05	
	SetAxis /W=PointChooserPanel#PointChooser left,(wavemin ($WaveList[Index]) - delta),(wavemax($WaveList[Index]) + delta)
	ReplaceWave /W=PointChooserPanel#PointChooser trace=$WaveList[IndexOld], $WaveList[Index]
	ReplaceWave /W=PointChooserPanel#PointChooser trace=$BaseWaveList[IndexOld], $BaseWaveList[Index]

	// paint background white (undecided), extra care here because we aren't dealing with integers
	if(abs(KeepWaveList[Index]) < 0.5)		
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,65535,65535)
		Button ApplyCurveFit,disable=2
		Button SaveCurveFit,disable=2

	// paint background redish (throwing it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[Index]) > 0.5)&&(KeepWaveList[Index] < 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,40863,40863)
		Button ApplyCurveFit,disable=2
		Button SaveCurveFit,disable=2

	// paint background greenish (keepin it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[Index]) > 0.5)&&(KeepWaveList[Index] > 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(40863,65535,40863)
		Button ApplyCurveFit,disable=0
		Button SaveCurveFit,disable=2
	
	endif

	//-----------------------------------
	// modifications to the wave lines
	//-----------------------------------
	// color the wave itself to indicate the point that localizer chose acording to WaveColorList
	ModifyGraph /W=PointChooserPanel#PointChooser zColor($WaveList[Index])={$WaveColorList[Index],*,*,directRGB} 

	// zero the wave, new trace 
	CurveFitWave = 0
	// set the type (dashed)
	ModifyGraph /W=PointChooserPanel#PointChooser lstyle(CurveFitWave)=3,rgb(CurveFitWave)=(0,0,65535)

	IndexOld = Index

	SetDataFolder root:
End

Function RunPreviousPlot()

	NVAR Index = root:PointIndex
	NVAR IndexOld = root:PointIndexOld
	wave/T WaveList = root:WaveNameArray
	wave/T BaseWaveList = root:BaselineWaveNameArray
	wave/T WaveColorList = root:WaveColorArray
	wave/I KeepWaveList = root:KeepWaveArray
	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color
	Wave CurveFitWave = root:SingleParticleTraceData:CurveFitWave

	IndexOld = Index

	// loop to end of list if at the very beginning
	if((Index <= 0) || ((Index-1) < 0)) // be careful, the index is a floating point number
		Index = dimsize(KeepWaveList,0) - 1
	else
		Index-=1
	Endif


	if(WavePlottedList[Index] == 0)	
		RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, Index)
	endif


	SetDataFolder root:SingleParticleTraceData

	// zoom the vertical axis properly
	// zoom out by 5 percent on the vertical data tp avoid hitting the top and bottom for ease of viewing 
	variable range = wavemax($WaveList[Index]) - wavemin ($WaveList[Index])
	variable delta = range * 0.05	
	SetAxis /W=PointChooserPanel#PointChooser left,(wavemin ($WaveList[Index]) - delta),(wavemax($WaveList[Index]) + delta)
	// Plot the previous wave
	ReplaceWave /W=PointChooserPanel#PointChooser trace=$WaveList[IndexOld], $WaveList[Index]
	ReplaceWave /W=PointChooserPanel#PointChooser trace=$BaseWaveList[IndexOld], $BaseWaveList[Index]

	// paint background white (undecided), extra care here because we aren't dealing with integers
	if(abs(KeepWaveList[Index]) < 0.5)	
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,65535,65535)
		Button ApplyCurveFit,disable=2
		Button SaveCurveFit,disable=2

	// paint background redish (throwin it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[Index]) > 0.5)&&(KeepWaveList[Index] < 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,40863,40863)
		Button ApplyCurveFit,disable=2
		Button SaveCurveFit,disable=2

	// paint background greenish (keepin it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[Index]) > 0.5)&&(KeepWaveList[Index] > 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(40863,65535,40863)
		Button ApplyCurveFit,disable=0
		Button SaveCurveFit,disable=2

	endif


	//-----------------------------------
	// modifications to the wave lines
	//-----------------------------------

	// color the wave itself to indicate the point that localizer chose acording to WaveColorList
	ModifyGraph /W=PointChooserPanel#PointChooser zColor($WaveList[Index])={$WaveColorList[Index],*,*,directRGB} 

	// zero the wave, new trace 
	CurveFitWave = 0
	// set the type (dashed)
	ModifyGraph /W=PointChooserPanel#PointChooser lstyle(CurveFitWave)=3,rgb(CurveFitWave)=(0,0,65535)

	IndexOld = Index

	SetDataFolder root:
End



Function RunChangePlot(NewPlot)
	variable NewPlot

	NVAR Index = root:PointIndex
	NVAR IndexOld = root:PointIndexOld
	Wave/T WaveList = root:WaveNameArray
	Wave/T BaseWaveList = root:BaselineWaveNameArray
	wave/T WaveColorList = root:WaveColorArray
	wave/I KeepWaveList = root:KeepWaveArray
	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color
	Wave CurveFitWave = root:SingleParticleTraceData:CurveFitWave
	

	if(WavePlottedList[NewPlot] == 0)	
		RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, NewPlot)
	endif
		
	SetDataFolder root:SingleParticleTraceData



	// zoom the vertical axis properly
	// zoom out by 5 percent on the vertical data tp avoid hitting the top and bottom for ease of viewing 
	variable range = wavemax($WaveList[NewPlot]) - wavemin ($WaveList[NewPlot])
	variable delta = range * 0.05	
	SetAxis /W=PointChooserPanel#PointChooser left,(wavemin ($WaveList[NewPlot]) - delta),(wavemax($WaveList[NewPlot]) + delta)
	ReplaceWave /W=PointChooserPanel#PointChooser trace=$WaveList[IndexOld], $WaveList[NewPlot]
	ReplaceWave /W=PointChooserPanel#PointChooser trace=$BaseWaveList[IndexOld], $BaseWaveList[NewPlot]

	// paint background white (undecided), extra care here because we aren't dealing with integers
	if(abs(KeepWaveList[NewPlot]) < 0.5)		
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,65535,65535)
		Button ApplyCurveFit,disable=2
		Button SaveCurveFit,disable=2

	// paint background redish (throwing it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[NewPlot]) > 0.5)&&(KeepWaveList[NewPlot] < 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(65535,40863,40863)
		Button ApplyCurveFit,disable=2
		Button SaveCurveFit,disable=2

	// paint background greenish (keepin it), extra care here because we aren't dealing with integers
	elseif((abs(KeepWaveList[NewPlot]) > 0.5)&&(KeepWaveList[NewPlot] > 0))
		ModifyGraph /W=PointChooserPanel#PointChooser gbRGB=(40863,65535,40863)
		Button ApplyCurveFit,disable=0
		Button SaveCurveFit,disable=2

	endif

	//-----------------------------------
	// modifications to the wave lines
	//-----------------------------------
	// color the wave itself to indicate the point that localizer chose acording to WaveColorList
	ModifyGraph /W=PointChooserPanel#PointChooser zColor($WaveList[NewPlot])={$WaveColorList[NewPlot],*,*,directRGB} 

	// zero the wave, new trace 
	CurveFitWave = 0
	// set the type (dashed)
	ModifyGraph /W=PointChooserPanel#PointChooser lstyle(CurveFitWave)=3,rgb(CurveFitWave)=(0,0,65535)

	IndexOld = NewPlot

	SetDataFolder root:
End

Function KeyboardWindowHook(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0 // 0 if we do not handle event, 1 if we handle it.


	switch(s.eventCode)
		case 11:  // Keyboard event

			if(s.keycode==28) //left arrow
				RunPreviousPlot()
				hookResult = 1

			elseif(s.keycode==29) //right arrow
				RunNextPlot()
				hookResult = 1

			elseif(s.keycode==30) //up arrow
				RunAcceptPlot()
				hookResult = 1

			elseif(s.keycode==31) //down arrow
				RunRejectPlot()
				hookResult = 1

			Endif

	endswitch
	return hookResult // If non-zero, we handled event and Igor will ignore it.
End










//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
// a series of data exporting functions 
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------


Function ExportCurveFits(CtrlName) : ButtonControl
	String CtrlName
	
	svar TempRawDataNames = root:SingleParticleTraceData:Statistics:StatWaveNames

	nvar WaveSizeData = root:SingleParticleTraceData:Statistics:StatWaveSize

	// check to see if there is a curve fit to include before doing so
	
	SetDataFolder root:SingleParticleTraceData:Statistics

		if(WaveSizeData > 1.5) // first point only contains column headers
			Save/J/B TempRawDataNames
		endif	

	SetDataFolder root:
End	


Function ExportBaselined(CtrlName) : ButtonControl
	String CtrlName

	// name of all the stored waves
	//wave/T WaveList = root:WaveNameArray
	Wave/T BaseWaveNameList = root:BaselineWaveNameArray
	wave WaveList_X = root:WaveNameArray_X
	wave WaveList_Y = root:WaveNameArray_Y
	// integer flag to indicate filter status
	wave/I KeepWaveList = root:KeepWaveArray

	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color

	

	variable itor = 0
	variable numberAccepted=0

	

	// wave to store the names of all the waves that were accepted
	String AcceptedWaves

	// get all of the accepted waves into one list
	for(itor = 0; itor < numpnts(BaseWaveNameList); itor++)
		// add one to the list if it is in the accepted pile
		if((abs(KeepWaveList[itor])>0.5)&&(KeepWaveList[itor] > 0))

			//-------------------------------------------------------
			// make sure the wave is created before manipulating it
			//-------------------------------------------------------
			if(WavePlottedList[itor] == 0)	
				RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, itor)
			endif

			SetDataFolder root:SingleParticleTraceData

			numberAccepted++
			// add wave name to growing string which contains references to the waves
			// duplicate the wave to add an 'X', 'Y' header to it
			if(numberAccepted < 1.5)

				Make /T/O/N=0 ModifiedWaveList
				// first create a wave listing the frames as the X-coordinate wave to start in the list
				// also add a header to indicate that the first two numbers in the following columns are the x,y coordinates
				Make /T/O/N=(numpnts($BaseWaveNameList[itor]) + 2) FrameList
				variable frameItor
				FrameList[0] = "X"
				FrameList[1] = "Y" 
				for(frameItor = 2; frameItor < dimsize(FrameList,0); frameItor++)
					FrameList[frameItor] = num2istr(frameItor - 2)
				endfor

				// now add the first wave name
				AcceptedWaves = "FrameList;" + BaseWaveNameList[itor] + "_A"
				// list for later extraneous wave deletion
				Redimension/N=(numberAccepted) ModifiedWaveList
				ModifiedWaveList[numberAccepted - 1] = BaseWaveNameList[itor] + "_A"
				// modify the wave with an x, y header
				Duplicate $BaseWaveNameList[itor], $(BaseWaveNameList[itor] + "_A")
				wave temp_edit = $(BaseWaveNameList[itor] + "_A")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit
				temp_edit[0] = WaveList_X[itor]
				temp_edit[1] = WaveList_Y[itor]

			elseif(numberAccepted > 1.5)

				AcceptedWaves = AcceptedWaves + ";" + BaseWaveNameList[itor] + "_A"
				// list for later extraneous wave deletion
				Redimension/N=(numberAccepted) ModifiedWaveList
				ModifiedWaveList[numberAccepted - 1] = BaseWaveNameList[itor] + "_A"
				Duplicate $BaseWaveNameList[itor], $(BaseWaveNameList[itor] + "_A")
				wave temp_edit = $(BaseWaveNameList[itor] + "_A")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit
				temp_edit[0] = WaveList_X[itor]
				temp_edit[1] = WaveList_Y[itor]

			endif

			SetDataFolder root:

		endif

	endfor

	// if there are accepted waves, modify for saving and then save them

	SetDataFolder root:SingleParticleTraceData

	if(numberAccepted > 0.5)
		Save/J/B AcceptedWaves
		// remove the modified waves
		for(itor = 0; itor < dimsize(ModifiedWaveList,0); itor++)
			killWaves $ModifiedWaveList[itor]
		endfor
		killWaves ModifiedWaveList
		killWaves FrameList
	endif

	SetDataFolder root:


End




Function ExportAccept(CtrlName) : ButtonControl
	String CtrlName

	// name of all the stored waves
	wave/T WaveList = root:WaveNameArray
	wave WaveList_X = root:WaveNameArray_X
	wave WaveList_Y = root:WaveNameArray_Y
	// integer flag to indicate filter status
	wave/I KeepWaveList = root:KeepWaveArray

	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color

	

	variable itor = 0
	variable numberAccepted=0

	

	// wave to store the names of all the waves that were accepted
	String AcceptedWaves

	// get all of the accepted waves into one list
	for(itor = 0; itor < dimsize(WaveList,0); itor++)
		// add one to the list if it is in the accepted pile
		if((abs(KeepWaveList[itor])>0.5)&&(KeepWaveList[itor] > 0))

			//-------------------------------------------------------
			// make sure the wave is created before manipulating it
			//-------------------------------------------------------
			if(WavePlottedList[itor] == 0)	
				RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, itor)
			endif

			SetDataFolder root:SingleParticleTraceData

			numberAccepted++
			// add wave name to growing string which contains references to the waves
			// duplicate the wave to add an 'X', 'Y' header to it
			if(numberAccepted < 1.5)

				Make /T/O/N=0 ModifiedWaveList
				// first create a wave listing the frames as the X-coordinate wave to start in the list
				// also add a header to indicate that the first two numbers in the following columns are the x,y coordinates
				Make /T/O/N=(dimsize($WaveList[itor],0) + 2) FrameList
				variable frameItor
				FrameList[0] = "X"
				FrameList[1] = "Y" 
				for(frameItor = 2; frameItor < dimsize(FrameList,0); frameItor++)
					FrameList[frameItor] = num2istr(frameItor - 2)
				endfor

				// now add the first wave name
				AcceptedWaves = "FrameList;" + WaveList[itor] + "_A"
				// list for later extraneous wave deletion
				Redimension/N=(numberAccepted) ModifiedWaveList
				ModifiedWaveList[numberAccepted - 1] = WaveList[itor] + "_A"
				// modify the wave with an x, y header
				Duplicate $WaveList[itor], $(WaveList[itor] + "_A")
				wave temp_edit = $(WaveList[itor] + "_A")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit

				//temp_edit[0] = WaveList_X[itor]
				//temp_edit[1] = WaveList_Y[itor]

				temp_edit[0] = X_Wave[itor]
				temp_edit[1] = Y_Wave[itor]

			elseif(numberAccepted > 1.5)

				AcceptedWaves = AcceptedWaves + ";" + WaveList[itor] + "_A"
				// list for later extraneous wave deletion
				Redimension/N=(numberAccepted) ModifiedWaveList
				ModifiedWaveList[numberAccepted - 1] = WaveList[itor] + "_A"
				Duplicate $WaveList[itor], $(WaveList[itor] + "_A")
				wave temp_edit = $(WaveList[itor] + "_A")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit

				//temp_edit[0] = WaveList_X[itor]
				//temp_edit[1] = WaveList_Y[itor]

				temp_edit[0] = X_Wave[itor]
				temp_edit[1] = Y_Wave[itor]

			endif

			SetDataFolder root:

		endif

	endfor

	// if there are accepted waves, modify for saving and then save them

	SetDataFolder root:SingleParticleTraceData

	if(numberAccepted > 0.5)
		Save/J/B AcceptedWaves
		// remove the modified waves
		for(itor = 0; itor < dimsize(ModifiedWaveList,0); itor++)
			killWaves $ModifiedWaveList[itor]
		endfor
		killWaves ModifiedWaveList
		killWaves FrameList
	endif

	SetDataFolder root:


End

Function ExportReject(CtrlName) : ButtonControl
	String CtrlName

	// name of all the stored waves
	wave/T WaveList = root:WaveNameArray
	wave WaveList_X = root:WaveNameArray_X
	wave WaveList_Y = root:WaveNameArray_Y
	// integer flag to indicate filter status
	wave/I KeepWaveList = root:KeepWaveArray

	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color


	variable itor = 0
	variable numberRejected=0

	

	// wave to store the names of all the waves that were accepted
	String RejectedWaves

	// get all of the rejected waves into one list
	for(itor = 0; itor < dimsize(WaveList,0); itor++)
		// add one to the list if it is in the accepted pile
		if((abs(KeepWaveList[itor])>0.5)&&(KeepWaveList[itor] < 0))

			//-------------------------------------------------------
			// make sure the wave is created before manipulating it
			//-------------------------------------------------------
			if(WavePlottedList[itor] == 0)	
				RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, itor)
			endif

			SetDataFolder root:SingleParticleTraceData

			numberRejected++
			// add wave name to growing string which contains references to the waves
			// duplicate the wave to add an 'X', 'Y' header to it
			if(numberRejected < 1.5)

				Make /T/O/N=0 ModifiedWaveList
				// first create a wave listing the frames as the X-coordinate wave to start in the list
				// also add a header to indicate that the first two numbers in the following columns are the x,y coordinates
				Make /T/O/N=(dimsize($WaveList[itor],0) + 2) FrameList
				variable frameItor
				FrameList[0] = "X"
				FrameList[1] = "Y" 
				for(frameItor = 2; frameItor < dimsize(FrameList,0); frameItor++)
					FrameList[frameItor] = num2istr(frameItor - 2)
				endfor

				// now add the first wave name
				RejectedWaves = "FrameList;" + WaveList[itor] + "_R"
				// list for later extraneous wave deletion
				Redimension/N=(numberRejected) ModifiedWaveList
				ModifiedWaveList[numberRejected - 1] = WaveList[itor] + "_R"
				// modify the wave with an x, y header
				Duplicate $WaveList[itor], $(WaveList[itor] + "_R")
				wave temp_edit = $(WaveList[itor] + "_R")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit

				//temp_edit[0] = WaveList_X[itor]
				//temp_edit[1] = WaveList_Y[itor]

				temp_edit[0] = X_Wave[itor]
				temp_edit[1] = Y_Wave[itor]

			elseif(numberRejected > 1.5)

				RejectedWaves = RejectedWaves + ";" + WaveList[itor] + "_R"
				// list for later extraneous wave deletion
				Redimension/N=(numberRejected) ModifiedWaveList
				ModifiedWaveList[numberRejected - 1] = WaveList[itor] + "_R"
				Duplicate $WaveList[itor], $(WaveList[itor] + "_R")
				wave temp_edit = $(WaveList[itor] + "_R")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit

				//temp_edit[0] = WaveList_X[itor]
				//temp_edit[1] = WaveList_Y[itor]

				temp_edit[0] = X_Wave[itor]
				temp_edit[1] = Y_Wave[itor]

			endif

			SetDataFolder root:
		
		endif

	endfor

	// if there are rejected waves, modify for saving and then save them

	SetDataFolder root:SingleParticleTraceData


	if(numberRejected > 0.5)
		Save/J/B RejectedWaves
		// remove the modified waves
		for(itor = 0; itor < dimsize(ModifiedWaveList,0); itor++)
			killWaves $ModifiedWaveList[itor]
		endfor
		killWaves ModifiedWaveList
		killWaves FrameList
	endif

	SetDataFolder root:


End

Function ExportUndecided(CtrlName) : ButtonControl
	String CtrlName


	// name of all the stored waves
	wave/T WaveList = root:WaveNameArray
	wave WaveList_X = root:WaveNameArray_X
	wave WaveList_Y = root:WaveNameArray_Y
	// integer flag to indicate filter status
	wave/I KeepWaveList = root:KeepWaveArray

	wave/I WavePlottedList = root:WavePlottedArray
	wave/I X_Wave = root:Filter_X
	wave/I Y_Wave = root:Filter_Y
	wave/I Frame_Wave = root:Filter_Color

	

	variable itor = 0
	variable numberNeither=0

	

	// wave to store the names of all the waves that were neither
	// this string list will be used to save the waves sequentially
	String NeitherWaves

	// get all of the accepted waves into one list
	for(itor = 0; itor < dimsize(WaveList,0); itor++)
		// add one to the list if it is in the neither pile
		if(abs(KeepWaveList[itor])<0.5)

			//-------------------------------------------------------
			// make sure the wave is created before manipulating it
			//-------------------------------------------------------
			if(WavePlottedList[itor] == 0)	
				RetrieveIntensityPlot(X_Wave, Y_Wave, Frame_Wave, itor)
			endif

			SetDataFolder root:SingleParticleTraceData

			numberNeither++
			// add wave name to growing string which contains references to the waves
			// duplicate the wave to add an 'X', 'Y' header to it
			if(numberNeither < 1.5)

				Make /T/O/N=0 ModifiedWaveList

				// first create a wave listing the frames as the X-coordinate wave to start in the list
				// also add a header to indicate that the first two numbers in the following columns are the x,y coordinates
				Make /T/O/N=(dimsize($WaveList[itor],0) + 2) FrameList
				variable frameItor
				FrameList[0] = "X"
				FrameList[1] = "Y" 

				// populate the first column with frames
				for(frameItor = 2; frameItor < dimsize(FrameList,0); frameItor++)
					FrameList[frameItor] = num2istr(frameItor - 2)
				endfor

				// now add the first wave name
				NeitherWaves = "FrameList;" + WaveList[itor] + "_N"
				// list for later extraneous wave deletion
				Redimension/N=(numberNeither) ModifiedWaveList
				ModifiedWaveList[numberNeither - 1] = WaveList[itor] + "_N"
				// modify the wave with an x, y header
				Duplicate $WaveList[itor], $(WaveList[itor] + "_N")
				wave temp_edit = $(WaveList[itor] + "_N")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit

				//temp_edit[0] = WaveList_X[itor]
				//temp_edit[1] = WaveList_Y[itor]

				temp_edit[0] = X_Wave[itor]
				temp_edit[1] = Y_Wave[itor]

			elseif(numberNeither > 1.5)

				NeitherWaves = NeitherWaves + ";" + WaveList[itor] + "_N"
				// list for later extraneous wave deletion
				Redimension/N=(numberNeither) ModifiedWaveList
				ModifiedWaveList[numberNeither - 1] = WaveList[itor] + "_N"
				Duplicate $WaveList[itor], $(WaveList[itor] + "_N")
				wave temp_edit = $(WaveList[itor] + "_N")
				// add in two points in the beginning of the new wave to store the x and y coordinates
				InsertPoints 0,2,temp_edit

				//temp_edit[0] = WaveList_X[itor]
				//temp_edit[1] = WaveList_Y[itor]

				temp_edit[0] = X_Wave[itor]
				temp_edit[1] = Y_Wave[itor]

			endif

			SetDataFolder root:


		endif

	endfor

	// if there are neither waves, modify for saving and then save them

	SetDataFolder root:SingleParticleTraceData

	if(numberNeither > 0.5)
		Save/J/B NeitherWaves
		// remove the modified waves
		for(itor = 0; itor < dimsize(ModifiedWaveList,0); itor++)
			killWaves $ModifiedWaveList[itor]
		endfor
		killWaves ModifiedWaveList
		killWaves FrameList
	endif

	SetDataFolder root:	

End

















//-------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------
// section associated with the two channel comparison functions
//-------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------



	
Function Get3by3BeamAverage()
	
	NVAR CurrentPointIndex = root:SubtractionInterface:CurrentPointIndex
	// global variable references for the image channels
	svar ChannelTwoReference
	svar ChannelOneReference
	// global variable references for the maximum and minumum for the plot cursor
	nvar MaxAxis = root:SubtractionInterface:MaxLeftAxis
	nvar MinAxis = root:SubtractionInterface:MinLeftAxis

	wave X_Points = root:WaveNameArray_X
	wave Y_Points = root:WaveNameArray_Y

	Wave ChannelOne=$ChannelOneReference
	Wave ChannelTwo=$ChannelTwoReference

	Wave Ch1_Average = root:SubtractionInterface:ChannelOne_Average
	Wave Ch2_Average = root:SubtractionInterface:ChannelTwo_Average

	variable itor_X, itor_Y

	Ch1_Average = 0

	// sum all of the waves in a three by three grid for channel one
	for(itor_X = -1; itor_X < 2; itor_X++)
		for(itor_Y = -1; itor_Y < 2; itor_Y++)
			ImageTransform/beam={(X_Points[CurrentPointIndex] + itor_X), (Y_Points[CurrentPointIndex] + itor_Y)} getbeam ChannelOne
			wave W_Beam
			Ch1_Average = Ch1_Average + W_Beam
		endfor
	endfor
	Ch1_Average = Ch1_Average / 9.0			

	Ch2_Average = 0

	// sum all of the waves in a three by three grid for channel two
	for(itor_X = -1; itor_X < 2; itor_X++)
		for(itor_Y = -1; itor_Y < 2; itor_Y++)
			ImageTransform/beam={(X_Points[CurrentPointIndex] + itor_X), (Y_Points[CurrentPointIndex] + itor_Y)} getbeam ChannelTwo
			wave W_Beam
			Ch2_Average = Ch2_Average + W_Beam
		endfor
	endfor
	Ch2_Average = Ch2_Average / 9.0			


	//--------------------------------------------------------------------------------
	// update the range for the traces as they are updated above and will have changed
	//--------------------------------------------------------------------------------
	// zoom out by 10 percent on the vertical data to avoid hitting the top and bottom for ease of viewing
	// get the maximum value for the two waves displayed
	variable maxVal = max((wavemax(Ch1_average)), (wavemax(Ch2_average)))
	variable minVal = min((wavemin(Ch1_average)), (wavemin(Ch2_average)))
	variable range = maxVal - minVal
	variable delta = range * 0.1 / 2.0
	MaxAxis = maxVal + delta
	MinAxis = minVal - delta
	SetAxis /W=ChannelSubtraction#IntensityTraces/Z left,MinAxis,MaxAxis 	
End






//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
// Functions associated with the importation and analysis of curve fits one one channel with respect to another channel
//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------




Function CompareBackgroundChannel()

	// set up a place to store all of the data and/or clean it up if it already exists (maybe not)
	if(DataFolderExists("root:CompareBackgroundChannel"))
		SetDataFolder root:CompareBackgroundChannel
		KillWaves/A/Z
	else
		NewDataFolder root:CompareBackgroundChannel
		SetDataFolder root:CompareBackgroundChannel
	endif


	// set up some important global variables	

	// references to the images wave
	String/G ChannelOneCompareReference 
	String/G ChannelTwoCompareReference
	String/G TextFilePathComparison

	NewPanel /W=(497,102,780,220) /N=CompareBackgroundChannels
	Button ChannelOne,pos={10.00,20.00},size={80.00,40.00},title="Channel One\nImage",proc=LoadChannelOne
	Button ChannelTwo,pos={100.00,20.00},size={80.00,40.00},title="Channel Two\nImage",proc=LoadChannelTwo
	Button CurveFit,pos={190.00,20.00},size={80.00,40.00},title="Curve Fit\nFile",proc=LoadCurveFit
	Button CancelTwoChannel,pos={100.00,90.00},size={80.00,20.00},title="Cancel",proc=CancelCompare		

	CheckBox OneCheck,pos={45.00,70.00},size={11.00,11.00},disable=2,title=""
	CheckBox OneCheck,value= 0
	CheckBox TwoCheck,pos={135.00,70.00},size={11.00,11.00},disable=2,title=""
	CheckBox TwoCheck,value= 0
	CheckBox FitCheck,pos={225.00,70.00},size={11.00,11.00},disable=2,title=""
	CheckBox FitCheck,value= 0

End

Function LoadChannelOne(CtrlName) : ButtonControl
	String CtrlName

	svar ChannelOneRef = root:CompareBackgroundChannel:ChannelOneCompareReference	
	svar ChannelTwoRef = root:CompareBackgroundChannel:ChannelTwoCompareReference	

	//---------------------------------
	// open the first channel tif file
	//---------------------------------
	Variable refNum
	String message = "Select a tiff file for channel one"
	String ChannelOneComparePath
	String fileFilters = "Tiff Files (*.tif,*.tiff):.tif,.tiff;"
	Open /D /R /F=fileFilters /M=message refNum
	ChannelOneComparePath = S_FileName

	// open as an image
	if(numtype(strlen(ChannelOneComparePath)) != 2)	
		ImageLoad/T=tiff/S=0/C=-1/Q ChannelOneComparePath

		// only go further if succesfully opened
		if(V_flag == 1)
			ChannelOneRef= StringFromList(0,s_waveNames)
			//-------------------------------------------------------
			// if loaded before channel two don't compare dimensions
			// set the checkbox flag to selected
			// otherwise compare
			// and set as required
			//-------------------------------------------------------	
			ControlInfo /W=CompareBackgroundChannels TwoCheck
			if(V_value == 0) // 0 if unselected 1 if selected	
				// just check the box, no other channel to compare to
				CheckBox OneCheck,value= 1
			elseif(V_value == 1)
				// the other file has been loaded, do comparison
				// get the waves
				Wave ChannelOne=$ChannelOneRef
				Wave ChannelTwo=$ChannelTwoRef

				// load the dimension information for both channels
				variable Ch1_width = dimsize(ChannelOne, 0)
				variable Ch1_height = dimsize(ChannelOne, 1)
				variable Ch1_frames = dimsize(ChannelOne, 2)
				variable Ch2_width = dimsize(ChannelTwo, 0)
				variable Ch2_height = dimsize(ChannelTwo, 1)
				variable Ch2_frames = dimsize(ChannelTwo, 2)

				// do the size comparison
				if( (Ch1_width == Ch2_width)&&(Ch1_height == Ch2_height)&&(Ch1_frames == Ch2_frames))
					//-------------------------------------------------------------------------------------
					// images match, if the curve fit already loaded then don't bother with the checkbox
					// go directly to the comparison operations window
					// otherwise hit the checkbox for the channel
					//-------------------------------------------------------------------------------------
					ControlInfo /W=CompareBackgroundChannels FitCheck

					if(V_value == 0) // 0 if unselected 1 if selected
						//--------------------------------------------------------
						// files are the same length, the fits have *not* loaded
						// check the box for the channel
						//--------------------------------------------------------
						CheckBox OneCheck,value= 1	
					elseif(V_value == 1) // selected
						//--------------------------------------------------------
						// files are the same length, the fits have been loaded
						// close file panel and open the comparison panel
						//--------------------------------------------------------
						CompareBackgroundChannelPanel()
						KillWindow/Z CompareBackgroundChannels
					endif

				else
					DoAlert /T="Image Load Error" 0, "Channel One Image Does Not Match\nDimensions of Channel Two"
				endif
			Endif
		Endif
	Endif


End

Function LoadChannelTwo(CtrlName) : ButtonControl
	String CtrlName

	svar ChannelOneRef = root:CompareBackgroundChannel:ChannelOneCompareReference	
	svar ChannelTwoRef = root:CompareBackgroundChannel:ChannelTwoCompareReference	

	//---------------------------------
	// open the second channel tif file
	//---------------------------------
	Variable refNum
	String message = "Select a tiff file for channel two"
	String ChannelTwoComparePath
	String fileFilters = "Tiff Files (*.tif,*.tiff):.tif,.tiff;"
	Open /D /R /F=fileFilters /M=message refNum
	ChannelTwoComparePath = S_FileName

	// open as an image
	if(numtype(strlen(ChannelTwoComparePath)) != 2)	
		ImageLoad/T=tiff/S=0/C=-1/Q ChannelTwoComparePath

		// only go further if succesfully opened
		if(V_flag == 1)
			ChannelTwoRef= StringFromList(0,s_waveNames)
			//-------------------------------------------------------
			// if loaded before channel two don't compare dimensions
			// set the checkbox flag to selected
			// otherwise compare
			// and set as required
			//-------------------------------------------------------	
			ControlInfo /W=CompareBackgroundChannels OneCheck
			if(V_value == 0) // 0 if unselected 1 if selected	
				// just check the box, no other channel to compare to
				CheckBox TwoCheck,value= 1
			elseif(V_value == 1)
				// the other file has been loaded, do comparison
				// get the waves
				Wave ChannelOne=$ChannelOneRef
				Wave ChannelTwo=$ChannelTwoRef

				// load the dimension information for both channels
				variable Ch1_width = dimsize(ChannelOne, 0)
				variable Ch1_height = dimsize(ChannelOne, 1)
				variable Ch1_frames = dimsize(ChannelOne, 2)
				variable Ch2_width = dimsize(ChannelTwo, 0)
				variable Ch2_height = dimsize(ChannelTwo, 1)
				variable Ch2_frames = dimsize(ChannelTwo, 2)

				// do the size comparison
				if( (Ch1_width == Ch2_width)&&(Ch1_height == Ch2_height)&&(Ch1_frames == Ch2_frames))
					//-------------------------------------------------------------------------------------
					// images match, if the curve fit already loaded then don't bother with the checkbox
					// go directly to the comparison operations window
					// otherwise hit the checkbox for the channel
					//-------------------------------------------------------------------------------------
					ControlInfo /W=CompareBackgroundChannels FitCheck

					if(V_value == 0) // 0 if unselected 1 if selected
						//--------------------------------------------------------
						// files are the same length, the fits have *not* loaded
						// check the box for the channel
						//--------------------------------------------------------
						CheckBox TwoCheck,value= 1	
					elseif(V_value == 1) // selected
						//--------------------------------------------------------
						// files are the same length, the fits have been loaded
						// close file panel and open the comparison panel
						//--------------------------------------------------------
						CompareBackgroundChannelPanel()
						KillWindow/Z CompareBackgroundChannels
					endif

				else
					DoAlert /T="Image Load Error" 0, "Channel One Image Does Not Match\nDimensions of Channel Two"
				endif
			Endif
		Endif
	Endif
End

Function LoadCurveFit(CtrlName) : ButtonControl
	String CtrlName

	Variable refNum
	String message = "Select a curve fit text file"
	svar TextFilePath = root:CompareBackgroundChannel:TextFilePathComparison
	String fileFilters = "Text Files (*.txt):.txt;"
	Open /D/R /F=fileFilters /M=message refNum
	TextFilePath = S_FileName

	// make sure the path exists
	if(numtype(strlen(TextFilePath)) != 2)
		//load up the text file with 8 waves
		LoadWave /A/H/J/K=1/O/W TextFilePath
		// there should be 8 waves
		if(V_flag == 8)
			// check for the existence of the proper waves
			Variable Wv1,Wv2,Wv3,Wv4,Wv5,Wv6,Wv7,Wv8
			Wv1 = WaveExists(root:CompareBackgroundChannel:XW)
			Wv2 = WaveExists(root:CompareBackgroundChannel:YW)
			Wv3 = WaveExists(root:CompareBackgroundChannel:Curve_Fit_Start)
			Wv4 = WaveExists(root:CompareBackgroundChannel:Curve_Fit_End)
			Wv5 = WaveExists(root:CompareBackgroundChannel:Event_Start)
			Wv6 = WaveExists(root:CompareBackgroundChannel:Event_End)
			Wv7 = WaveExists(root:CompareBackgroundChannel:Event_Height)
			Wv8 = WaveExists(root:CompareBackgroundChannel:Event_Mean_Height)

			if(Wv1 && Wv2 && Wv3 && Wv4 && Wv5 && Wv6 && Wv7 && Wv8)
				//-----------------------------------------------
				//now check to see if the image waves are loaded
				//if both aren't, just hit checkbox and return
				//-----------------------------------------------
				variable OneTest, TwoTest
				ControlInfo /W=CompareBackgroundChannels OneCheck
				OneTest = V_value
				ControlInfo /W=CompareBackgroundChannels OneCheck
				TwoTest = V_value
				if(OneTest && TwoTest) 
					CompareBackgroundChannelPanel()
					KillWindow/Z CompareBackgroundChannels
				else 
					CheckBox FitCheck,value= 1
				endif

			else
				DoAlert /T="File Load Error" 0, "This does not appear to be a curve fit file"
			Endif

		Endif
	Endif

End

Function CancelCompare(CtrlName) : ButtonControl
	String CtrlName
	SetDataFolder root:CompareBackgroundChannel
	KillWaves/A/Z
	SetDataFolder root:
	KillWindow/Z CompareBackgroundChannels
End




Function CompareBackgroundChannelPanel()
	//------------------------------------------------------------------------------
	// first things first, do the calculations before we mess around with fancy
	// interfaces.
	// These were user-determined points so do it all at once, it shouldn't matter
	// much time wise
	//------------------------------------------------------------------------------
	Wave XWave = root:CompareBackgroundChannel:XW

	// references to the created 3x3 intensity traces
	Make /T/N=(numpnts(XWave)) ComparisonWaveChannelOne
	Make /T/N=(numpnts(XWave)) ComparisonWaveChannelTwo

	// waves to store the integration wave
	Make /T/N=(numpnts(XWave)) IntegrationWaveChannelOne
	Make /T/N=(numpnts(XWave)) IntegrationWaveChannelTwo

	// waves to store the numerical outcome of the integration,
	// the ratio of the integration, and the mean value
	// of the line at the integration interval
	Make /N=(numpnts(XWave)) IntegrationValChannelOne
	Make /N=(numpnts(XWave)) IntegrationValChannelTwo

	Make /N=(numpnts(XWave)) RatioChannelOneChannelTwo
	Make /N=(numpnts(XWave)) RatioChannelTwoChannelOne

	Make /N=(numpnts(XWave)) MeanValueChannelOne
	Make /N=(numpnts(XWave)) MeanValueChannelTwo

	ComparisonGet3by3BeamAverageAll()
	BaseIntRatMVComp()

	ExpCalcCompVal()

	
End


Function ExpCalcCompVal()
	// data to store
	Wave XWave = root:CompareBackgroundChannel:XW
	Wave YWave = root:CompareBackgroundChannel:YW
	wave IntegrationStart = root:CompareBackgroundChannel:Event_Start
	wave IntegrationEnd = root:CompareBackgroundChannel:Event_End
	wave IntegrationValOne = root:CompareBackgroundChannel:IntegrationValChannelOne
	wave IntegrationValTwo = root:CompareBackgroundChannel:IntegrationValChannelTwo
	wave RatioOneTwo = root:CompareBackgroundChannel:RatioChannelOneChannelTwo
	wave RatioTwoOne = root:CompareBackgroundChannel:RatioChannelTwoChannelOne
	wave MeanOne = root:CompareBackgroundChannel:MeanValueChannelOne
	wave MeanTwo = root:CompareBackgroundChannel:MeanValueChannelTwo


	// Size of the waves (data plus 1 for the header)
	
	variable CompWaveSize = numpnts(XWave) + 1

	// list of the waves to store in the text file
	String/G ComparisonWaveNames = "Comp_X;Comp_Y;Comp_Start;Comp_End;Comp_Ch1Integration;Comp_Ch2Integration;Comp_Ch1Ch2Ratio;Comp_Ch2Ch1Ratio;Comp_Ch1Mean;Comp_Ch2Mean"

	// create text waves to store 

	Make /O/T/N=(CompWaveSize) Comp_X
	Make /O/T/N=(CompWaveSize) Comp_Y
	Make /O/T/N=(CompWaveSize) Comp_Start
	Make /O/T/N=(CompWaveSize) Comp_End
	Make /O/T/N=(CompWaveSize) Comp_Ch1Integration
	Make /O/T/N=(CompWaveSize) Comp_Ch2Integration
	Make /O/T/N=(CompWaveSize) Comp_Ch1Ch2Ratio
	Make /O/T/N=(CompWaveSize) Comp_Ch2Ch1Ratio
	Make /O/T/N=(CompWaveSize) Comp_Ch1Mean
	Make /O/T/N=(CompWaveSize) Comp_Ch2Mean
	
	// put in headers
	Comp_X[0] = "X"
	Comp_Y[0] = "Y"
	Comp_Start[0] = "Event Start"
	Comp_End[0] = "Event End"
	Comp_Ch1Integration[0] = "Channel 1 Integration"
	Comp_Ch2Integration[0] = "Channel 2 Integration"
	Comp_Ch1Ch2Ratio[0] = "Channel 1 / Channel 2 Integration"
	Comp_Ch2Ch1Ratio[0] = "Channel 2 / Channel 1 Integration"
	Comp_Ch1Mean[0] = "Channel 1 Mean Value"
	Comp_Ch2Mean[0] = "Channel 2 Mean Value"

	// fill in the text waves with the numerical values
	String NumericalValue
	variable Itor
	for(Itor = 1; Itor <= numpnts(XWave); Itor++)
		sprintf NumericalValue, "%.10f", XWave[itor - 1]
		Comp_X[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", YWave[itor - 1]
		Comp_Y[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", IntegrationStart[itor - 1]
		Comp_Start[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", IntegrationEnd[itor - 1]
		Comp_End[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", IntegrationValOne[itor - 1]
		Comp_Ch1Integration[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", IntegrationValTwo[itor - 1]
		Comp_Ch2Integration[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", RatioOneTwo[itor - 1]
		Comp_Ch1Ch2Ratio[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", RatioTwoOne[itor - 1]
		Comp_Ch2Ch1Ratio[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", MeanOne[itor - 1]
		Comp_Ch1Mean[Itor] = NumericalValue
		sprintf NumericalValue, "%.10f", MeanTwo[itor - 1]
		Comp_Ch2Mean[Itor] = NumericalValue
	Endfor

	// now save the data
	Save/J/B ComparisonWaveNames

End


	




Function BaseIntRatMVComp()

	wave/T Ch1WaveNameList = root:CompareBackgroundChannel:ComparisonWaveChannelOne
	wave/T Ch2WaveNameList = root:CompareBackgroundChannel:ComparisonWaveChannelTwo
	wave/T Ch2IntegrationNameList = root:CompareBackgroundChannel:IntegrationWaveChannelOne
	wave/T Ch1IntegrationNameList = root:CompareBackgroundChannel:IntegrationWaveChannelOne
	wave IntegrationValOne = root:CompareBackgroundChannel:IntegrationValChannelOne
	wave IntegrationValTwo = root:CompareBackgroundChannel:IntegrationValChannelTwo
	wave RatioOneTwo = root:CompareBackgroundChannel:RatioChannelOneChannelTwo
	wave RatioTwoOne = root:CompareBackgroundChannel:RatioChannelTwoChannelOne
	wave MeanOne = root:CompareBackgroundChannel:MeanValueChannelOne
	wave MeanTwo = root:CompareBackgroundChannel:MeanValueChannelTwo

	wave IntegrationStart = root:CompareBackgroundChannel:Event_Start
	wave IntegrationEnd = root:CompareBackgroundChannel:Event_End

	variable Point_Itor

	For(Point_Itor = 0; Point_Itor < dimsize(Ch1WaveNameList,0); Point_Itor++)
		//--------------------------------------------------------
		// first calculate the baseline and replace in the file
		//--------------------------------------------------------

		wave Ch1Plot = $Ch1WaveNameList[Point_Itor]

		// get the baseline channel one
		duplicate /O Ch1Plot Ch1Baseline
		smooth /M=0 (floor(numPnts(Ch1Baseline)/3)), Ch1Baseline
		//subtract it
		Ch1Plot = Ch1Plot - Ch1Baseline

		wave Ch2Plot = $Ch2WaveNameList[Point_Itor]

		// get the baseline channel two
		duplicate /O Ch2Plot Ch2Baseline
		smooth /M=0 (floor(numPnts(Ch2Baseline)/3)), Ch2Baseline
		//subtract it
		Ch2Plot = Ch2Plot - Ch2Baseline

		//----------------------------------------------------------------------------
		// second calculate the integration of the region determined by the curve fit
		//----------------------------------------------------------------------------

		// first duplicate the original channel one and crop it at the same time
		Duplicate /O/R=[floor(IntegrationStart[Point_Itor]),floor(IntegrationEnd[Point_Itor])] $Ch1WaveNameList[Point_Itor] CroppedCh1
		// integrate it
		Integrate CroppedCh1 /D=IntegrateReturnValueCh1
		//input the final value of the integration wave
		IntegrationValOne[Point_Itor] = IntegrateReturnValueCh1[numpnts(IntegrateReturnValueCh1) - 1]
		//set the X scale to start at the integration start point for the cropped wave
		//useful for plotting
		SetScale/P x IntegrationStart[Point_Itor],1,"", IntegrateReturnValueCh1
		// duplicate the wave and store the name
		String IntegralWaveCh1 = Ch1WaveNameList[Point_Itor] + "Integral" 
		Duplicate /O IntegrateReturnValueCh1 $IntegralWaveCh1
		Ch1IntegrationNameList[Point_Itor] = IntegralWaveCh1

		// then duplicate channel two and crop it at the same time
		Duplicate /O/R=[floor(IntegrationStart[Point_Itor]),floor(IntegrationEnd[Point_Itor])] $Ch2WaveNameList[Point_Itor] CroppedCh2
		// integrate it
		Integrate CroppedCh2 /D=IntegrateReturnValueCh2
		//input the final value of the integration wave
		IntegrationValTwo[Point_Itor] = IntegrateReturnValueCh2[numpnts(IntegrateReturnValueCh2) - 1]
		//set the X scale to start at the integration start point for the cropped wave
		//useful for plotting
		SetScale/P x IntegrationStart[Point_Itor],1,"", IntegrateReturnValueCh2
		// duplicate the wave and store the name
		String IntegralWaveCh2 = Ch2WaveNameList[Point_Itor] + "Integral" 
		Duplicate /O IntegrateReturnValueCh2 $IntegralWaveCh2
		Ch2IntegrationNameList[Point_Itor] = IntegralWaveCh2

		//--------------------------------------------------------------------
		// Third calculate the mean over the cropped region of the curve fit
		//--------------------------------------------------------------------

		MeanOne[Point_Itor] = mean(CroppedCh1)
		MeanTwo[Point_Itor] = mean(CroppedCh2)

		//----------------------------------------------
		// Fourth calculate the ratios of the integrals
		//----------------------------------------------

		RatioOneTwo[Point_Itor] = IntegrationValOne[Point_Itor] / IntegrationValTwo[Point_Itor]
		RatioTwoOne[Point_Itor] = IntegrationValTwo[Point_Itor] / IntegrationValOne[Point_Itor]

	endfor


End


Function ComparisonGet3by3BeamAverageAll()
	
	svar ChannelOneRef = root:CompareBackgroundChannel:ChannelOneCompareReference	
	svar ChannelTwoRef = root:CompareBackgroundChannel:ChannelTwoCompareReference	
	wave/T Ch1WaveNameList = root:CompareBackgroundChannel:ComparisonWaveChannelOne
	wave/T Ch2WaveNameList = root:CompareBackgroundChannel:ComparisonWaveChannelTwo
	wave X_Points = root:CompareBackgroundChannel:XW
	wave Y_Points = root:CompareBackgroundChannel:YW


	// make a place to store the 3x3 averages waves
	variable Ch1_Frames = dimsize($ChannelOneRef, 2)
	Make /O/N=(Ch1_Frames) Ch1_Average	
	Make /O/N=(Ch1_Frames) Ch2_Average


	variable Point_Itor, X_Itor, Y_Itor
	For(Point_Itor = 0; Point_Itor < numpnts(X_Points); Point_Itor++)

		// sum all of the waves in a three by three grid for channel one
		Ch1_Average = 0
		for(X_Itor = -1; X_Itor < 2; X_Itor++)
			for(Y_Itor = -1; Y_Itor < 2; Y_Itor++)
				ImageTransform/beam={(X_Points[Point_Itor] + X_Itor), (Y_Points[Point_Itor] + Y_Itor)} getbeam $ChannelOneRef
				wave W_Beam
				Ch1_Average = Ch1_Average + W_Beam
			Endfor
		Endfor
		Ch1_Average = Ch1_Average / 9.0	

		// sum all of the waves in a three by three grid for channel two
		Ch2_Average = 0
		for(X_Itor = -1; X_Itor < 2; X_Itor++)
			for(Y_Itor = -1; Y_Itor < 2; Y_Itor++)
				ImageTransform/beam={(X_Points[Point_Itor] + X_Itor), (Y_Points[Point_Itor] + Y_Itor)} getbeam $ChannelTwoRef
				wave W_Beam
				Ch2_Average = Ch2_Average + W_Beam
			Endfor
		Endfor
		Ch2_Average = Ch2_Average / 9.0	

		// duplicate and store the wave names in a reference array

		String Ch1WaveName = "X"+num2str(X_Points[Point_Itor])+"Y"+num2str(Y_Points[Point_Itor])+"Ch1"
		String Ch2WaveName = "X"+num2str(X_Points[Point_Itor])+"Y"+num2str(Y_Points[Point_Itor])+"Ch2"
		Duplicate /O Ch1_Average $Ch1WaveName
		Duplicate /O Ch2_Average $Ch2WaveName
		Ch1WaveNameList[Point_Itor] = Ch1WaveName
		Ch2WaveNameList[Point_Itor] = Ch2WaveName

	Endfor
End


//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------




























































































Function GenerateScatterPlot(positions, xSize, ySize, pixelSize)
	wave positions
	variable xSize, ySize, pixelSize
	
	Assert(WaveExists(positions))
	Assert((xSize > 0) && (ySize > 0))
	
	variable nPositions = DimSize(positions, 0)
	variable nImages = positions[nPositions - 1][0] - positions[0][0] + 1
	variable pixelSizeUM = pixelSize / 1000
	
	Make /D/O/N=(nPositions) Scatter_positions_X, Scatter_positions_Y, Scatter_positions_Color
	wave Scatter_positions_X
	wave Scatter_positions_Y
	wave Scatter_positions_Color
	
	// get the indices of the columns containing the x and y coordinates
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1))
		Abort "The positions passed to GenerateScatterPlot() do not appear to contain any (x,y) information"
	endif
	
	Scatter_positions_X = positions[p][xCol]
	Scatter_positions_Y = positions[p][yCol]
	Scatter_positions_Color = positions[p][0]
	
	
	if (pixelSize != 0)
		Scatter_positions_X *= pixelSizeUM
		Scatter_positions_Y *= pixelSizeUM
		SetScale /P x, 0, 1, "um", Scatter_positions_X
		SetScale /P x, 0, 1, "um", Scatter_positions_Y
	else
		SetScale /P x, 0, 1, "pixel", Scatter_positions_X
		SetScale /P x, 0, 1, "pixel", Scatter_positions_Y
	endif
End

Function Generate3DScatterPlot(positions, xSize, ySize, pixelSize, colorScaleString)
	wave positions
	variable xSize, ySize, pixelSize
	string colorScaleString
	
	Assert(WaveExists(positions))
	Assert((xSize > 0) && (ySize > 0))
	
	variable nPositions = DimSize(positions, 0)
	variable nImages = positions[nPositions - 1][0] - positions[0][0] + 1
	
	variable pixelSizeUM = pixelSize / 1000
	
	Make /D/O/N=(nPositions) Scatter_positions_X
	Make /D/O/N=(nPositions) Scatter_positions_Y
	Make /D/O/N=(nPositions) Scatter_positions_Z
	wave Scatter_positions_X
	wave Scatter_positions_Y
	wave Scatter_positions_Z
	
	// get the indices of the columns containing the x and y coordinates
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1) || (zCol == -1))
		Abort "The positions passed to GenerateScatterPlot() do not appear to contain any (x,y,z) information"
	endif
	
	Scatter_positions_X = positions[p][xCol]
	Scatter_positions_Y = positions[p][yCol]
	Scatter_positions_Z = positions[p][zCol]
	
	if (pixelSize != 0)
		Scatter_positions_X *= pixelSizeUM
		Scatter_positions_Y *= pixelSizeUM
		Scatter_positions_Z *= pixelSizeUM
		SetScale /P x, 0, 1, "um", Scatter_positions_X
		SetScale /P x, 0, 1, "um", Scatter_positions_Y
		SetScale /P x, 0, 1, "um", Scatter_positions_Z
	else
		SetScale /P x, 0, 1, "pixel", Scatter_positions_X
		SetScale /P x, 0, 1, "pixel", Scatter_positions_Y
		SetScale /P x, 0, 1, "pixel", Scatter_positions_Z
	endif
	
	// determine range of z values
//	MatrixOP /FREE W_zPos = col(positions, zCol)
//	WaveStats /Q/M=1 W_zPos
//	variable minZ = V_min
//	variable maxZ = V_max
//	ColorTab2Wave $colorScaleString
//	wave M_Colors
//	Make /O/N=(nPositions, 3) /W/U Scatter_positions_color
//	wave Scatter_positions_color = Scatter_positions_color
//	variable nColors = DimSize(M_Colors, 0)
//	variable colorDelta = (maxZ - minZ) / nColors
//	Scatter_positions_color = M_Colors[floor((positions[p][zCol] - minZ) / (max(1e-6, maxZ - minZ + 1e-6)) * nColors)][q]
//	KillWaves /Z M_Colors
End

ThreadSafe Function /WAVE MakeAccumulatedImage(pos, xSize, ySize, binSize, [pixelSize])
	wave pos
	variable xSize, ySize, binSize, pixelSize
	
	variable nPoints = DimSize(pos, 0)
	variable xWidth = ceil(xSize / binSize)
	variable yWidth = ceil(ySize / binSize)
	
	variable px, py, i, xCol, yCol, zCol
	
	// get the indices of the columns containing the x and y coordinates
	getColumnsForEmitterPositions(pos, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1))
		Print "The positions passed to MakeAccumulatedImage() do not appear to contain any (x,y) information"
		return $""
	endif
	
	Make /D/FREE/N=(xWidth, yWidth) M_AccumulatedImage
	if (pixelSize == 0)	// no pixel size specified
		SetScale /I x, 0, xSize - 1, "pixel", M_AccumulatedImage
		SetScale /I y, 0, ySize - 1, "pixel", M_AccumulatedImage
	else		// pixel size specified, rescale the nm pixel size into um
		SetScale /I x, 0, (xSize - 1) * pixelSize / 1000, "um", M_AccumulatedImage
		SetScale /I y, 0, (ySize - 1) * pixelSize / 1000, "um", M_AccumulatedImage
	endif
	
	M_AccumulatedImage = 0
	
	for (i = 0; i < nPoints; i+=1)
		px = floor(pos[i][xCol] / binSize)
		py = floor(pos[i][yCol] / binSize)
		
		if ((px < 0) || (px >= xWidth) || (py < 0) || (py >= yWidth))
			continue
		endif
		
		M_AccumulatedImage[px][py] += 1
	endfor
	
	return M_AccumulatedImage
End





Function /WAVE ParticleTrackToScatterPairs(W_TracksWave)
	wave /WAVE W_TracksWave
	
	variable nTracks = DimSize(W_TracksWave, 0)
	variable nPositions = NumPositionsInTrack(W_TracksWave)
	
	Make /FREE/N=(nPositions + nTracks, 3) M_ScatterCoordinates
	variable i, xCol, yCol, zCol, offset = 0, nPosInTrack
	for (i = 0; i < nTracks; i+=1)
		wave theseTracks = W_TracksWave[i]
		nPosInTrack = DimSize(theseTracks, 0)
		GetColumnsForEmitterPositions(theseTracks, xCol, yCol, zCol)
		M_ScatterCoordinates[offset, offset + nPosInTrack - 1][0] = theseTracks[p - offset][xCol]
		M_ScatterCoordinates[offset, offset + nPosInTrack - 1][1] = theseTracks[p - offset][yCol]
		M_ScatterCoordinates[offset, offset + nPosInTrack - 1][2] = i
		M_ScatterCoordinates[offset + nPosInTrack] = NaN
		offset += nPosInTrack + 1
	endfor
	
	return M_ScatterCoordinates
End

Function /WAVE GetColorsForScatterCoordinates(M_ScatterCoordinates, M_Colors)
	wave M_ScatterCoordinates
	wave M_Colors
	
	variable nColors = DimSize(M_Colors, 0)
	WaveStats /Q/M=1 M_ScatterCoordinates
	variable nTracks = V_numNaNs
	variable nCoordinates = DimSize(M_ScatterCoordinates, 0)
	
	Make /FREE/W/U/N=(nCoordinates, 3) M_TrackColors
	variable i, trackIndex = 0
	for (i = 0; i < nCoordinates; i+=1)
		if (NumType(M_ScatterCoordinates[i][0]) == 2)
			M_TrackColors[i][] = NaN
			trackIndex += 1
		else
			M_TrackColors[i][] = M_Colors[floor(nColors / nTracks * trackIndex)][q]
		endif
	endfor
	
	return M_TrackColors
End

Function NumPositionsInTrack(W_TracksWave)
	wave /WAVE W_TracksWave
	
	variable nTracks = DimSize(W_TracksWave, 0)
	variable i, nPositionsInTrack = 0
	for (i = 0; i < nTracks; i+=1)
		wave theseTracks = W_TracksWave[i]
		nPositionsInTrack += DimSize(theseTracks, 0)
	endfor
	
	return nPositionsInTrack
End

Function ConsolidateSameEmitters(positions, posOutputName, maxPositionDifference, maxBlinking, [ minLocalizationsForError, maxLocalizationsForError])
	wave positions
	string posOutputName
	variable maxPositionDifference, maxBlinking
	variable minLocalizationsForError, maxLocalizationsForError
	
	Assert(WaveExists(positions))
	Assert(maxPositionDifference >= 0)
	Assert(maxBlinking >= 0)
	
	// Go over all the positions in the positions wave and look for emitters that show up in
	// different frames but might actually be the same spot. The two emitters have to be close
	// enough to be counted as the same molecule (maxPositionDifference). To allow for blinking
	// a molecule can disappear from the images for a number of frames up to maxBlinking
	// and still be accepted as the same emitter
	
	// the algorithm works as follows: every new position is temporarily added to a 'waiting list'
	// that contains the index into positions where it is listed, its x and y position, and a counter
	// that is initially set to maxBlinking
	
	// when moving on to a new frame every position on the waiting list is checked to see if a point
	// within a distance of maxPositionDifference is present in the new frame. If so then the number of
	// frames associated with the position on the waiting list is incremented by one, the counter is reset
	// to maxBlinking, and the matching point is deleted from the positions wave.
	// If not then the counter is decremented. If the counter becomes negative then the point is deleted .
	
	// the format of the waiting list is index of first appearance / x position / y position / counter
	
	// in its original form the algorithm actively deleted points that were the same as a previous emitter
	// during the processing. But this introduces the overhead of copying the wave elements each time.
	// therefore in its current form the algorithm does not delete points directly, but works on a duplicate
	// of the positions where the entries to be deleted are marked by a negative frame number
	
	variable minLocForError, maxLocForError
	if (ParamIsDefault(minLocalizationsForError))
		minLocForError = 2
	else
		minLocForError = minLocalizationsForError
	endif
	if (ParamIsDefault(maxLocalizationsForError))
		maxLocForError = inf
	else
		maxLocForError = maxLocalizationsForError
	endif
	if (minLocForError > maxLocForError)
		Abort "Invalid limits for localization error"
	endif
	
	variable nFramesPresentColumn, xCol, yCol, zCol, intensityColumn, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol
	GetColumnForNFramesPresent(positions, nFramesPresentColumn)	// handle different types of positions
	GetColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	GetColumnForIntegratedIntensity(positions, intensityColumn)
	GetColumnsForLocalizationError(positions, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol)
	
	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(positions, xSize, ySize, pixelSize)
	
	// check if these positions have not been consolidated before
	// if so then the results are unpredictable
	ImageTransform /G=(nFramesPresentColumn) sumCol, positions
	if (V_value != DimSize(positions, 0))
		DoAlert 1, "These positions look as if they have already been combined, or were created using an old version of this program. In the latter case, would you like me to correct this and continue?"
		if (V_flag == 1)	// yes
			positions[][nFramesPresentColumn] = 1
		else
			return 0
		endif
	endif
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer

	// group the emitters
	wave /WAVE W_GroupedEmitters = GroupEmitters(positions, CombinedPosEstimator_Mean, maxPositionDifference, maxBlinking)
	variable nGroups = DimSize(W_GroupedEmitters, 0)
	variable i, localizationStdDev, trackLength
	
	// allocate an output positions wave with the right number of entries
	Make /O/N=(nGroups, DimSize(positions, 1)) /D $posOutputName
	wave consolidatedPositions = $posOutputName
	
	// make room for some other statistics
	Make /N=(nGroups) /D /FREE W_LocalizationErrors
	Make /N=(nGroups) /D /FREE W_TrackLength
	
	for (i = 0; i < nGroups; i+=1)
		LocalizerProgressWindow("Combining", i, nGroups)
		wave currentEmitters = W_GroupedEmitters[i]
		wave combinedEmitters = CombineEmittersIntoPosition(currentEmitters, localizationStdDev, trackLength)
		consolidatedPositions[i][] = combinedEmitters[0][q]
		W_LocalizationErrors[i] = localizationStdDev
		W_TrackLength[i] = trackLength
	endfor
	
	Make /N=0 /D /O root:Packages:Localizer:W_LocalizationErrorHist
	Histogram /B=4 W_LocalizationErrors, root:Packages:Localizer:W_LocalizationErrorHist
	// strategy for the track length histogram: no more than 256 bins. If no molecules survive longer than 255 frames,
	// limit the histogram to that range.
	variable nPointsInTrackLengthHistogram = min(1024, WaveMax(W_TrackLength) - 1)
	Make /N=(nPointsInTrackLengthHistogram) /D /O root:Packages:Localizer:W_TrackLengthHist
	SetScale /I x, 1, WaveMax(W_TrackLength), root:Packages:Localizer:W_TrackLengthHist
	Histogram /B=2 W_TrackLength, root:Packages:Localizer:W_TrackLengthHist
	
	// set the wave scaling of the statistics to match the pixel size if it was specified
	GetCCDDimensionsFromPositions(positions, xSize, ySize, pixelSize)
	if (pixelSize != 0)
		SetScale /P x, 0, DimDelta(root:Packages:Localizer:W_LocalizationErrorHist, 0) * pixelSize, root:Packages:Localizer:W_LocalizationErrorHist
	endif
	
	// add info on the calculation params to the wave note
	string waveNote = Note(positions) + "CONSOLIDATE MAX SHIFT:" + num2str(maxPositionDifference) + ";" + "CONSOLIDATE MAX BLINKING:" + num2str(maxBlinking) + ";"
	Note /K consolidatedPositions, waveNote

	// perform a quality check
	ImageTransform /G=(nFramesPresentColumn) sumCol, consolidatedPositions
	if (DimSize(positions, 0) != V_Value)
		Abort "The calculation appears to have gone wrong, the results should not be trusted"
	endif
	
	// display the histograms
	DoWindow /F LocalizationErrorViewer
	if (V_flag != 1)
		Display /K=1 /N=LocalizationErrorViewer root:Packages:Localizer:W_LocalizationErrorHist
		if (pixelSize == 0)
			Label Bottom, "Estimated localization error (pixels)"
		else
			Label Bottom, "Estimated localization error (nm)"
		endif
		Label Left, "Occurrence"
	endif
	
	DoWindow /F TrackLengthViewer
	if (V_flag != 1)
		Display /K=1 /N=TrackLengthViewer root:Packages:Localizer:W_TrackLengthHist
		Label Bottom, "Track length (frames)"
		Label Left, "Occurrence"
		AutoPositionWindow /R=LocalizationErrorViewer TrackLengthViewer
	endif
End

Function ParticleTracking(positions, posOutputName, maxPositionJump, maxBlinking, [ minTrackLength, maxTrackLength])
	wave positions
	string posOutputName
	variable maxPositionJump, maxBlinking
	variable minTrackLength, maxTrackLength
	
	Assert(WaveExists(positions))
	Assert(maxPositionJump >= 0)
	Assert(maxBlinking >= 0)
	
	// perform single particle tracking, with possible limitation on the length of an
	// acceptable track
	// This output of this function is a wave containing wave references
	// each wave reference contains a single track
	// the tracks are sorted into ascending frame of first occurrence
	
	// TODO: 
	// 1. the analysis ignores the possibility that more than one positions may be localized
	// within maxPositionJump. In this case, the emitter that the algorithm will choose is currently
	// undefined
	// 2. the analysis ignores the possibility that the particle may diffuse while it is not emitting
	
	variable localMinTrackLength, localMaxTrackLength
	if (ParamIsDefault(minTrackLength))
		localMinTrackLength = 2
	else
		localMinTrackLength = minTrackLength
	endif
	if (ParamIsDefault(maxTrackLength))
		localMaxTrackLength = inf
	else
		localMaxTrackLength = maxTrackLength
	endif
	if (localMinTrackLength > localMaxTrackLength)
		Abort "Invalid limits for localization error"
	endif
	
	variable nFramesPresentColumn, xCol, yCol, zCol, intensityColumn, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol
	GetColumnForNFramesPresent(positions, nFramesPresentColumn)	// handle different types of positions
	GetColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	GetColumnForIntegratedIntensity(positions, intensityColumn)
	GetColumnsForLocalizationError(positions, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol)
	
	// check if these positions have not been consolidated before
	// if so then there's no meaning in performing this
	ImageTransform /G=(nFramesPresentColumn) sumCol, positions
	if (V_value != DimSize(positions, 0))
		Abort "These positions have been consolidated. There is no point in continuing"
	endif
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer

	// group the emitters
	wave /WAVE W_GroupedEmitters = GroupEmitters(positions, CombinedPosEstimator_LastPos, maxPositionJump, maxBlinking)
	
	// filter the emitter tracks to make sure that they are within the acceptable bounds
	variable i
	variable nGroups = DimSize(W_GroupedEmitters, 0)
	variable nGroupsRemaining = nGroups
	// first see how many emitters will remain
	for (i = 0; i < nGroups; i+=1)
		wave currentEmitter = W_GroupedEmitters[i]
		if ((DimSize(currentEmitter, 0) < localMinTrackLength) || (DimSize(currentEmitter, 0) > localMaxTrackLength))
			nGroupsRemaining -= 1
		endif
	endfor
	
	Make /FREE /WAVE /N=(nGroupsRemaining) W_FilteredEmitters
	variable offset
	for (i = 0; i < nGroups; i+=1)
		wave currentEmitter = W_GroupedEmitters[i]
		if ((DimSize(currentEmitter, 0) < localMinTrackLength) || (DimSize(currentEmitter, 0) > localMaxTrackLength))
			continue
		endif
		W_FilteredEmitters[offset] = W_GroupedEmitters[i]
		offset += 1
	endfor
	
	// no more modification is needed, just make a copy of the positions and add a wavenote
	string waveNote = Note(positions) + "PARTICLETRACKING MAX SHIFT:" + num2str(maxPositionJump) + ";" + "PARTICLETRACKING MAX BLINKING:" + num2str(maxBlinking) + ";"
	Note /K W_FilteredEmitters, waveNote
	Duplicate /WAVE /O W_FilteredEmitters, $posOutputName
End

Function /WAVE GroupEmitters(positions, PositionEstimator, maxPositionDifference, maxBlinking)
	wave positions
	FUNCREF CombinedPosEstimator_Prototype PositionEstimator
	variable maxPositionDifference, maxBlinking
	
	// given a list of positions, group these into subsets that combine positions occurring in subsequent
	// frames and that are not shifted in position by more than maxPositionDifference
	// also allow the emitters to blink for maxBlinking number of frames before being dropped from consideration
	
	variable nFramesPresentColumn, xCol, yCol, zCol, intensityColumn, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol
	GetColumnForNFramesPresent(positions, nFramesPresentColumn)	// handle different types of positions
	GetColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	GetColumnForIntegratedIntensity(positions, intensityColumn)
	GetColumnsForLocalizationError(positions, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol)
	
	// maintain a copy of the string containing the format specification for these positions
	// each grouped set of emitters needs this string
	string positionsFormatStr = "LOCALIZATION METHOD:" + StringByKey("LOCALIZATION METHOD", Note(positions)) + ";"
	
	// Make a copy of the positions
	Duplicate /FREE /O positions, M_TempConsolidatedPositions
	
	// set up some waves that will provide temporary storage
	Make /N=0 /WAVE /FREE W_WaitingList
	// M_WaitingList is a wave containing wave references. Each of these references points to
	// a 2D wave, that contains those localizations that are thought to arise from the same particle.
	// each entry is a full localization, so contains the same number of columns as the positions wave
	Make /N=0 /I /FREE W_RemainingFrames
	// contains the number of frames this emiter can be off while still being counted as active
	
	Make /N=0 /FREE /D W_ClosestEmitterDistances
	Make /N=0 /I /FREE /D W_ClosestEmitterIndex
	// these waves will be used internally by the algorithm
	// to allow the globally closest emitter - previous location
	// to be found
	// W_ClosestEmitterDistances will contain the distance from the corresponding entry in the waiting
	// list to the closest emitter in the current frame. If there is no such emitter then it will be set to NaN
	// W_ClosestEmitterIndex contains the index of the closest emitter in this frame to the corresponding
	// emitter in the waiting list. If there is no such emitter then it is set to -1
	
	Make /N=0 /WAVE /FREE W_GroupedEmitters
	// the result of this function
	// each entry is a reference to a wave containing one or more positions
	// that are considered to belong together
	// every position shows up somewhere, and only once, in this wave
	// though it may be present in a one-row wave or together with a bunch of others
	
	// check if the positions have been consolidated. This should not be the case since otherwise
	// there is no real reason to run this function. This is also so ExtractPositionsInFrame() does
	// not have to perform this check each time.
	ImageTransform /G=(nFramesPresentColumn) sumCol, positions
	variable positionsHaveBeenConsolidated = V_value != DimSize(positions, 0)
	if (positionsHaveBeenConsolidated != 0)
		Abort "GroupEmitters() being called on consolidated positions"
	endif
	
	variable n, err, i
	variable startingFrame = positions[0][0]	// warning: magic numbers
	variable endingFrame = positions[DimSize(positions,0) - 1][0]
	variable nFramesInCalculation = endingFrame - startingFrame + 1
	
	variable lastX, lastY, nPointsInThisTrajectory, localizationStdDev, trackLength, minIndex
	variable nearestX, nearestY, index, distance, offset, nPointsToBeDeleted = 0
	variable currentXPos, currentYPos, currentNFrames, currentXUncertainty, currentYUncertainty
	variable previouslyMatchedEmitter
	string trackWaveName
	
	for (n = startingFrame; n <= endingFrame; n+=1)
		// the calculation can be slow occasionally, so add a progress bar
		// but only if there is a sufficiently large number of frames
		if (nFramesInCalculation > 10)
			LocalizerProgressWindow("Grouping", n - startingFrame, endingFrame - startingFrame + 1)
		endif
		
		// get all emitters localized in this frame
		wave /WAVE M_Extracted = ExtractPositionsInFrame(positions, n, positionsHaveBeenConsolidated=positionsHaveBeenConsolidated)
		wave M_emittersInCurrentFrame = M_Extracted[0]
		wave W_emittersInCurrentFrameIndices = M_Extracted[1]
		
		// for every emitter in the waiting list we need to check if there is a matching emitter in the 
		// current frame. One of the side issues with that is that more than one of the emitters in the
		// waiting list can match the same one in the current frame. So iterate over this procedure
		// in such a way that the closest emitter is taken at each point
		Redimension /N=(DimSize(W_WaitingList, 0)) W_ClosestEmitterDistances, W_ClosestEmitterIndex
		W_ClosestEmitterDistances = inf
		W_ClosestEmitterIndex = -1
		previouslyMatchedEmitter = NaN
		do
			
			for (i = 0; i < DimSize(W_WaitingList, 0); i+=1)
				if (NumType(W_ClosestEmitterDistances[i]) == 2)
					// every emitter that has been matched already, or that has no corresponding match,
					// will have its entry in W_ClosestEmitterDistances set to NaN
					// so these points need to be skipped
					continue
				endif
				
				if ((W_ClosestEmitterIndex[i] != previouslyMatchedEmitter) && (NumType(previouslyMatchedEmitter) != 2))
					// in all iterations of the do-while loop following the first there is no need to recalculate any of the matching
					// emitters, except when this emitter matched to the same one as the one matched in the previous iteration
					
					// however, if the matching index of the current emitter is larger than the previously matched one, then we need to
					// decrement it by one. This compensates for the deletion of the emitter from the extracted positions list upon a
					// successful match
					if (W_ClosestEmitterIndex[i] > previouslyMatchedEmitter)
						W_ClosestEmitterIndex[i] -= 1
					endif
					
					continue
				endif
				
				// get the point in the current frame that is the closest to the last observation of this emitter
				// in the waiting list
				wave currentEmitter = W_WaitingList[i]
				
				// the measure of the last position of the emitter depends on what we want to do with the results
				// so it is parameterized using PositionEstimator. This function takes care of filling in
				// lastX and lastY
				PositionEstimator(currentEmitter, xCol, yCol, lastX, lastY)
				
				err = returnPointNearestFitPositions(M_emittersInCurrentFrame, lastX, lastY, nearestX, nearestY, index)
				if (err != 0)
					// some error looking for the nearest point
					// there are probably no more emitters in M_emittersInCurrentFrame
					W_ClosestEmitterDistances[i] = NaN
					W_ClosestEmitterIndex[i] = -1
					continue
				endif
				
				distance = sqrt((nearestX - lastX)^2 + (nearestY - lastY)^2)
				
				if (distance > maxPositionDifference)
					// the point on the waiting list is not present in the current frame
					W_ClosestEmitterDistances[i] = NaN
					W_ClosestEmitterIndex[i] = -1
				
				else // the point on the waiting list is also present in this frame
					// add the information to the intermediate waves
					W_ClosestEmitterDistances[i] = distance
					W_ClosestEmitterIndex[i] = index
				endif
			endfor
			
			// W_ClosestEmitterDistances now contains all distances between the waiting list and the closest matching positions
			// select the minimum distance, i.e. the closest possible match, and evaluate that
			// leave the other matching emitters to subsequent iterations of the do-while loop
			// unfortunately WaveStats throws an error if there are no points in W_ClosestEmitterDistances
			// (which occurs when the waiting list is empty)
			if (DimSize(W_WaitingList, 0) == 0)
				minIndex = -1
			else
				WaveStats /Q /M=1 W_ClosestEmitterDistances
				minIndex = V_minLoc
			endif
			if (minIndex == -1)
				// if W_ClosestEmitterDistances only contains NaN then WaveStats sets V_minLoc to -1
				// time to move on to the next frame
				break
			endif
			
			// if we're here then we have a new matching point
			// add it to the output
			wave currentEmitter = W_WaitingList[minIndex]
			nPointsInThisTrajectory = DimSize(currentEmitter, 0)
			index = W_ClosestEmitterIndex[minIndex]
			Assert(index != -1)
					
			Redimension /N=(nPointsInThisTrajectory + 1, -1) currentEmitter
			currentEmitter[nPointsInThisTrajectory][] = M_emittersInCurrentFrame[index][q]
			W_ClosestEmitterDistances[minIndex] = NaN
			previouslyMatchedEmitter = index
					
			// delete the matching point in the frame from the extracted positions list
			DeletePoints index, 1, M_emittersInCurrentFrame
			DeletePoints index, 1, W_emittersInCurrentFrameIndices
		while (1)
		
		// W_ClosestEmitterDistances must now contain only NaN's
		if (DimSize(W_ClosestEmitterDistances, 0) > 0)
			WaveStats /M=1 /Q W_ClosestEmitterDistances
			Assert(V_numNaNs == DimSize(W_ClosestEmitterDistances, 0))
		endif
		
		// update the blinking counters
		// for this the algorithm needs to know whether the emitters were successfully matched or not
		// this cannot be checked using W_ClosestEmitterDistances since it's NaN everywhere
		// instead make use of W_ClosestEmitterIndex: if it is -1 then there was no successful match
		for (i = 0; i < DimSize(W_WaitingList, 0); i+=1)
			if (W_ClosestEmitterIndex[i] == -1)
				W_RemainingFrames[i] -= 1
			else
				W_RemainingFrames[i] = maxBlinking
			endif
		endfor
		
		// any emitter in the waiting list that has a negative value in W_RemainingFrames is finished
		// the corresponding entry in currentEmitter represents a finished group of positions
		// add it as an entry to W_GroupedEmitters
		for (i = 0; i < DimSize(W_WaitingList, 0); i+=1)
			if (W_RemainingFrames[i] < 0)
				wave currentEmitter = W_WaitingList[i]
				Redimension /N=(DimSize(W_GroupedEmitters, 0) + 1) W_GroupedEmitters
				W_GroupedEmitters[DimSize(W_GroupedEmitters, 0) - 1] = currentEmitter
				DeletePoints i, 1, W_WaitingList
				DeletePoints i, 1, W_RemainingFrames
				i -= 1
			endif
		endfor
		
		// any emitters that remain in the current frame are new, add them to the waiting list
		offset = DimSize(W_WaitingList, 0)
		Redimension /N=(offset + DimSize(M_emittersInCurrentFrame, 0)) W_WaitingList
		Redimension /N=(offset + DimSize(M_emittersInCurrentFrame, 0)) W_RemainingFrames
		
		for (i = 0; i < Dimsize(M_emittersInCurrentFrame, 0); i+=1)
			Make /O/D/FREE /N=(1, DimSize(positions, 1)) M_NewEmitter
			// downstream routines (such as when combining these emitters for consolidation
			// must know what positions format this wave is in
			Note M_NewEmitter, positionsFormatStr
			M_NewEmitter = positions[W_emittersInCurrentFrameIndices[i]][q]
			W_WaitingList[i + offset] = M_NewEmitter
			W_RemainingFrames[offset + i] = maxBlinking
		endfor
		
	endfor
	
	// if we're here than we've looped over all the frames in the movie
	// it's possible that some emitters were 'on' throughout the movie
	// in that case they are still in W_WaitingList and need to be added
	// without regard to their blinking allowance
	for (i = 0; i < DimSize(W_WaitingList, 0); i+=1)
		wave currentEmitter = W_WaitingList[i]
		Redimension /N=(DimSize(W_GroupedEmitters, 0) + 1) W_GroupedEmitters
		W_GroupedEmitters[DimSize(W_GroupedEmitters, 0) - 1] = currentEmitter
		DeletePoints i, 1, W_WaitingList
		DeletePoints i, 1, W_RemainingFrames
		i -= 1
	endfor
	
	
	// no more emitters may be in the waiting list now
	Assert(DimSize(W_WaitingList, 0) == 0)
	
	// sort the grouped emitters such that the first frame they occur in
	// is in order
	// this is somewhat involved
	variable nEntries = DimSize(W_GroupedEmitters, 0)
	Make /O/N=(nEntries) /I/U /FREE W_FirstFrames
	Make /O/N=(nEntries) /I/U /FREE W_Indices
	W_Indices = p
	for (i = 0; i < nEntries; i+=1)
		wave currentEmitter = W_GroupedEmitters[i]
		W_FirstFrames[i] = currentEmitter[0][0]
	endfor
	Sort W_FirstFrames, W_Indices
	Duplicate /O/FREE /WAVE W_GroupedEmitters, W_TempGroupedEmitters
	W_GroupedEmitters = W_TempGroupedEmitters[W_Indices[p]]
	
	// perform a quality check on the algorithm
	// the total number of positions must not have been changed
	variable nGroupedPositions = 0;
	for (i = 0; i < nEntries; i+=1)
		wave currentEmitter = W_GroupedEmitters[i]
		nGroupedPositions += DimSize(currentEmitter, 0)
	endfor
	
	if (DimSize(positions, 0) != nGroupedPositions)
		Abort "The calculation appears to have gone wrong, the results should not be trusted"
	endif
	
	return W_GroupedEmitters
End

Function CombinedPosEstimator_Prototype(combinedPositions, xCol, yCol, xPos, yPos)
	wave combinedPositions
	variable xCol, yCol
	variable &xPos
	variable &yPos
	// to decide whether emitters belong together, we need some kind of metric
	// for the actual position of an emitter already fitted over several frames
	// this metric depends on what we want to do: for simple consolidation it should
	// probably be the average of all previous positions,
	// while for particle tracking it should really be the last observed positions
	
	// the current solution is to parametrize these different behaviors by providing different functions that
	// can then be passed to GroupEmitters
	// this is the prototype, which should not be called
	Abort "CombinedPosEstimator_Prototype should not be called"
End

Function CombinedPosEstimator_Mean(combinedPositions, xCol, yCol, xPos, yPos)
	wave combinedPositions
	variable xCol, yCol
	variable &xPos
	variable &yPos
	// estimate the actual position of the emitter as the mean of all observations
	variable nPositions = DimSize(combinedPositions, 0)
	
	ImageTransform /G=(xCol) sumCol, combinedPositions
	xPos = V_Value / nPositions
	ImageTransform /G=(yCol) sumCol, combinedPositions
	yPos = V_Value / nPositions
End

Function CombinedPosEstimator_LastPos(combinedPositions, xCol, yCol, xPos, yPos)
	wave combinedPositions
	variable xCol, yCol
	variable &xPos
	variable &yPos
	// estimate the actual position of the emitter as the mean of all observations
	variable nPositions = DimSize(combinedPositions, 0)
	
	xPos = combinedPositions[nPositions - 1][xCol]
	yPos = combinedPositions[nPositions - 1][yCol]
End

Function /WAVE RemoveOutliers(removeFromThesePos, testAgainstThesePos, radius, nPosInRadius)
	wave removeFromThesePos, testAgainstThesePos
	variable radius, nPosInRadius
	
	variable nPositionsToTest = DimSize(removeFromThesePos, 0)
	variable nPositionsToTestAgainst = DimSize(testAgainstThesePos, 0)
	variable lookingAtSamePositions = WaveRefsEqual(removeFromThesePos, testAgainstThesePos)
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(testAgainstThesePos, xCol, yCol, zCol)
	
	MatrixOP /FREE W_minX = minVal(col(testAgainstThesePos, xCol))
	MatrixOP /FREE W_minY = minVal(col(testAgainstThesePos, yCol))
	MatrixOP /FREE W_maxX = maxVal(col(testAgainstThesePos, xCol))
	MatrixOP /FREE W_maxY = maxVal(col(testAgainstThesePos, yCol))
	variable minX = W_minX[0], maxX = W_maxX[0]
	variable minY = W_minY[0], maxY = W_maxY[0]
	
	// strategy: bin all positions into a 2D histogram. Then use this histogram to quickly calculate
	// the number of neighbors. We assume that a pixel size of radius / 10 is small enough.
	variable histogramPixelSize = radius / 10
	variable nHistogramRows = ceil((maxX - minX) / histogramPixelSize) + 1
	variable nHistogramCols = ceil((maxY - minY) / histogramPixelSize) + 1
	Make /FREE /N=(nHistogramRows, nHistogramCols) /I/U M_NeighborHistogram
	SetScale /P x, minX, histogramPixelSize, M_NeighborHistogram
	SetScale /P y, minY, histogramPixelSize, M_NeighborHistogram
	FastOP M_NeighborHistogram = 0
	
	variable i, pp, qq
	for (i = 0; i < nPositionsToTestAgainst; i+=1)
		pp = round((testAgainstThesePos[i][xCol] - minX) / histogramPixelSize)
		qq = round((testAgainstThesePos[i][yCol] - minY) / histogramPixelSize)
		M_NeighborHistogram[pp][qq] += 1
	endfor
	
	Duplicate /FREE removeFromThesePos, filteredPositions
	
	variable nThreads = min(ThreadProcessorCount, nPositionsToTest)
	variable tgID = ThreadGroupCreate(nThreads), threadStatus
	Make /FREE/N=1/D W_Progress = 0, W_Abort=0
	Duplicate /FREE removeFromThesePos, filteredPositions
	for (i = 0; i < nThreads; i+=1)
		ThreadStart tgID, i, FindNuclei_Worker(i, nThreads, filteredPositions, M_NeighborHistogram, radius, nPosInRadius, lookingAtSamePositions, W_Progress, W_Abort)
	endfor
	
	// now wait for the threads to finish and occasionally report progress
	for ( ; ; )
		threadStatus = ThreadGroupWait(tgID, 100)
		if (threadStatus == 0)
			// all threads finished
			break
		endif
		// progress notification
		if (nPositionsToTest > 1000)
			W_Abort[0] = LocalizerProgressWindow("Mapping", W_Progress[0], nPositionsToTest, allowAbort=0)
		endif
	endfor
	
	threadStatus = ThreadGroupWait(tgID, inf)
	threadStatus = ThreadGroupRelease(tgID)
	
	if (W_Abort[0] != 0)
		// user aborted
		return $""
	endif
	
	// now get rid of all entries that contain a NaN
	variable offset = 0
	for (i = 0; i < nPositionsToTest; i+=1)
		if (NumType(filteredPositions[i]) == 2)
			continue
		endif
		filteredPositions[offset][] = filteredPositions[i][q]
		offset += 1
	endfor
	
	Redimension /N=(offset, -1) filteredPositions
	
	return filteredPositions
End

ThreadSafe Function FindNuclei_Worker(threadIndex, nThreads, PositionsToFilter, M_NeighborHistogram, radius, nPosInRadius, comparingSamePositions, W_Progress, W_Abort)
	variable nThreads, threadIndex
	wave PositionsToFilter, M_NeighborHistogram
	variable radius, nPosInRadius, comparingSamePositions
	wave W_Progress, W_Abort
	
	variable nPositions = DimSize(PositionsToFilter, 0)
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(PositionsToFilter, xCol, yCol, zCol)
	
	variable i, nPositionsProcessed = 0, nNeighbors
	for (i = threadIndex; i < nPositions; i += nThreads)
		if (mod(nPositionsProcessed, 100) == 0)
			W_Progress[0] += nPositionsProcessed
			nPositionsProcessed = 0
			if (W_Abort[0] != 0)
				return 0
			endif
		endif
		nNeighbors = CountNeighbors(PositionsToFilter[i][xCol], PositionsToFilter[i][yCol], M_NeighborHistogram, radius)
		if (comparingSamePositions)
			nNeighbors -= 1		// because we have included the point itself
		endif
		if (nNeighbors < nPosInRadius)
			PositionsToFilter[i][] = NaN
		endif
		nPositionsProcessed += 1
	endfor
	
End

ThreadSafe Function CountNeighbors(xx, yy, M_NeighborHistogram, radius)
	variable xx, yy
	wave M_NeighborHistogram
	variable radius
	
	variable nRowsHistogram = DimSize(M_NeighborHistogram, 0)
	variable nColsHistogram = DimSize(M_NeighborHistogram, 1)
	variable minX = DimOffset(M_NeighborHistogram, 0)
	variable minY = DimOffset(M_NeighborHistogram, 1)
	variable histPixelSize = DimDelta(M_NeighborHistogram, 0)
	variable radiusSq = radius * radius
	
	variable startP = limit(floor((xx - radius - minX) / histPixelSize), 0, nRowsHistogram - 1)
	variable endP = limit(ceil((xx + radius - minX) / histPixelSize), 0, nRowsHistogram - 1)
	variable startQ = limit(floor((yy - radius - minY) / histPixelSize), 0, nColsHistogram - 1)
	variable endQ = limit(ceil((yy + radius - minY) / histPixelSize), 0, nColsHistogram - 1)
	
	variable xHist, yHist, pp, qq
	variable nPosInRadius = 0
	for (qq = startQ; qq <= endQ; qq+=1)
		for (pp = startP; pp <= endP; pp+=1)
			xHist = minX + pp * histPixelSize
			yHist = minY + qq * histPixelSize
			if ((xHist - xx)^2 + (yHist - yy)^2 <= radiusSq)
				nPosInRadius += M_NeighborHistogram[pp][qq]
			endif
		endfor
	endfor
	
	return nPosInRadius
End

Function AutoCleanRegistrationPos(M_CombinedEmitterPositions, nPosInNeighborhood)
	wave M_CombinedEmitterPositions
	variable nPosInNeighborhood
	
	// strategy: for each match in M_CombinedEmitterPositions, combare the angle and distance of the match with the averages of its
	// nPosInNeighborhood nearest neighbors. If the difference is larger than maxAngleDiff or maxLengthDiff the point is rejected.
	
	variable nMatches = DimSize(M_CombinedEmitterPositions, 0)
	
	if (nPosInNeighborhood >= nMatches - 1)
		return 0
	endif
	
	Make /FREE/N=(nMatches)/D W_MedianDiff = 0, W_DistanceDiff = 0
	variable nThreads = min(ThreadProcessorCount, nMatches)
	Make /FREE/N=(nThreads) W_Dummy
	MultiThread W_Dummy = AutoCleanRegistrationWorker(p, nThreads, M_CombinedEmitterPositions, nPosInNeighborhood, W_MedianDiff)
	variable threshold
	WaveStats /Q W_MedianDiff
	W_MedianDiff -= V_avg
	threshold = 2 * V_sdev
	
	variable offset = 0, i
	for (i = 0; i < nMatches; i+=1)
		if (abs(W_MedianDiff[i]) > threshold)
			continue
		endif
		M_CombinedEmitterPositions[offset][] = M_CombinedEmitterPositions[i][q]
		offset += 1
	endfor
	Redimension /N=(offset, -1) M_CombinedEmitterPositions
End

ThreadSafe Function AutoCleanRegistrationWorker(threadIndex, nThreads, M_CombinedEmitterPositions, nPosInNeighborhood, W_MedianDiff)
	variable threadIndex, nThreads
	wave M_CombinedEmitterPositions
	variable nPosInNeighborhood
	wave W_MedianDiff
	
	variable nMatches = DimSize(M_CombinedEmitterPositions, 0)
	variable nMatchesPerThread = round(nMatches / nThreads)
	variable thisThreadFirstMatch = limit(threadIndex * nMatchesPerThread, 0, nMatches - 1)
	variable thisThreadLastMatch = limit((threadIndex + 1) * nMatchesPerThread - 1, 0, nMatches - 1)
	
	variable i, j, avgDistance = 0, avgAngle, thisDistance, thisAngle, nearestDistance, nearestAngle, closestIndex, dx, dy
	variable /C complexVal
	Make /FREE/N=(nPosInNeighborhood)/D W_ClosestMatchIndices, W_DX, W_DY
	for (i = thisThreadFirstMatch; i <= thisThreadLastMatch; i+=1)
		dx = M_CombinedEmitterPositions[i][2] - M_CombinedEmitterPositions[i][0]
		dy = M_CombinedEmitterPositions[i][3] - M_CombinedEmitterPositions[i][1]
		FindNClosestMatches(i, M_CombinedEmitterPositions, W_ClosestMatchIndices)
		W_DX = M_CombinedEmitterPositions[W_ClosestMatchIndices[p]][2] - M_CombinedEmitterPositions[W_ClosestMatchIndices[p]][0]
		W_DY = M_CombinedEmitterPositions[W_ClosestMatchIndices[p]][3] - M_CombinedEmitterPositions[W_ClosestMatchIndices[p]][1]
		W_MedianDiff[i] = sqrt((dx - StatsMedian(W_DX))^2 + (dy - StatsMedian(W_DY))^2)
	endfor
	
	return 0
End

ThreadSafe Function DisplacementToAngleAndDistance(x1, y1, x2, y2, angle, distance)
	variable x1, y1, x2, y2
	variable &angle, &distance
	
	variable /C complexVal = r2polar(cmplx(x2 - x1, y2 - y1))
	distance = real(complexVal)
	angle = imag(complexVal)
End

ThreadSafe Function FindNClosestMatches(closestToThisIndex, M_CombinedEmitterPositions, W_ClosestMatchIndices)
	variable closestToThisIndex
	wave M_CombinedEmitterPositions, W_ClosestMatchIndices
	
	variable nMatches = DimSize(M_CombinedEmitterPositions, 0)
	variable nClosestMatches = DimSize(W_ClosestMatchIndices, 0)
	variable xx = M_CombinedEmitterPositions[closestToThisIndex][0]
	variable yy = M_CombinedEmitterPositions[closestToThisIndex][1]
	
	Make /FREE/N=(nClosestMatches) /D W_ClosestSqDistances = inf
	
	variable i, sqDistance
	for (i = 0; i < nMatches; i += 1)
		if (i == closestToThisIndex)
			continue
		endif
		sqDistance = (M_CombinedEmitterPositions[i][0] - xx)^2 + (M_CombinedEmitterPositions[i][1] - yy)^2
		WaveStats /Q/M=1 W_ClosestSqDistances
		if (sqDistance < V_max)
			W_ClosestMatchIndices[V_maxRowLoc] = i
			W_ClosestSqDistances[V_maxRowLoc] = sqDistance
		endif
	endfor
End

Function /WAVE CreateRegistration_poly(M_EmitterPositions)
	wave M_EmitterPositions
	
	variable nEmitters = DimSize(M_EmitterPositions, 0)
	
	// start with x shift, then do y
	Make /FREE/N=(nEmitters)/D W_xPositionsChannel1 = M_EmitterPositions[p][0]
	Make /FREE/N=(nEmitters)/D W_yPositionsChannel1 = M_EmitterPositions[p][1]
	Make /FREE/N=(nEmitters)/D W_xPositionsChannel2 = M_EmitterPositions[p][2]
	Make /FREE/N=(nEmitters)/D W_yPositionsChannel2 = M_EmitterPositions[p][3]
	
	Make /FREE/N=(12)/D W_RegistrationMap
	
	MatrixOP /FREE W_ShiftInX = W_xPositionsChannel2 - W_xPositionsChannel1
	MatrixOP /FREE W_ShiftInY = W_yPositionsChannel2 - W_yPositionsChannel1
	
	CurveFit /N=1/W=2/Q poly2D 2, W_shiftInX /X=W_xPositionsChannel1 /Y=W_yPositionsChannel1
	wave W_Coef
	W_RegistrationMap[0,5] = W_Coef[p]
	
	CurveFit /N=1/W=2/Q poly2D 2, W_ShiftInY /X=W_xPositionsChannel1 /Y=W_yPositionsChannel1
	wave W_Coef
	W_RegistrationMap[6,11] = W_Coef[p - 6]
	
	return W_RegistrationMap
End

Function /WAVE CreateRegistration(M_EmitterPositions, nReferencesInLocalMean)
	wave M_EmitterPositions	// wave containing a row for every calibration emitter, using the format x (ch1) / y (ch1) / x (ch2) / y (ch2)
	variable nReferencesInLocalMean
	
	assert(DimSize(M_EmitterPositions, 0) >= nReferencesInLocalMean)
	assert(DimSize(M_EmitterPositions, 1) == 4)
	
	variable nEmitters = DimSize(M_EmitterPositions, 0)
	
	Make /FREE/N=(nEmitters)/D W_xPositions = M_EmitterPositions[p][0]
	Make /FREE/N=(nEmitters)/D W_yPositions = M_EmitterPositions[p][1]
	Make /FREE/N=(nEmitters)/I/U W_SortedIndices
	
	Make /FREE/N=(nReferencesInLocalMean)/D W_xValuesForFit, W_yValuesForFit, W_shiftInX, W_shiftInY
	
	Make /FREE/N=(nEmitters , 17)/D M_RegistrationMap
	// M_RegistrationMap has format
	// xLoc of ref / yLoc of ref / 6 coeffs for x poly / 6 coeffs for y poly / xLoc of centroid / yLoc of centroid / distance to furthest reference point for centroid of this poly
	// (see calculation of Rn below eq 21 in http://www.cs.wright.edu/~agoshtas/IVC88.pdf, except that distances are measured from the centroid)
	
	MatrixOP /FREE W_allShiftInX = col(M_EmitterPositions, 2) - col(M_EmitterPositions, 0)
	MatrixOP /FREE W_allShiftInY = col(M_EmitterPositions, 3) - col(M_EmitterPositions, 1)
	
	variable i, xLoc, yLoc, centroidXLoc, centroidYLoc, j
	for (i = 0; i < nEmitters; i+=1)
		xLoc = M_EmitterPositions[i][0]
		yLoc = M_EmitterPositions[i][1]
		centroidXLoc = 0
		centroidYLoc = 0
		
		M_RegistrationMap[i][0] = xLoc
		M_RegistrationMap[i][1] = yLoc
		
		// find the 6 nearest emitters (which will include the point under consideration)
		MatrixOP /FREE W_Distances = sqrt(magSqr(W_xPositions - xLoc) + magSqr(W_yPositions - yLoc))
		MakeIndex W_Distances, W_SortedIndices
		W_xValuesForFit = M_EmitterPositions[W_SortedIndices[p]][0]
		W_yValuesForFit = M_EmitterPositions[W_SortedIndices[p]][1]
		W_shiftInX = W_allShiftInX[W_SortedIndices[p]]
		W_shiftInY = W_allShiftInY[W_SortedIndices[p]]
		
		centroidXLoc = sum(W_xValuesForFit) / DimSize(W_xValuesForFit, 0)
		centroidYLoc = sum(W_yValuesForFit) / DimSize(W_yValuesForFit, 0)
		
		// perform a 2nd order 2D poly fit for the shift in x
		CurveFit /N=1/W=2/Q poly2D 2, W_shiftInX /X=W_xValuesForFit /Y=W_yValuesForFit
		wave W_Coef
		M_RegistrationMap[i][2,7] = W_Coef[q - 2]
		
		// perform a 2nd order 2D poly fit for the shift in y
		CurveFit /N=1/W=2/Q poly2D 2, W_shiftInY /X=W_xValuesForFit /Y=W_yValuesForFit
		// todo: handle failed fits
		wave W_Coef
		M_RegistrationMap[i][8,13] = W_Coef[q - 8]
		
		M_RegistrationMap[i][14] = centroidXLoc
		M_RegistrationMap[i][15] = centroidYLoc
		
		// calculate the  distance from the centroid to the farthest reference point
		MatrixOP /FREE W_MaxDistance = maxVal(sqrt(magSqr(W_xValuesForFit - centroidXLoc) + magSqr(W_yValuesForFit - centroidYLoc)))
		M_RegistrationMap[i][16] = W_MaxDistance[0]
	endfor
	
	return M_RegistrationMap
End

Function GuessSystematicShiftInPositions(pos1, pos2, maxJump, dx, dy)
	wave pos1, pos2
	variable maxJump, &dx, &dy
	// dx and dy start with initial guess, then get updated to estimated values
	
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(pos1, xCol, yCol, zCol)
	
	// For now this routine looks only at the first frame in each of the waves
	// and ignores the rest
	variable firstFrameIndex1 = pos1[0][0]
	variable firstFrameIndex2 = pos2[0][0]
	wave /WAVE extractedPositionsInFrame1= ExtractPositionsInFrame(pos1, firstFrameIndex1)
	wave /WAVE extractedPositionsInFrame2 = ExtractPositionsInFrame(pos2, firstFrameIndex2)
	
	wave extractedPositions1 = extractedPositionsInFrame1[0]
	wave extractedPositions2 = extractedPositionsInFrame2[0]
	
	variable nEmittersIn1 = DimSize(extractedPositions1, 0)
	variable nEmittersIn2 = DimSize(extractedPositions2, 0)
	Duplicate /O /FREE pos1, M_AppendedEmitters	// duplicate takes care of wave note
	Redimension /N=(nEmittersIn1 + nEmittersIn2, -1) M_AppendedEmitters
	M_AppendedEmitters[0, nEmittersIn1 - 1] = extractedPositions1[p][q]	// extractedPositions1 is never modified
	
	variable shiftX = dx, shiftY = dy
	variable estimatedShiftInThisRunX = dx, estimatedShiftInThisRunY = dy
	variable nIterations = 0, nMatchedEmitters, i
	do
		if (nIterations > 20)
			break
		endif
		
		extractedPositions2[][xCol] -= estimatedShiftInThisRunX
		extractedPositions2[][yCol] -= estimatedShiftInThisRunY
		
		// the grouping functions that already exist are designed for consolidation or particle tracking
		// so make a new positions wave, in which pos2 is appended to pos1
		M_AppendedEmitters[nEmittersIn1, ] = extractedPositions2[p - nEmittersIn1][q]
		M_AppendedEmitters[][0] = 0
		M_AppendedEmitters[nEmittersIn1, ][0] = 1
		
		wave /WAVE W_GroupedEmitters = GroupEmitters(M_AppendedEmitters, CombinedPosEstimator_LastPos, maxJump, 0)
		nMatchedEmitters = 0
		for (i = 0; i < DimSize(W_GroupedEmitters, 0); i+=1)
			wave thisWave = W_GroupedEmitters[i]
			if (DimSize(thisWave, 0) == 2)
				nMatchedEmitters += 1
			endif
		endfor
		
		Make /FREE/N=(nMatchedEmitters)/D W_xShift, W_yShift
		variable offset = 0
		for (i = 0; i < DimSize(W_GroupedEmitters, 0); i+=1)
			wave thisWave = W_GroupedEmitters[i]
			if (DimSize(thisWave, 0) == 2)
				W_xShift[offset] = thisWave[1][xCol] - thisWave[0][xCol]
				W_yShift[offset] = thisWave[1][yCol] - thisWave[0][yCol]
				offset += 1
			endif
		endfor
		
		estimatedShiftInThisRunX = mean(W_xShift)
		estimatedShiftInThisRunY = mean(W_yShift)
		
		shiftX += estimatedShiftInThisRunX
		shiftY += estimatedShiftInThisRunY
		
		nIterations += 1
	while ((estimatedShiftInThisRunX > 0.1) || (estimatedShiftInThisRunY > 0.1))
	
	dx = shiftX
	dy = shiftY
	
End

ThreadSafe Function ApplyRegistrationMap_Worker(threadIndex, nThreads, positions, M_RegistrationMap, M_CorrectedPositions, W_Progress, W_Abort)
	variable threadIndex, nThreads
	wave positions, M_RegistrationMap
	wave M_CorrectedPositions	// corrected result goes here
	wave W_Progress	// a wave containing a single point that will be incremented by all worker functions to show progress.
						// intrinsically unsafe, but since it's only for progress reporting there's no problem if the incrementing
						// goes wrong from time to time.
	wave W_Abort
	
	variable nReferencePoints = DimSize(M_RegistrationMap, 0) 
	variable nPositions = DimSize(positions, 0)
	
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	
	Make /FREE/N=(nReferencePoints)/D W_ReferenceX = M_RegistrationMap[p][0]
	Make /FREE/N=(nReferencePoints)/D W_ReferenceY = M_RegistrationMap[p][1]
	
	Make /FREE/N=6 /D W_PolyParamsX, W_PolyParamsY
	
	variable i, j, dx, dy, xLoc, yLoc, weight, summedWeights, distanceToCentroid
	variable Rn
	for (i = threadIndex; i < nPositions; i+=nThreads)
		if (W_Abort[0] != 0)
			return 0
		endif
		
		xLoc = positions[i][xCol]
		yLoc = positions[i][yCol]
		
		summedWeights = 0
		dx = 0
		dy = 0
		// go over all the reference points, and use all of them except those that are too far away (determined by Rn)
		for (j = 0; j < nReferencePoints; j+=1)
			// is this reference point too far to be of use?
			// unlike the original paper, base this on the distance to the centroid of the poly
			distanceToCentroid = sqrt((M_RegistrationMap[j][14] - xLoc)^2 + (M_RegistrationMap[j][15] - yLoc)^2)
			Rn = M_RegistrationMap[j][16]
			if (distanceToCentroid > Rn)
				continue
			endif
			
			weight = 1 - 3 * (distanceToCentroid / Rn)^2 + 2 * (distanceToCentroid / Rn)^3
			W_PolyParamsX = M_RegistrationMap[j][p + 2]
			W_PolyParamsY = M_RegistrationMap[j][p + 8]
			
			dx += weight * Poly2D(W_PolyParamsX, xLoc, yLoc)
			dy += weight * Poly2D(W_PolyParamsY, xLoc, yLoc)
			summedWeights += weight
		endfor
		
		dx /= summedWeights
		dy /= summedWeights
		
		M_CorrectedPositions[i][xCol] += dx
		M_CorrectedPositions[i][yCol] += dy
		
		W_Progress[0] += 1
	endfor
End

Function /WAVE ApplyRegistrationMap(positions, M_RegistrationMap)
	wave positions, M_RegistrationMap
	
	variable nPositions = DimSize(positions, 0)
	
	// allocate the output
	Duplicate /FREE positions, M_CorrectedPositions
	// and the wave used to report progress
	Make /FREE/D/N=1 W_Progress = 0, W_Abort = 0
	
	// set up the threads
	// assume each thread needs at least a 500 positions to be worthwhile
	variable nThreads = min(floor(nPositions / 500), ThreadProcessorCount)
	if (nThreads == 0)
		nThreads = 1
	endif
	
	variable tgID = ThreadGroupCreate(nThreads)
	
	// divide the positions to correct over the threads
	variable i
	for (i = 0; i < nThreads; i+=1)
		ThreadStart tgID, i, ApplyRegistrationMap_Worker(i, nThreads, positions, M_RegistrationMap, M_CorrectedPositions, W_Progress, W_Abort)
	endfor
	
	// now wait for the threads to finish and occasionally report progress
	variable threadStatus, nFinished
	for (;;)
		threadStatus = ThreadGroupWait(tgID, 100)
		if (threadStatus == 0)
			// all threads finished
			break
		endif
		
		// progress notification
		if (nPositions > 2500)
			nFinished = W_Progress[0]
			W_Abort[0] = LocalizerProgressWindow("Mapping", nFinished, nPositions, allowAbort=0)
		endif
	endfor
	
	threadStatus = ThreadGroupRelease(tgID)
	
	return M_CorrectedPositions
End

Function /WAVE ApplyRegistrationMap_Poly(positions, M_RegistrationMap)
	wave positions, M_RegistrationMap
	
	variable nPositions = DimSize(positions, 0)
	
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	
	// poly coefficients
	Make /FREE/N=6 /D W_xPoly, W_yPoly
	W_xPoly = M_RegistrationMap[p]
	W_yPoly = M_RegistrationMap[p + 6]
	
	// allocate the output
	Duplicate /FREE positions, M_CorrectedPositions
	
	variable i, dx, dy
	for (i = 0; i < nPositions; i+=1)
		dx = poly2D(W_xPoly, positions[i][xCol], positions[i][yCol])
		dy = poly2D(W_yPoly, positions[i][xCol], positions[i][yCol])
		
		M_CorrectedPositions[i][xCol] += dx
		M_CorrectedPositions[i][yCol] += dy
	endfor
	
	return M_CorrectedPositions
End

Function /WAVE Do3DCalibrationWork_SinglePlane(W_Positions, W_FocusDepths)
	wave /WAVE W_Positions
	wave W_FocusDepths
	
	DFREF packageFolder = root:Packages:Localizer
	
	variable nMeasurements = DimSize(W_Positions, 0)
	Assert(DimSize(W_FocusDepths, 0) == nMeasurements)
	
	// fetch the columns needed to extract information
	variable rotationCol, stdDev1Col, stdDev2Col
	wave firstPositions = W_Positions[0]
	GetColumnForRotation(firstPositions, rotationCol)
	GetColumnsForFittedWidth(firstPositions, stdDev1Col, stdDev2Col)
	if ((rotationCol == -1) || (stdDev1Col == -1) || (stdDev2Col == -1))
		Abort "All positions need to be localized using ellipsoidal fitting for astigmatism"
	endif
	
	// the user may not have entered the focus depths in order. So sort them from lowest to highest stage position
	Make /FREE /D /N=(nMeasurements) W_SortIndex
	MakeIndex W_FocusDepths, W_SortIndex
	IndexSort  W_SortIndex, W_Positions, W_FocusDepths
	
	// create waves that will visualize the measured axis ratios and theta
	Make /N=(nMeasurements) /O/D packageFolder:W_FocusPositions, packageFolder:W_StdDev1, packageFolder:W_StdDev2, packageFolder:W_RotationAngle
	wave W_FocusPositions = packageFolder:W_FocusPositions
	wave W_StdDev1 = packageFolder:W_StdDev1
	wave W_StdDev2 = packageFolder:W_StdDev2
	wave W_RotationAngle = packageFolder:W_RotationAngle
	// wave that will serve as weight and masks for the fit
	Make /N=(nMeasurements) /D /O /D W_FitWeight, W_FitMask
	
	// determine the angles of deformation of the PSF. This is determined by the rotation of the cylindrical lens
	variable i, j, nPositionsInThisPlane, avgTheta
	for (i = 0; i < nMeasurements; i+=1)
		avgTheta = 0
		wave thesePositions = W_Positions[i]
		nPositionsInThisPlane = DimSize(thesePositions, 0)
		
		for (j = 0; j < nPositionsInThisPlane; j+=1)
			// theta is always the angle to the longest axis
			avgTheta += (thesePositions[j][stdDev1Col] >= thesePositions[j][stdDev2Col]) ? thesePositions[j][rotationCol] : ((thesePositions[j][rotationCol] >= 0) ? thesePositions[j][rotationCol] - pi/2 : thesePositions[j][rotationCol] + pi/2)
		endfor
		avgTheta /= nPositionsInThisPlane
		W_RotationAngle[i] = avgTheta
	endfor
	
	// now fit these angles to a sigmoidal function to determine the principal angles (should be a difference of pi/2 relative to one another)
	variable principalAngle1
	variable principalAngle2
	variable V_FitError = 0
	CurveFit /Q/W=2/N=1 sigmoid, W_RotationAngle /X=W_FocusDepths
	if (V_FitError == 0)
		// sigmoidal fit succeeded, use its values
		wave W_Coef
		principalAngle1 = W_Coef[0]
		principalAngle2 = W_Coef[1] + W_Coef[0]
	else
		// sigmoidal fit failed, simply use the first and last value
		principalAngle1 = W_RotationAngle[0]
		principalAngle2 = W_RotationAngle[nMeasurements -1]
	endif
	
	// now loop over the positions again, and extract the shape of the fitted PSFs along the principal angles
	variable avgStdDev1, avgStdDev2, angleDiff
	for (i = 0; i < nMeasurements; i+=1)
		avgStdDev1 = 0
		avgStdDev2 = 0
		wave thesePositions = W_Positions[i]
		nPositionsInThisPlane = DimSize(thesePositions, 0)
		
		for (j = 0; j < nPositionsInThisPlane; j+=1)
			angleDiff = thesePositions[j][rotationCol] - principalAngle1
			avgStdDev1 += sqrt((thesePositions[j][stdDev1Col] * cos(angleDiff))^2 + (thesePositions[j][stdDev2Col] * sin(angleDiff))^2)
			avgStdDev2 += sqrt((thesePositions[j][stdDev1Col] * sin(angleDiff))^2 + (thesePositions[j][stdDev2Col] * cos(angleDiff))^2)
		endfor
		avgStdDev1 /= nPositionsInThisPlane
		avgStdDev2 /= nPositionsInThisPlane
		
		W_StdDev1[i] = avgStdDev1
		W_StdDev2[i] = avgStdDev2
		W_FocusPositions[i] = W_FocusDepths[i]
	endfor
	
	Make /FREE/N=13 /D W_Calibration
	W_Calibration[0] = principalAngle1
	W_Calibration[1] = WaveMin(W_FocusDepths)
	W_Calibration[2] = WaveMax(W_FocusDepths)
	
	// fit the determined standard deviations to a defocusing fit function
	Make /FREE/O/N=(5) W_DefocusFitParams
	W_DefocusFitParams[0] = 1.5
	W_DefocusFitParams[2] = 400
	W_DefocusFitParams[3] = 0
	W_DefocusFitParams[4] = 0
	// do the first fit with A and B fixed to zero, then repeat with A free
	WaveStats /Q /M=1 W_StdDev1
	W_DefocusFitParams[0] = V_min
	W_DefocusFitParams[1] = W_FocusPositions[V_minLoc]
	FuncFit /W=2/N=1/Q /H="00011" AstigDefocusFitFunction, W_DefocusFitParams, W_StdDev1 /X=W_FocusPositions
	W_DefocusFitParams[2] = Abs(W_DefocusFitParams[2])	// negative DoF can arise if A is zero
	FuncFit /W=2/N=1/Q /H="00000" AstigDefocusFitFunction, W_DefocusFitParams, W_StdDev1 /X=W_FocusPositions
	W_Calibration[3, 7] = W_DefocusFitParams[p - 3]
	
	W_DefocusFitParams[3, 4] = 0
	WaveStats /Q /M=1 W_StdDev2
	W_DefocusFitParams[0] = V_min
	W_DefocusFitParams[1] = W_FocusPositions[V_minLoc]
	FuncFit /W=2/N=1/Q /H="00011" AstigDefocusFitFunction, W_DefocusFitParams, W_StdDev2 /X=W_FocusPositions
	W_DefocusFitParams[2] = Abs(W_DefocusFitParams[2])
	FuncFit /W=2/N=1/Q /H="00000" AstigDefocusFitFunction, W_DefocusFitParams, W_StdDev2 /X=W_FocusPositions
	W_Calibration[8, 12] = W_DefocusFitParams[p - 8]
	
	return W_Calibration
End

Function AstigDefocusFitFunction(w,z) : FitFunc
	Wave w
	Variable z

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(z) = w0 * sqrt(1 + ((z - z0) / DoF)^2 + A * ((z - z0) / DoF)^3 + B * ((z - z0) / DoF)^4)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ z
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = w0
	//CurveFitDialog/ w[1] = z0
	//CurveFitDialog/ w[2] = DoF
	//CurveFitDialog/ w[3] = A
	//CurveFitDialog/ w[4] = B

	return w[0] * sqrt(1 + ((z - w[1]) / w[2])^2 + w[3] * ((z - w[1]) / w[2])^3 + w[4] * ((z - w[1]) / w[2])^4)
End

Function /WAVE Do3DCalibrationWork_MatrixBeads(W_Positions, W_FocusDepths)
	wave /WAVE W_Positions
	wave W_FocusDepths
	
	DFREF packageFolder = root:Packages:Localizer
	
	variable nMeasurements = DimSize(W_Positions, 0)
	Assert(DimSize(W_FocusDepths, 0) == nMeasurements)
	
	// the user may not have entered the focus depths in order. So sort them from lowest to highest stage position
	Make /FREE /D /N=(nMeasurements) W_SortIndex
	MakeIndex W_FocusDepths, W_SortIndex
	IndexSort  W_SortIndex, W_Positions, W_FocusDepths
	
	// combine the positions into a single wave, with each position at a different frame number
	variable nPositionsTotal = 0
	variable i
	for (i = 0; i < nMeasurements; i+=1)
		wave thesePositions = W_Positions[i]
		nPositionsTotal += DimSize(thesePositions, 0)
	endfor
	
	Make /N=(nPositionsTotal, DimSize(thesePositions, 1)) /FREE /D M_CombinedPositions
	variable offset = 0, nPositionsInThisSet
	for (i = 0; i < nMeasurements; i+=1)
		wave thesePositions = W_Positions[i]
		nPositionsInThisSet = DimSize(thesePositions, 0)
		M_CombinedPositions[offset, offset + nPositionsInThisSet - 1][] = thesePositions[p - offset][q]
		M_CombinedPositions[offset, offset + nPositionsInThisSet - 1][0] = i
		offset += nPositionsInThisSet
	endfor
	
	Note /K M_CombinedPositions, note(W_Positions[0])
	
	// fetch the columns needed to extract information
	variable rotationCol, stdDev1Col, stdDev2Col
	GetColumnForRotation(M_CombinedPositions, rotationCol)
	GetColumnsForFittedWidth(M_CombinedPositions, stdDev1Col, stdDev2Col)
	
	if ((rotationCol == -1) || (stdDev1Col == -1) || (stdDev2Col == -1))
		Abort "All positions need to be localized using ellipsoidal fitting for astigmatism"
	endif
	
	// now group positions that are close in space, and assume that they derive from the same particle
	wave /WAVE W_GroupedPositions = GroupEmitters(M_CombinedPositions, CombinedPosEstimator_Mean, 5, 0)
	
	variable nGroups = DimSize(W_GroupedPositions, 0)
	// determine the principal angles of the deformation (determined by the orientation of the cylindrical lens)
	// do this by looping over every grouped emitter, determining the angle between the horizontal axis and the axis
	// associated with the longest standard deviation, and fitting the trend as a function of focus depth to a sigmoid function
	// then average the limiting values of the sigmoid function (and apply some limits to ensure that the fit is reasonable)
	
	Duplicate /FREE /WAVE W_GroupedPositions, W_ValidGroupedPositions
	variable nPositionsInThisGroup, principalAngle1 = 0, principalAngle2 = 0, nValidGroups = 0, nValidPositions = 0
	for (i = nGroups - 1; i >= 0; i-=1)	// reversed so we can delete points from W_ValidGroupedPositions without worrying about index changes
		wave thisGroup = W_GroupedPositions[i]
		nPositionsInThisGroup = DimSize(thisGroup, 0)
		// reject groups with insufficient emitters
		if (nPositionsInThisGroup < 8)
			DeletePoints i, 1, W_ValidGroupedPositions
			continue
		endif
		Make /O /D /N=(nPositionsInThisGroup) W_Theta
		W_Theta = (thisGroup[p][stdDev1Col] >= thisGroup[p][stdDev2Col]) ? thisGroup[p][rotationCol] : ((thisGroup[p][rotationCol] >= 0) ? thisGroup[p][rotationCol] - pi/2 : thisGroup[p][rotationCol] + pi/2)
		
		variable V_FitError = 0
		CurveFit /W=2/N=1/Q sigmoid, W_Theta
		if (V_FitError != 0)	// fit failed for some reason
			DeletePoints i, 1, W_ValidGroupedPositions
			continue
		endif
		wave W_Coef
		// quality control: the midpoint must be somewhere in between the data range, and the principal angles should be pi/2 apart
		variable midPointLocation = W_Coef[2]
		variable lowerMidPointLimit = 0.2 * (nPositionsInThisGroup - 1)	// assume wave scaling start = 0, delta = 2 for W_Theta
		variable upperMidPointLimit = 0.8 * (nPositionsInThisGroup - 1)	// assume wave scaling start = 0, delta = 2 for W_Theta
		if ((midPointLocation < lowerMidPointLimit) || (midPointLocation > upperMidPointLimit))
			DeletePoints i, 1, W_ValidGroupedPositions
			continue
		endif
		variable angle1 = W_Coef[0]
		variable angle2 = W_Coef[1] + W_Coef[0]
		if ((abs(angle1 - angle2) - pi / 2) > pi / 2 * 0.3)
			DeletePoints i, 1, W_ValidGroupedPositions
			continue
		endif
		
		principalAngle1 += angle1
		principalAngle2 += angle2
		nValidGroups += 1
		nValidPositions += nPositionsInThisGroup
	endfor
	principalAngle1 /= nValidGroups
	principalAngle2 /= nValidGroups
	Assert(nValidGroups == DimSize(W_ValidGroupedPositions, 0))
	
	// W_ValidGroupedPositions now contains those grouped positions that were accepted above. Continue with those only.
	Make /O/N=(nValidPositions) /D packageFolder:W_AxisRatiosX, packageFolder:W_AxisRatiosY, packageFolder:W_RotationAngle
	wave W_AxisRatiosX = packageFolder:W_AxisRatiosX
	wave W_AxisRatiosY = packageFolder:W_AxisRatiosY
	wave W_RotationAngle = packageFolder:W_RotationAngle
	offset = 0
	for (i = 0; i < nValidGroups; i+=1)
		wave thisGroup = W_ValidGroupedPositions[i]
		nPositionsInThisGroup = DimSize(thisGroup, 0)
		
		// determine the standard deviations of each fitted ellipsoid along the principal axes
		Make /FREE /N=(nPositionsInThisGroup) /D W_ThisEmitterAxisRatiosY, W_ThisEmitterAxisRatiosX, W_principalStdDev1, W_principalStdDev2
		W_principalStdDev1 = abs(thisGroup[p][stdDev1Col] * cos(thisGroup[p][rotationCol] - principalAngle1) - thisGroup[p][stdDev2Col] * sin(thisGroup[p][rotationCol] - principalAngle1))
		W_principalStdDev2 = abs(thisGroup[p][stdDev1Col] * sin(thisGroup[p][rotationCol] - principalAngle1) + thisGroup[p][stdDev2Col] * cos(thisGroup[p][rotationCol] - principalAngle1))
		W_ThisEmitterAxisRatiosY = W_principalStdDev1[p] / W_principalStdDev2[p]
		W_ThisEmitterAxisRatiosX = W_FocusDepths[thisGroup[p][0]]
		
		// determine the center of focus for this emitter so that all emitters can be synchronized to the same focus depth
		V_FitError = 0
		CurveFit /W=2/N=1/Q sigmoid, W_ThisEmitterAxisRatiosY /X=W_ThisEmitterAxisRatiosX
		if (V_FitError != 0)
			continue
		endif
		wave W_coef
		W_coef[2] = 0	// midpoint location
		
		W_AxisRatiosY[offset, offset + nPositionsInThisGroup - 1] = W_ThisEmitterAxisRatiosY[p - offset]
		W_AxisRatiosX[offset, offset + nPositionsInThisGroup - 1] = W_ThisEmitterAxisRatiosX[p - offset] - W_coef[2]		// subtract midpoint location
		W_RotationAngle[offset, offset + nPositionsInThisGroup - 1] = thisGroup[p - offset][rotationCol]
		
		offset += nPositionsInThisGroup
	endfor
	
	Redimension /N=(offset) W_AxisRatiosY, W_AxisRatiosX
	
	// now fit the combination of all the emitters
	CurveFit /W=2/N=1/Q sigmoid, W_AxisRatiosY /X=W_AxisRatiosX
	wave W_coef
	
	W_AxisRatiosX -= W_coef[2]	// midpoint location
	W_coef[2] = 0
	
	Make /FREE/N=5 /D W_Calibration
	W_Calibration[1,4] = W_coef[p - 1]
	W_Calibration[0] = principalAngle1
	
	return W_Calibration
End

Function /WAVE Apply3DCalibrationMap(positions, W_CalibrationParams)
	wave positions, W_CalibrationParams
	
	variable nPositions = DimSize(positions, 0)
	
	// can work only with astigmatic positions
	if (NumberByKey("LOCALIZATION METHOD", note(positions)) != LOCALIZATION_ELLIPSGAUSS_ASTIG)
		Abort "ApplyCalibrationMap() requires astigmatic positions"
	endif
	
	variable rotationCol, stdDev1Col, stdDev2Col
	GetColumnForRotation(positions, rotationCol)
	GetColumnsForFittedWidth(positions, stdDev1Col, stdDev2Col)
	variable pixelSize = NumberByKey("X PIXEL SIZE", note(W_CalibrationParams))
	
	Duplicate /FREE positions, M_MappedPositions
	string waveNote = note(positions)
	waveNote = ReplaceNumberByKey("LOCALIZATION METHOD", waveNote, LOCALIZATION_ASTIG_3D)
	if (NumType(NumberByKey("X PIXEL SIZE", waveNote)) == 2)
		waveNote += "X PIXEL SIZE:" + num2str(pixelSize) + ";"
	endif
	Note /K M_MappedPositions, waveNote
	
	variable outputXCol, outputYCol, outputZCol
	variable outputXDeviation, outputYDeviation, outputZDeviation
	GetColumnsForEmitterPositions(M_MappedPositions, outputXCol, outputYCol, outputZCol)
	GetColumnsForLocalizationError(M_MappedPositions, outputXDeviation, outputYDeviation, outputZDeviation)
	
	variable principalAngle1 = W_CalibrationParams[0]
	variable minZPos = W_CalibrationParams[1]
	variable maxZPos = W_CalibrationParams[2]
	
	Make /N=(15) /FREE /D W_CalibrationAndExpParams
	W_CalibrationAndExpParams[0,12] = W_CalibrationParams[p]
	
	variable offset = 0
	variable i, angleDiff, projStdDev1, projStdDev2, zDepth
	for (i = 0; i < nPositions; i+=1)
		if ((mod(i, 1e3) == 0) && (nPositions > 25e3))
			LocalizerProgressWindow("Mapping", i, nPositions)
		endif
		angleDiff = positions[i][rotationCol] - principalAngle1
		projStdDev1 = sqrt((positions[i][stdDev1Col] * cos(angleDiff))^2 + (positions[i][stdDev2Col] * sin(angleDiff))^2)
		projStdDev2 = sqrt((positions[i][stdDev1Col] * sin(angleDiff))^2 + (positions[i][stdDev2Col] * cos(angleDiff))^2)
		
		W_CalibrationAndExpParams[13] = projStdDev1
		W_CalibrationAndExpParams[14] = projStdDev2
		
		Optimize /L=(minZPos) /H=(maxZPos) /A=0 /Q AstigStdDevCostFunction, W_CalibrationAndExpParams
		if (V_flag != 0)
			continue
		endif
		zDepth = V_minloc
		
		if ((zDepth < minZPos) || (zDepth > maxZPos))
			continue
		endif
		zDepth /= pixelSize
		
		M_MappedPositions[offset][] = positions[i][q]
		M_MappedPositions[offset][outputZCol] = zDepth
		M_MappedPositions[offset][outputZDeviation] = 0
		offset += 1
	endfor
	
	Redimension /N=(offset, -1) M_MappedPositions
	
	return M_MappedPositions
End

Function AstigStdDevCostFunction(W_CalibrationAndExpParams, proposedZDepth)
	wave W_CalibrationAndExpParams
	variable proposedZDepth
	
	// format of W_CalibrationAndExpParams is
	// principal angle
	// min z depth in calibration
	// max z depth in calibration
	// defocus params for x (5)
	// defocus params for y (5)
	// observed stdDev1
	// observed stdDev2
	variable minZ = W_CalibrationAndExpParams[1]
	variable maxZ = W_CalibrationAndExpParams[2]
	variable stdDev1 = W_CalibrationAndExpParams[13]
	variable stdDev2 = W_CalibrationAndExpParams[14]
	
	if ((proposedZDepth < minZ) || (proposedZDepth > maxZ))
		return inf
	endif
	
	Make /N=5 /FREE /D W_DefocusParams
	W_DefocusParams = W_CalibrationAndExpParams[p + 3]
	variable expectedStdDev1 = AstigDefocusFitFunction(W_DefocusParams, proposedZDepth)
	W_DefocusParams = W_CalibrationAndExpParams[p + 8]
	variable expectedStdDev2 = AstigDefocusFitFunction(W_DefocusParams, proposedZDepth)
	
	return sqrt((stdDev1 - expectedStdDev1)^2 + (stdDev2 - expectedStdDev2)^2)
End

Function /WAVE CombineEmittersIntoPosition(positions, localizationStdDev, trackLength)
	wave positions
	variable &localizationStdDev
	variable &trackLength
	// this function is called during the emitter consolidation
	// it takes a set of fitted positions from different frames,
	// thought to arise from the same emitter
	// these localizations are combined into a single one and returned
	
	// additionally the following is performed:
	// the localization coordinates are obtained as an intensity-weighted average
	// the PSF widths are combined as an intensity-weighted average
	// the integrated intensity is the sum of all the intensities
	// the localization uncertainty is set to the average divided by the square root of the number of positions
	
	// the consolidation also provides some estimate of the localization error by looking at the shifts in position that occur
	// in what is assumed to be the same emitter. This is returned by reference in localizationStdDev.
	// while this is not perfect by any means, it does provide a seemingly unbiased estimate
	// it also returns the number of positions that were combined via trackLength
	variable nPositions = DimSize(positions, 0)
	trackLength = nPositions
	variable i
	
	if (nPositions == 0)
		return positions
	endif
	
	variable nFramesPresentColumn, xCol, yCol, zCol, intensityColumn, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol
	variable xWidthCol, yWidthCol
	GetColumnForNFramesPresent(positions, nFramesPresentColumn)	// handle different types of positions
	GetColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	GetColumnsForFittedWidth(positions, xWidthCol, yWidthCol)
	GetColumnForIntegratedIntensity(positions, intensityColumn)
	GetColumnsForLocalizationError(positions, xUncertaintyCol, yUncertaintyCol, zUncertaintyCol)
	
	// make the output wave with the combined positions
	Make /N=(1, DimSize(positions, 1)) /FREE /D M_Combined
	M_Combined[0][] = positions[0][q]
	
	Make /FREE/D/N=(nPositions) W_Intensities
	if (intensityColumn >= 0)
		W_Intensities = positions[p][intensityColumn]
	else
		W_Intensities = 1
	endif
	variable summedIntensities = Sum(W_Intensities)
	
	// positions
	MatrixOP /FREE W_avgX = sum(W_Intensities * col(positions, xCol)) / summedIntensities
	M_Combined[0][xCol] = W_avgX[0]
	MatrixOP /FREE W_avgY = sum(W_Intensities * col(positions, yCol)) / summedIntensities
	M_Combined[0][yCol] = W_avgY[0]
	if (zCol >= 0)
		MatrixOP /FREE W_avgZ = sum(W_Intensities * col(positions, zCol)) / summedIntensities
		M_Combined[0][zCol] = W_avgZ[0]
	endif
	
	// fitted widths
	if (xWidthCol >= 0)
		MatrixOP /FREE W_avgXWidth = sum(W_Intensities * col(positions, xWidthCol)) / summedIntensities
		M_Combined[0][xWidthCol] = W_avgXWidth[0]
	endif
	if (yWidthCol >= 0)
		MatrixOP /FREE W_avgYWidth = sum(W_Intensities * col(positions, yWidthCol)) / summedIntensities
		M_Combined[0][yWidthCol] = W_avgYWidth[0]
	endif
	
	// summed intensity
	if (intensityColumn >= 0)
		M_Combined[0][intensityColumn] = summedIntensities
	endif
	
	// localization uncertainty
	if (xUncertaintyCol >= 0)
		MatrixOP /FREE W_avgXUncertainty = sum(W_Intensities * col(positions, xUncertaintyCol)) / summedIntensities
		M_Combined[0][xUncertaintyCol] = W_avgXUncertainty[0] / sqrt(nPositions)
	endif
	if (yUncertaintyCol >= 0)
		MatrixOP /FREE W_avgYUncertainty = sum(W_Intensities * col(positions, yUncertaintyCol)) / summedIntensities
		M_Combined[0][yUncertaintyCol] = W_avgYUncertainty[0] / sqrt(nPositions)
	endif
	if (zUncertaintyCol >= 0)
		MatrixOP /FREE W_avgZUncertainty = sum(W_Intensities * col(positions, zUncertaintyCol)) / summedIntensities
		M_Combined[0][zUncertaintyCol] = W_avgZUncertainty[0] / sqrt(nPositions)
	endif
	
	// number of frames the emitter is present
	M_Combined[0][nFramesPresentColumn] = nPositions
	
	// now try to get the localization estimate
	// we don't do intensity weighing
	MatrixOP /FREE W_summedVariance = varCols(col(positions, xCol)) + varCols(col(positions, yCol))
	localizationStdDev = sqrt(W_summedVariance[0][0])
	
	return M_Combined
End

Function /WAVE CalculateDrift_ExplicitMarkers(positions, markers, maxStepSize)
	wave positions, markers
	variable maxStepSize
	
	variable nPositions = DimSize(positions, 0)
	
	Assert(WaveExists(positions))
	Assert(WaveExists(markers))
	Assert(maxStepSize >= 0)
	
	variable nMarkers = DimSize(markers, 0)
	variable nTotalFrames = positions[nPositions - 1][0] + 1
	variable startFrame = positions[0][0]
	variable i
	variable currentFrame
	variable driftX, driftY, nearestX, nearestY, index
	variable nMarkersRecovered
	variable positionsHaveBeenConsolidated = 0	// we know this because it is checked by one of the calling routines
	
	Make /N=(nTotalFrames)/FREE/D W_RelativeDriftX = 0, W_RelativeDriftY = 0, W_nMarkersUsed = 0
	
	Make /D/O/N=(nMarkers, 2) /FREE M_previousPositions	// the positions of the markers in the previous frame, not taking drift into account
	M_previousPositions = markers[p][q]	// initialize the marker positions to the given values
	
	for (currentFrame = 0; currentFrame < nTotalFrames; currentFrame += 1)
		LocalizerProgressWindow("Estimating", currentFrame + 1, nTotalFrames)
		// extract the positions found in the current frame
		wave /WAVE M_Extracted = ExtractPositionsInFrame(positions, currentFrame, positionsHaveBeenConsolidated=positionsHaveBeenConsolidated)
		wave extractedPositions = M_Extracted[0]
		
		if (DimSize(extractedPositions, 0) == 0)	// no positions in the current frame
			continue
		endif
		
		driftX = 0
		driftY = 0
		nMarkersRecovered = 0
		for (i = 0; i < nMarkers; i += 1)
				returnPointNearestFitPositions(extractedPositions, M_previousPositions[i][0], M_previousPositions[i][1], nearestX, nearestY, index)
				if (sqrt((nearestX - M_previousPositions[i][0])^2 + (nearestY - M_previousPositions[i][1])^2) > maxStepSize)
					// it looks like this particular marker was not recovered in the fit in this frame
					// we shouldn't use it
					continue
				else	// the marker was recovered
					nMarkersRecovered += 1
					driftX += nearestX - M_previousPositions[i][0]
					driftY += nearestY - M_previousPositions[i][1]
					M_previousPositions[i][0] = nearestX
					M_previousPositions[i][1] = nearestY
				endif
		endfor
		
		if (nMarkersRecovered == 0)
			continue
		endif
		
		// take the average of the individual drifts as the absolute drift
		// and store the relative drift compared to the previous frame
		W_RelativeDriftX[currentFrame] = driftX / nMarkersRecovered
		W_RelativeDriftY[currentFrame] = driftY / nMarkersRecovered
		W_nMarkersUsed[currentFrame] += nMarkersRecovered
	endfor
	
	// convert the relative drift to an absolute drift
	Duplicate /FREE W_RelativeDriftX, W_AbsoluteDriftX
	Duplicate /FREE W_RelativeDriftY, W_AbsoluteDriftY
	W_AbsoluteDriftX[0] = W_RelativeDriftX[0]
	W_AbsoluteDriftY[0] = W_RelativeDriftY[0]
	W_AbsoluteDriftX[1,] = W_AbsoluteDriftX[p - 1] + W_RelativeDriftX[p]
	W_AbsoluteDriftY[1,] = W_AbsoluteDriftY[p - 1] + W_RelativeDriftY[p]
	
	// combine the positions waves into a single wave, and return the number of marker positions in a 2-point wave wave
	Make /N=(nTotalFrames, 3)/FREE/D M_DriftParameters
	M_DriftParameters[][0] = p	// time points of the associated drift value
	M_DriftParameters[][1] = W_AbsoluteDriftX[p]
	M_DriftParameters[][2] = W_AbsoluteDriftY[p]
	
	Make /N=2 /WAVE /FREE W_Result
	W_Result[0] = M_DriftParameters
	W_Result[1] = W_nMarkersUsed
	
	return W_Result
End

Function /WAVE DriftCorrection_automatic(positions, maxStepSize, minTrackLength, maxBlinking)
	wave positions
	variable maxStepSize, minTrackLength, maxBlinking
	
	assert(minTrackLength >= 1)
	assert(maxStepSize >= 0)
	assert(maxBlinking >= 0)
	
	// get the indices of the columns containing the x and y coordinates
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	
	variable nPositions = DimSize(positions, 0)
	variable nTotalFrames = positions[nPositions - 1][0] + 1
	
	Make /N=(nTotalFrames)/FREE/D W_RelativeDriftX = 0, W_RelativeDriftY = 0, W_nMarkersUsed = 0
	
	wave /WAVE W_GroupedEmitters = GroupEmitters(positions, CombinedPosEstimator_LastPos, maxStepSize, maxBlinking)
	
	variable nGroups = DimSize(W_GroupedEmitters, 0), nPositionsInThisGroup
	variable i, driftX, driftY, j
	for (i = 0; i < nGroups; i+=1)
		LocalizerProgressWindow("Estimating", i + 1, nGroups)
		wave thisGroup = W_GroupedEmitters[i]
		nPositionsInThisGroup = DimSize(thisGroup, 0)
		
		if (nPositionsInThisGroup < minTrackLength)
			continue
		endif
		
		for (j = 1; j < nPositionsInThisGroup; j+=1)
			W_RelativeDriftX[thisGroup[j][0]] += thisGroup[j][xCol] - thisGroup[j - 1][xCol]
			W_RelativeDriftY[thisGroup[j][0]] += thisGroup[j][yCol] - thisGroup[j - 1][yCol]
			W_nMarkersUsed[thisGroup[j][0]] += 1
		endfor
	endfor
	
	W_RelativeDriftX /= W_nMarkersUsed[p]
	W_RelativeDriftY /= W_nMarkersUsed[p]
	
	// replace NaN's
	MatrixOP /FREE W_RelativeDriftX = ReplaceNaNs(W_RelativeDriftX, 0)
	MatrixOP /FREE W_RelativeDriftY = ReplaceNaNs(W_RelativeDriftY, 0)
	
	// convert the relative drift to an absolute drift
	Duplicate /FREE W_RelativeDriftX, W_AbsoluteDriftX
	Duplicate /FREE W_RelativeDriftY, W_AbsoluteDriftY
	W_AbsoluteDriftX[0] = W_RelativeDriftX[0]
	W_AbsoluteDriftY[0] = W_RelativeDriftY[0]
	W_AbsoluteDriftX[1,] = W_AbsoluteDriftX[p - 1] + W_RelativeDriftX[p]
	W_AbsoluteDriftY[1,] = W_AbsoluteDriftY[p - 1] + W_RelativeDriftY[p]
	
	// combine the waves and return them in a 2-point wave wave also containing the number of markers used in each frame
	Make /N=(nTotalFrames, 3)/FREE/D M_DriftParameters
	M_DriftParameters[][0] = p
	M_DriftParameters[][1] = W_AbsoluteDriftX[p]
	M_DriftParameters[][2] = W_AbsoluteDriftY[p]
	
	Make /N=2 /FREE /WAVE W_Result
	W_Result[0] = M_DriftParameters
	W_Result[1] = W_nMarkersUsed
	
	return W_Result
End

Function /WAVE DriftCorrection_SubImage(positions, minPositionsPerImage, subImagePixelSize, useMultipleThreads)
	wave positions
	variable minPositionsPerImage, subImagePixelSize, useMultipleThreads
	
	variable nPositions = DimSize(positions, 0)
	variable firstFrame = positions[0][0]
	variable lastFrame = positions[nPositions - 1][0]
	variable nFrames = lastFrame - firstFrame + 1
	if (nFrames <= 1)
		Abort "Drift correction needs more than one frame containing localized emitters"
	endif
	
	// make a wave that contains the number of positions in every frame
	Make /N=(nFrames) /D /FREE W_PositionsPerFrame = 0
	SetScale /P x, firstFrame, 1, W_PositionsPerFrame
	
	variable i
	for (i = 0; i < nPositions; i += 1)
		W_PositionsPerFrame[positions[i][0] - firstFrame] += 1
	endfor
	
	// also make a wave that contains the indices at which the positions for a particular frame start
	Make /N=(nFrames) /D /FREE W_FramePositionsIndices
	SetScale /P x, firstFrame, 1, W_FramePositionsIndices
	W_FramePositionsIndices[0] = 0
	W_FramePositionsIndices[1,] = W_FramePositionsIndices[p - 1] + W_PositionsPerFrame[p- 1]
	
	// find out how many frames need to be binned together to satisfy minPhotonsPerImage
	// the calculation is simple-minded: we start from no binning and work our way upwards from there
	// also check if there is no way to meet the binning
	if (nPositions <= minPositionsPerImage)
		Abort "Not enough localized positions to meet the required positions per subimage"
	endif
	
	Make /N=0 /D /FREE W_PhotonsPerBin = 0
	variable nFramesToBin, nSubImages
	for (nFramesToBin = 1; ; nFramesToBin += 1)
		W_PhotonsPerBin = 0
		nSubImages = floor(nFrames / nFramesToBin)
		Redimension /N=(nSubImages) W_PhotonsPerBin
		for (i = 0; i < nSubImages * nFramesToBin; i += 1)
			W_PhotonsPerBin[floor(i / nFramesToBin)] += W_PositionsPerFrame[i]
		endfor
		
		if (WaveMin(W_PhotonsPerBin) > minPositionsPerImage)
			break
		endif
	endfor
	
	// Print a summary of the calculation
	Printf "Binning %d frames together, %d subimages in total\r", nFramesToBin, nSubImages
	
	// nFramesToBin is what we need to use in the calculation
	Make /N=(nSubImages + 1) /FREE /D W_FirstPositionsInBinnedFrames	// last bin contains the index to the position just after the last position to be included
	SetScale /P x, FirstFrame, nFramesToBin, W_FirstPositionsInBinnedFrames
	W_FirstPositionsInBinnedFrames[0, nSubImages - 1] = W_FramePositionsIndices[p * nFramesToBin]
	W_FirstPositionsInBinnedFrames[nSubImages] = (nFramesToBin * nSubImages == nFrames) ? nPositions :  W_FramePositionsIndices[nSubImages * nFramesToBin]
	
	// associate a time value (in acquisitions) with every of the binned frames. Assume that every binned frame is lies at a time point equal to the average of
	// all time points of the positions that it contains.
	Make /N=(nSubImages) /D /FREE W_BinnedFrameTimePoints
	MatrixOP /FREE W_PositionTimePoints = col(positions, 0)
	for (i = 0; i < nSubImages; i+=1)
		WaveStats /Q /M=1 /R=[W_FirstPositionsInBinnedFrames[i], W_FirstPositionsInBinnedFrames[i+1] - 1] W_PositionTimePoints
		W_BinnedFrameTimePoints[i] = V_avg
	endfor
	
	// do the actual calculation
	variable delta = 1
	wave M_SubImageRelativeDrift = EstimateSubimageDriftDelta(positions, subImagePixelSize, W_FirstPositionsInBinnedFrames, delta, useThreads=useMultipleThreads)
	variable nDriftEstimates = DimSize(M_SubImageRelativeDrift, 0)
	
	// convert the drift into the required output format
	Make /FREE/D/N=(nDriftEstimates + 1, 3) M_SubImageAbsoluteDrift
	M_SubImageAbsoluteDrift[0][] = 0
	M_SubImageAbsoluteDrift[1,][0] = W_BinnedFrameTimePoints(M_SubImageRelativeDrift[p-1][0])
	M_SubImageAbsoluteDrift[1,][1,2] = M_SubImageAbsoluteDrift[p - 1][q] + M_SubImageRelativeDrift[p - 1][q]
	
	return M_SubImageAbsoluteDrift
End

Function /WAVE EstimateSubimageDriftDelta(positions, subImagePixelSize, W_FirstPositionsInBinnedFrames, delta, [useThreads])
	wave positions
	variable subImagePixelSize
	wave  W_FirstPositionsInBinnedFrames
	variable delta, useThreads	// if 1 then drift between frame i and (i+1) at point i, if 2 then drift between i and (i+2) and so on. Must be positive
	
	if (ParamIsDefault(useThreads))
		useThreads = 0
	endif
	
	variable nSubImages = DimSize(W_FirstPositionsInBinnedFrames, 0) - 1
	variable nDriftEstimates = nSubImages - delta
	
	if ((delta <= 0) || (nDriftEstimates <= 0))
		Abort "Invalid delta provided to EstimateSubimageDriftDelta()"
	endif
	
	// used as output
	Make /FREE/N=(nDriftEstimates, 3) /D M_SubImageRelativeDrift = 0
	
	// used for the threading
	Make /N=(nDriftEstimates) /FREE /D W_SubImage1Indices, W_SubImage2Indices
	W_SubImage1Indices = p
	W_SubImage2Indices = p + delta
	
	// split up the calculations so multiple threads can have a go
	variable nThreads = min(ThreadProcessorCount, nDriftEstimates)
	if (useThreads == 0)
		nThreads = 1
	endif
	
	variable tgID = ThreadGroupCreate(nThreads)
	Make /FREE /D/N=1 W_Abort = 0, W_Progress = 0
	
	variable i
	for (i = 0; i < nThreads; i+=1)
		ThreadStart tgID, i,  EstimateDriftBetweenSubImageRng(i, nThreads, positions, W_subImage1Indices, W_subImage2Indices, subImagePixelSize, W_FirstPositionsInBinnedFrames, M_SubImageRelativeDrift, W_Progress, W_Abort)
	endfor
	
	variable threadsRunning, status
	for ( ; ; )
		threadsRunning = ThreadGroupWait(tgID, 100)
		if (threadsRunning == 0)
			break
		endif
		status = LocalizerProgressWindow("Estimating", W_Progress[0], nSubImages, allowAbort=0)
		if (status != 0)
			W_Abort[0] = 1
		endif
	endfor
	variable dummy = ThreadGroupRelease(tgID)
	
	return M_SubImageRelativeDrift
End

ThreadSafe Function EstimateDriftBetweenSubImageRng(threadIndex, nThreads, positions, W_subImage1Indices, W_subImage2Indices, subImagePixelSize, W_BinnedFramePositionsIndices, M_DriftOutput, W_Progress, W_Abort)
	variable threadIndex, nThreads
	wave positions
	wave W_subImage1Indices, W_subImage2Indices
	variable subImagePixelSize
	wave W_BinnedFramePositionsIndices
	wave M_DriftOutput, W_Progress, W_Abort
	
	variable nComparisons = DimSize(W_subImage1Indices, 0)
	variable i, driftX, driftY, estimationErr
	for (i = threadIndex; i < nComparisons; i += nThreads)
		if (W_Abort[0] != 0)
			return 0
		endif
		
		// find the indices of the positions needed in this iteration
		variable subImage1Index = W_subImage1Indices[i]
		variable subImage2Index = W_subImage2Indices[i]
		variable firstPositionIn1, lastPositionIn1
		variable firstPositionIn2, lastPositionIn2
		firstPositionIn1 = W_BinnedFramePositionsIndices[subImage1Index]
		lastPositionIn1 = W_BinnedFramePositionsIndices[subImage1Index + 1] - 1
		firstPositionIn2 = W_BinnedFramePositionsIndices[subImage2Index]
		lastPositionIn2 = W_BinnedFramePositionsIndices[subImage2Index + 1] - 1
		
		estimationErr = EstimateDriftBetweenSubImages(positions, subImagePixelSize, firstPositionIn1, lastPositionIn1, firstPositionIn2, lastPositionIn2, driftX, driftY)
		
		M_DriftOutput[i][0] = (subImage1Index + subImage2Index) / 2
		if (!estimationErr)
			M_DriftOutput[i][1] = driftX
			M_DriftOutput[i][2] = driftY
		else
			M_DriftOutput[i][1,2] = 0
		endif
		
		W_Progress[0] += 1
	endfor
End

ThreadSafe Function EstimateDriftBetweenSubImages(positions, subImagePixelSize, firstPositionIn1, lastPositionIn1, firstPositionIn2, lastPositionIn2, driftX, driftY)
	wave positions
	variable subImagePixelSize, firstPositionIn1, lastPositionIn1, firstPositionIn2, lastPositionIn2
	variable &driftX, &driftY
	
	// extract the positions into their own waves
	Make /FREE/N=(lastPositionIn1 - firstPositionIn1 + 1, DimSize(positions, 1)) M_Positions1
	Make /FREE/N=(lastPositionIn2 - firstPositionIn2 + 1, DimSize(positions, 1)) M_Positions2
	Note /K M_Positions1, note(positions)
	Note /K M_Positions2, note(positions)
	M_Positions1 = positions[firstPositionIn1 + p][q]
	M_Positions2 = positions[firstPositionIn2 + p][q]
	
	// and make accumulated images for both
	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(positions, xSize, ySize, pixelSize)
	wave M_AccumulatedImage1 = MakeAccumulatedImage(M_Positions1, xSize, ySize, subImagePixelSize)
	wave M_AccumulatedImage2 = MakeAccumulatedImage(M_Positions2, xSize, ySize, subImagePixelSize)
	
	// subtract the mean of both images so we can do a padded FFT
	WaveStats /Q/M=1 M_AccumulatedImage1
	variable mean1 = V_avg
	WaveStats /Q/M=1 M_AccumulatedImage2
	variable mean2 = V_avg
	M_AccumulatedImage1 -= mean1
	M_AccumulatedImage2 -= mean2
	
	variable nFFTRows = NearestConvenientFFTSize(2 * DimSize(M_AccumulatedImage1, 0))
	variable nFFTCols = NearestConvenientFFTSize(2 * DimSize(M_AccumulatedImage1, 1))
	FFT /PAD={nFFTRows, nFFTCols} /DEST=image1_FFT M_AccumulatedImage1
	FFT /PAD={nFFTRows, nFFTCols} /DEST=image2_FFT M_AccumulatedImage2
	
	// perform PCM
	variable /C smallNum = cmplx(0.01, 0)
	MatrixOP /FREE RR = (image1_FFT * conj(image2_FFT)) / mag(image1_FFT * image2_FFT + smallNum)
	KillWaves /Z M_AccumulatedImage1, M_AccumulatedImage2, image1_FFT, image2_FFT
	IFFT /DEST=RR_IFFT RR
	ImageTransform swap, RR_IFFT
	MatrixFilter /N=3 gauss, RR_IFFT
	SetScale /P x, -nFFTRows / 2 * subImagePixelSize, subImagePixelSize, RR_IFFT
	SetScale /P y, -nFFTCols / 2 * subImagePixelSize, subImagePixelSize, RR_IFFT
	
	// determine a subrange over which to fit. Assume that the drift will be described well when the fit runs over about 10 CCD pixels,
	// clamped to the size of the image
	variable fitWindowSize = round(5 / subImagePixelSize)
	variable fitStartRow = round((-5 - DimOffset(RR_IFFT, 0)) / DimDelta(RR_IFFT, 0))
	variable fitEndRow = round((5 - DimOffset(RR_IFFT, 0)) / DimDelta(RR_IFFT, 0))
	variable fitStartCol = round((-5 - DimOffset(RR_IFFT, 1)) / DimDelta(RR_IFFT, 1))
	variable fitEndCol = round((5 - DimOffset(RR_IFFT, 1)) / DimDelta(RR_IFFT, 1))
	fitStartRow = max(fitStartRow, 0)
	fitEndRow = min(fitEndRow, DimSize(RR_IFFT, 0) - 1)
	fitStartCol = max(fitStartCol, 0)
	fitEndCol = min(fitEndRow, DimSize(RR_IFFT, 1) - 1)
	
	variable V_FitError = 0
	CurveFit /Q /N=1 /W=2 Gauss2D, RR_IFFT[fitStartRow, fitEndRow][fitStartCol,fitEndCol]
	wave W_Coef
	if ((V_fitError != 0) || !Within(W_Coef[2], -5, 5) || !Within(W_Coef[4], -5, 5))
		driftX = NaN
		driftY = NaN
		return -1
	endif
	driftX = -W_Coef[2]
	driftY = -W_Coef[4]
	
	KillWaves /Z RR_IFFT
	return 0
End

ThreadSafe Function NearestConvenientFFTSize(baseSize)
	variable baseSize
	
	variable maxPowerOf2 = ceil(log(baseSize) / log(2))	// log2
	variable maxPowerOf3 = ceil(log(baseSize) / log(3))	// log3
	variable maxPowerOf5 = ceil(log(baseSize) / log(5))	// log5
	Make /FREE/D/N=((maxPowerOf2 + 1) * (maxPowerOf3 + 1) * (maxPowerOf5 + 1)) W_ConvenientSizes
	variable i, j, k, offset = 0, size
	for (i = 0; i <= maxPowerOf2; i+=1)
		for (j = 0; j <= maxPowerOf3; j+=1)
			for (k = 0; k <= maxPowerOf5; k+=1)
				size = 2^i * 3^j * 5^k
				if (mod(size, 2) != 0)
					continue
				endif
				W_ConvenientSizes[offset] = size
				offset += 1
			endfor
		endfor
	endfor
	
	Redimension /N=(offset) W_ConvenientSizes
	
	Sort W_ConvenientSizes, W_ConvenientSizes
	for (i = 0; i < DimSize(W_ConvenientSizes, 0); i+=1)
		if (W_ConvenientSizes[i] >= baseSize)
			return W_ConvenientSizes[i]
		endif
	endfor
	
	return NaN
End

Function /WAVE DriftCorrection_AvgImage(positions, nPositionsToGroup)
	wave positions
	variable nPositionsToGroup
	
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	
	variable nPositions = DimSize(positions, 0)
	
	variable nAverages = floor(nPositions / nPositionsToGroup)
	Make /FREE/N=(nAverages, 2) /D M_AveragePositions
	Make /FREE/N=(nAverages) /D W_AverageFrameLocations
	
	// use multiple threads
	variable nThreads = ThreadProcessorCount
	if (nThreads > nAverages)
		nThreads = nAverages
	endif
	
//	Make /FREE /D/N=1 W_Abort = 0, W_Progress = 0
//	DriftCorrection_AvgImage_MT(positions, nPositionsToGroup, M_AveragePositions, W_AverageFrameLocations, 0, nAverages - 1, W_Progress, W_Abort)
	
	variable tgID = ThreadGroupCreate(nThreads)
	
	Make /FREE /D/N=1 W_Abort = 0, W_Progress = 0
	variable i
	for (i = 0; i < nThreads; i+=1)
		ThreadStart tgID, i,  DriftCorrection_AvgImage_MT(i, nThreads, positions, nPositionsToGroup, M_AveragePositions, W_AverageFrameLocations, W_Progress, W_Abort)
	endfor
	
	variable threadsRunning, status
	for ( ; ; )
		threadsRunning = ThreadGroupWait(tgID, 100)
		if (threadsRunning == 0)
			break
		endif
		status = LocalizerProgressWindow("Calculating", W_Progress[0], nAverages, allowAbort=0)
		if (status != 0)
			W_Abort[0] = 1
		endif
	endfor
	variable dummy = ThreadGroupRelease(tgID)
	
//	MatrixOP /O W_AvgX = col(M_AveragePositions, 0)
//	MatrixOP /O W_AvgY = col(M_AveragePositions, 1)
	
	// modify the obtained averages so that they are expressed as deviations from the beginning
	M_AveragePositions[1,][] = M_AveragePositions[p][q] - M_AveragePositions[0][q]
	M_AveragePositions[0][] = 0
	
	// and combine them into a single wave
	Make /FREE /D /N=(nAverages, 3) M_DriftEstimates
	M_DriftEstimates[][0] = W_AverageFrameLocations[p]
	M_DriftEstimates[][1,2] = M_AveragePositions[p][q - 1]
	
	return M_DriftEstimates
	
End

ThreadSafe Function DriftCorrection_AvgImage_MT(threadIndex, nThreads, positions, nPositionsToGroup, M_AveragePositions, W_AverageFrameLocations, W_Progress, W_Abort)
	variable threadIndex, nThreads
	wave positions
	variable nPositionsToGroup
	wave M_AveragePositions, W_AverageFrameLocations		// output
	wave W_Progress, W_Abort
	
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	variable firstFrame = positions[0][0]
	variable nPositions = DimSize(positions, 0)
	variable positionsHaveBeenConsolidated = 0	// we know this because one of the earlier routines checks for this
	variable nAverages = floor(nPositions / nPositionsToGroup)
	
	variable i, j, avgX, avgY, avgFrame
	variable firstPositionInThisAverage
	
	for (i = threadIndex; i < nAverages; i += nThreads)
		if (W_Abort[0]  != 0)
			return 0
		endif
		avgX = 0
		avgY = 0
		firstPositionInThisAverage = i * nPositionsToGroup
		
		// use WaveStats to calculate the averages, using knowledge of the layout of Igor's data in memory
		WaveStats /Q /M=1 /R=[firstPositionInThisAverage, firstPositionInThisAverage + nPositionsToGroup - 1] positions
		avgFrame = V_avg
		
		WaveStats /Q /M=1 /R=[firstPositionInThisAverage + xCol * nPositions, firstPositionInThisAverage + nPositionsToGroup - 1 + xCol * nPositions] positions
		avgX = V_avg
		
		WaveStats /Q /M=1 /R=[firstPositionInThisAverage + yCol * nPositions, firstPositionInThisAverage + nPositionsToGroup - 1 + yCol * nPositions] positions
		avgY = V_avg
		
		M_AveragePositions[i][0] = avgX
		M_AveragePositions[i][1] = avgY
		W_AverageFrameLocations[i] = avgFrame
		
		W_Progress[0] += 1
	endfor
End

Function /WAVE ApplyDriftCorrection(positions, W_DriftX, W_DriftY)
	wave positions, W_DriftX, W_DriftY
	
	variable nFrames = positions[DimSize(positions, 0) - 1][0] - positions[0][0] + 1
	variable nPositions = DimSize(positions, 0)
	assert(DimSize(W_DriftX, 0) == DimSize(W_DriftY, 0))
	assert(DimSize(W_DriftX, 0) == nFrames)
	
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)

	Duplicate /FREE positions, correctedPositions
	if (nPositions < 1e5)
		correctedPositions[][xCol] -= W_DriftX(correctedPositions[p])
		correctedPositions[][yCol] -= W_DriftY(correctedPositions[p])
	else
		MultiThread correctedPositions[][xCol] -= W_DriftX(correctedPositions[p])
		MultiThread correctedPositions[][yCol] -= W_DriftY(correctedPositions[p])
	endif
	
	return correctedPositions
End

Function /WAVE SampleProbability(signal, nums)
	wave signal
	variable nums
	
	variable nRows = DimSize(signal, 0)
	variable nCols = DimSize(signal, 1)
	Make/O/N=(nums,2) points
	
	Duplicate /FREE signal, W_SortedSignal
	Redimension /N=(nRows * nCols) W_SortedSignal
	Sort W_SortedSignal, W_SortedSignal
	variable minSignal = W_SortedSignal[ceil(0.50 * nRows * nCols)]
	variable maxSignal = W_SortedSignal[floor(0.999 * nRows * nCols)]
	
	variable numFound, pp, qq, fractionalDistance
	for (numFound = 0; numFound < nums; )
		pp = round(unoise(nRows - 11)) + 5
		qq = round(unoise(nCols - 11)) + 5
		if (!Within(signal[pp][qq], minSignal, maxSignal))
			continue
		endif
		fractionalDistance = (maxSignal - signal[pp][qq]) / (maxSignal - minSignal)
		if (unoise(1) < 2 * (StatsNormalCDF(1 - fractionalDistance, 0, 0.4) - 0.5))
			if (!HaveFoundPoint(pp, qq, points, numFound))
				points[numFound][0] = pp
				points[numFound][1] = qq
				numFound += 1
			endif
		endif
	endfor
	
	return points
End

Static Function HaveFoundPoint(pp, qq, points, numFound)
	variable pp, qq
	wave points
	variable numFound
	
	variable i
	for (i = 0; i < numFound; i+=1)
		if ((points[i][0] == pp) && (points[i][1] == qq))
			return 1
		endif
	endfor
	
	return 0
End

Function /WAVE CalculateWeightsForVirtualPixel(W_AverageSignal, M_AverageCovarianceMatrix)
	Wave W_AverageSignal						// average signal for each pixel combination
	wave M_AverageCovarianceMatrix		// average covariance matrix for each pixel combination
	
	Variable nPixelCombinations = DimSize(M_AverageCovarianceMatrix,0)
	Make /N=(nPixelCombinations) /FREE/D W_WeightsForThisVirtualPixel
	if (nPixelCombinations == 1)
		W_WeightsForThisVirtualPixel[0] = 1
	ElseIf (nPixelCombinations == 2)
		MAKE/FREE/N=(nPixelCombinations - 1) MatrixB
		MAKE/FREE/N=(nPixelCombinations - 1, nPixelCombinations - 1) matrixA
		matrixB[] = M_AverageCovarianceMatrix[0][0] / W_AverageSignal[0] - M_AverageCovarianceMatrix[p+1][0] / W_AverageSignal[p+1]
		matrixA[][] = M_AverageCovarianceMatrix[p+1][q+1] / W_AverageSignal[p+1] - M_AverageCovarianceMatrix[0][q+1] / W_AverageSignal[0]
		W_WeightsForThisVirtualPixel[0] = 1
		W_WeightsForThisVirtualPixel[1] = matrixB[0] / matrixA[0][0]
	Else
		MAKE/FREE/N=(nPixelCombinations - 1) MatrixB
		MAKE/FREE/N=(nPixelCombinations - 1, nPixelCombinations - 1) matrixA
		matrixB[] = M_AverageCovarianceMatrix[0][0] / W_AverageSignal[0] - M_AverageCovarianceMatrix[p+1][0] / W_AverageSignal[p+1]
		matrixA[][] = M_AverageCovarianceMatrix[p+1][q+1] / W_AverageSignal[p+1] - M_AverageCovarianceMatrix[0][q+1] / W_AverageSignal[0]
		MatrixLinearSolve /O matrixA matrixB
		W_WeightsForThisVirtualPixel[0] = 1
		W_WeightsForThisVirtualPixel[1,*] = matrixB[p-1]
		KillWaves /Z W_IPIV	
	EndIf
	return W_WeightsForThisVirtualPixel
End

Function CovarianceMatrixForVirtualPixel(M_CombinationJackknifeValues, M_CovarianceMatrices, point)
	WAVE M_CombinationJackknifeValues
	WAVE M_CovarianceMatrices				// covariance matrix saved here
	Variable point
	
	Variable nCombinationsForThisPixel = dimsize(M_CombinationJackknifeValues, 1)
	Variable nImages = dimsize(M_CombinationJackknifeValues, 0)
	Variable i
	Variable nThreads = min(nImages, ThreadProcessorCount)
	MAKE /FREE /N=(nCombinationsForThisPixel) avs = 0
	MAKE /FREE /N=(nCombinationsForThisPixel,nCombinationsForThisPixel,nThreads) M_temps = 0
	For (i=0; i < nImages; i+=1)
		avs[] += M_CombinationJackknifeValues[i][p]
	EndFor
	avs /= nImages
	
	Make /FREE/N=(nThreads) /B W_Dummy
	MultiThread W_Dummy = CovarianceWorker(p, nThreads, M_temps, M_CombinationJackknifeValues, avs)
	MatrixOP /FREE temps = sumBeams(M_temps)

	temps /= nImages - 1
	M_CovarianceMatrices[point][][] = temps[q][r]
End

ThreadSafe Function CovarianceWorker(threadIndex, nThreads, M_temps, M_CombinationJackknifeValues, avs)
	variable threadIndex, nThreads
	wave M_temps, M_CombinationJackknifeValues, avs
	
	variable nFrames = DimSize(M_CombinationJackknifeValues, 0)
	variable i
	for (i = threadIndex; i < nFrames; i+=nThreads)
		M_temps[][][threadIndex] += (M_CombinationJackknifeValues[i][p] - avs[p])*(M_CombinationJackknifeValues[i][q] - avs[q])
	endfor
	
	return 0
End

// Calculates optimal weights
Function /WAVE PrepareWeights(comb,order,points,input)
	Variable comb	// selects number of pixel combinations to consider
	variable order		// order of the SOFI calculation
	WAVE points	// pixels to sample
	wave input		// 3D wave containing full acquired dataset (xyt)

	variable nImages = min(DimSize(input, 2), 800)
	variable totpoints = dimsize(points,0)
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct progressStruct
	progressStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
	
	// determine which pixel combinations will be used
	SOFIPixelCombinations /comb=(comb) order
	WAVE M_SOFIPixelcombinations
	Variable SOFIPixelcombinations_Size = dimsize(M_SOFIPixelcombinations,0)
	
	MAKE /FREE/N=(order,order) combocount = 0		// number of pixel combinations that give rise to each virtual pixel
	Variable i,j,point,c
	
	// determine how many pixel combinations give rise to each virtual pixel
	variable pixelOffsetX, pixelOffsetY
	For (i=0; i < SOFIPixelcombinations_Size; i+=1)
		pixelOffsetX = M_SOFIPixelcombinations[i][0]
		pixelOffsetY = M_SOFIPixelcombinations[i][1]
		combocount[pixelOffsetX][pixelOffsetY] += 1
	EndFor
	Variable maxNPixelCombinations = wavemax(combocount)
	
	// create helper waves so combinations can be easily looked up in M_SOFIPixelcombinations
	MAKE /FREE/N=(order,order,maxNPixelCombinations) comboguide = 0	// links between pixel combinations for virtual pixels and their index in M_SOFIPixelcombinations
	MAKE /FREE/N=(order,order) comborunning = 0		// temporary used to fill comboguide
	For(i=0; i < SOFIPixelcombinations_Size; i+=1)
		pixelOffsetX = M_SOFIPixelcombinations[i][0]
		pixelOffsetY = M_SOFIPixelcombinations[i][1]
		comboguide[pixelOffsetX][pixelOffsetY][comborunning[pixelOffsetX][pixelOffsetY]] = i
		comborunning[pixelOffsetX][pixelOffsetY] += 1
	EndFor
	
	MAKE/O/N=(SOFIPixelcombinations_Size) W_PixelCombinationWeights	// output from this function
	// set up intermediate waves
	Make /N=(order, order) /WAVE/FREE M_VirtualPixelJackknifes, M_VirtualPixelSignals, M_VirtualPixelCovarianceMat
	variable nCombinationsForThisPixel
	For(i=0;i<order;i+=1)
		For(j=0;j<order;j+=1)
			nCombinationsForThisPixel = combocount[i][j]
			MAKE /FREE /N=(nImages,nCombinationsForThisPixel) W_ThisPixelJackknife
			M_VirtualPixelJackknifes[i][j] = W_ThisPixelJackknife
			MAKE /FREE /N=(totpoints,nCombinationsForThisPixel) M_CombinationSignals
			M_VirtualPixelSignals[i][j] = M_CombinationSignals
			MAKE /FREE /N=(totpoints,nCombinationsForThisPixel,nCombinationsForThisPixel) M_CovarianceMatrices
			M_VirtualPixelCovarianceMat[i][j] = M_CovarianceMatrices
		EndFor
	EndFor
	
	// extract the selected pixels and their environment into a single large wave
	Make /O/N=(totPoints * 5, 5, nImages) M_CombinedSubstack
	for (point = 0; point < totpoints; point+=1)
		M_CombinedSubstack[point * 5, (point+1) * 5 - 1][] = input[p - (point * 5) +points[point][0] - 2][q+points[point][1] - 2][r]
	endfor
	
	// SOFI and SOFI jackknife
	variable doAbort
	Make /FREE/N=(totpoints, order, order, maxNPixelCombinations) /D M_Signal
	Make /FREE/N=(totpoints, order, order, maxNPixelCombinations) /WAVE M_JackknifeValues
	For(c=0; c < maxNPixelCombinations; c+=1)
		doAbort = LocalizerProgressWindow("Jackknife", c + 1, maxNPixelCombinations, allowAbort = 0)
		if (doAbort)
			return $""
		endif
		W_PixelCombinationWeights = 0
		For(i=0; i < order; i+=1)
			For(j=0; j < order; j+=1)
				If (c < combocount[i][j])	// have combinations left for this virtual pixel?
					W_PixelCombinationWeights[comboguide[i][j][c]] = 1
				Endif
			EndFor
		EndFor
		
		NewSOFI /WGHT=W_PixelCombinationWeights /ORDR=(order) /JACK /PXCR=0 /COMB=(comb) /DEST=M_SOFI_Weights /PROG=progressStruct "M_CombinedSubstack"
		WAVE M_SOFI = M_SOFI_Weights
		WAVE M_SOFI_jack = M_SOFI_Weights_jack
		
		for (point = 0; point < totpoints; point+=1)
			For(i=0; i < order; i+=1)
				For(j=0; j < order; j+=1)
					If (c < combocount[i][j])
						variable rowOffset = point * order * 5
						M_Signal[point][i][j][c] = M_SOFI[rowOffset+i][j]
						MatrixOP /FREE W_ThisJack = beam(M_SOFI_jack, rowOffset+i, j)
						M_JackknifeValues[point][i][j][c] = W_ThisJack
					Endif
				EndFor
		EndFor
		endfor
	EndFor
	KillWaves /Z M_CombinedSubstack
	
	// start looping over points
	For (point=0; point < totpoints; point+=1)
		doAbort = LocalizerProgressWindow("Pixel", point, totpoints, allowAbort = 0)
		if (doAbort)
			return $""
		endif
		
		For(c=0; c < maxNPixelCombinations; c+=1)
			For(i=0; i < order; i+=1)
				For(j=0; j < order; j+=1)
					If (c < combocount[i][j])
						WAVE CombinationSignals = M_VirtualPixelSignals[i][j]
						CombinationSignals[point][c] = M_Signal[point][i][j][c]
						WAVE M_CombinationJackknifeValues = M_VirtualPixelJackknifes[i][j]
						WAVE thisJackKnife = M_JackknifeValues[point][i][j][c]
						if (DimSize(M_CombinationJackknifeValues, 0) != DimSize(thisJackKnife, 0))
							Redimension /N=(DimSize(thisJackKnife, 0), -1) M_CombinationJackknifeValues	// needed in case the SOFI jackknife calculation didn't use all images
						endif
						M_CombinationJackknifeValues[][c] = thisJackKnife[p]
					Endif
				EndFor
			EndFor
		EndFor
		
		// make covariance matrix for samplepoint
		For (i=0; i < order; i+=1)
			For (j=0; j < order; j+=1)
				WAVE M_CombinationJackknifeValues = M_VirtualPixelJackknifes[i][j]
				WAVE M_CovarianceMatrices = M_VirtualPixelCovarianceMat[i][j]
				CovarianceMatrixForVirtualPixel(M_CombinationJackknifeValues, M_CovarianceMatrices, point)
			EndFor
		EndFor	
	EndFor
	
	For (i=0; i < order; i+=1)
		For (j=0; j < order; j+=1)
			nCombinationsForThisPixel = combocount[i][j]
			WAVE M_CombinationSignals = M_VirtualPixelSignals[i][j]
			WAVE M_CovarianceMatrices = M_VirtualPixelCovarianceMat[i][j]
			WAVE M_CombinationJackknifeValues = M_VirtualPixelJackknifes[i][j]
			
			// calculate average covariance matrix and signal wave for this virtual pixel
			MAKE/FREE/N=(nCombinationsForThisPixel) W_AverageSignal
			MAKE/FREE/N=(nCombinationsForThisPixel, nCombinationsForThisPixel) M_AverageCovarianceMatrix
			For (Point=0 ; point < totpoints; point+=1)
				W_AverageSignal += M_CombinationSignals[point][p]
				M_AverageCovarianceMatrix += M_CovarianceMatrices[point][p][q]
			EndFor
			W_AverageSignal /= totpoints
			M_AverageCovarianceMatrix /= totpoints
			
			//calculate the optimal weights for this virtual pixel
			WAVE W_WeightsForThisVirtualPixel = CalculateWeightsForVirtualPixel(W_AverageSignal, M_AverageCovarianceMatrix)
			// put weights in correct place
			variable k
			For (k=0; k < combocount[i][j]; k+=1)
				W_PixelCombinationWeights[comboguide[i][j][k]] = W_WeightsForThisVirtualPixel[k]
			EndFor
		EndFor
	EndFor
	KillWaves /Z M_SOFI, M_SOFI_jack, M_SOFIPixelCombinations
	
	return W_PixelCombinationWeights
End

Function AddEmitter(image, integral, width, xPos, yPos)
	wave image
	variable integral, width, xPos, yPos
	
	SubtractEmitter(image, -1 * integral, width, xPos, yPos)
End

Function SubtractEmitter(image, integral, width, xPos, yPos)
	wave image
	variable integral, width, xPos, yPos
	
	Assert(WaveExists(image))
	Assert((width > 0) && (xPos >= 0) && (yPos >= 0))
	
	// given image, subtract the contribution of an emitter centered at the given position with the given characteristics
	// we don't make an effort to check if the emitter is really there or not, we just assume that it is
	
	variable startP, endP
	variable startQ, endQ
	variable i, j
	variable xSize = DimSize(image, 0)
	variable ySize = DimSize(image, 1)
	variable value, distance, width_exp, amplitude
	
	startP = floor(xPos - 4 * width)
	endP = ceil(xPos + 4 * width)
	startQ = floor(yPos - 4 * width)
	endQ = ceil(yPos + 4 * width)
	
	startP = (startP >= 0) ? startP : 0
	startQ = (startQ >= 0) ? startQ : 0
	endP = (endP < xSize) ? endP : xSize - 1
	endQ = (endQ < ySize) ? endQ : ySize - 1
	
	amplitude = integral / (2 * pi * width * width)
	width_exp = 2 * width^2
	
	for (i = startP; i <= endP; i+=1)
		for (j = startQ; j <= endQ; j+=1)
			distance = (i - xPos)^2 + (j - yPos)^2
			value = amplitude * exp(-distance / width_exp)
			image[i][j] -= value
		endfor
	endfor
	
End

ThreadSafe Function GetCCDDimensionsFromPositions(positions, xSize, ySize, pixelSize)
	wave positions
	variable &xSize, &ySize, &pixelSize
	
	// given a wave containing positions, try to deduce the x and y size (in pixels) of the original data based on metadata
	// if the pixel size is known then extract that as well, otherwise return 0
	// the primary focus of this function is to isolate upstream functions from having to manipulate the wavenote directly
	
	xSize = str2num(StringByKey("X SIZE", note(positions)))
	ySize = str2num(StringByKey("Y SIZE", note(positions)))
	
	if ((NumType(xSize) == 2) || (NumType(ySize) == 2))	
		// we didn't find the xSize and/or ySize from the waveNote. 
		// Try to get the x- and y-sizes from the variables defined when loading the images
		Print "Unable to deduce the dimensions of the original images. Try re-analyzing the data"
		variable xCol, yCol, zCol
		GetColumnsForEmitterPositions(positions, xCol, yCol, zCol)
		MatrixOP /FREE W_MaxX = maxVal(col(positions, xCol))
		MatrixOP /FREE W_MaxY = maxVal(col(positions, yCol))
		xSize = ceil(W_MaxX[0] + 0.5)
		ySize = ceil(W_MaxY[0] + 0.5)
		pixelSize = 0
		return 0
	endif
	
	if (numtype(NumberByKey("X PIXEL SIZE", note(positions))) != 2)
		pixelSize = NumberByKey("X PIXEL SIZE", note(positions))
	else
		pixelSize = 0
	endif
End

ThreadSafe Function GetMinMaxCoordsFromPositions(pos, minX, maxX, minY, maxY)
	wave pos
	variable &minX, &maxX, &minY, &maxY
	
	variable nPositions = DimSize(pos, 0)
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(pos, xCol, yCol, zCol)
	
	WaveStats /Q/M=1 /R=[xCol * nPositions, (xCol + 1) * nPositions - 1] pos
	minX = V_min
	maxX = V_max
	WaveStats /Q/M=1 /R=[yCol * nPositions, (yCol + 1) * nPositions - 1] pos
	minY = V_min
	maxY = V_max
End

ThreadSafe Function GetColumnsForEmitterPositions(pos, xCol, yCol, zCol)
	wave pos
	variable &xCol, &yCol, &zCol
	
	// The different types of localization routines in the software all return their own types of parameters for every positions
	// depending on what kind of information is offered by the routine
	// This means that the column containing the x, y, and possibly z positions could be located in different columns of the positions wave
	// this function will inspect the type of positions and return the column indices containing that information by reference
	// if a particular type does not support one of these then -1 will be returned
	
	// the kind of localization used in generating the positions must be set in the wave note
	
	// special case: if LOCALIZATION METHOD is not set, then assume that the positions originate from a symmetric 2D Gaussian
	// this might be the case if the user is using positions fitted with a previous version of the package
	if (NumType(NumberByKey("LOCALIZATION METHOD", note(pos))) == 2)	// NaN
		xCol = 3
		yCol = 4
		zCol = -1
		// append the information for future use
		Note /NOCR pos, "LOCALIZATION METHOD:" + num2str(LOCALIZATION_GAUSS_FITTING) + ";"
	endif
	
	switch (NumberByKey("LOCALIZATION METHOD", note(pos)))
		case LOCALIZATION_GAUSS_FITTING:
			xCol = 3
			yCol = 4
			zCol = -1
			break
		case LOCALIZATION_GAUSS_FITTING_FIX:
			xCol = 2
			yCol = 3
			zCol = -1
			break
		case LOCALIZATION_CENTROID:
			xCol = 1
			yCol = 2
			zCol = -1
			break
		case LOCALIZATION_MULTIPLICATION:
			xCol = 2
			yCol = 3
			zCol = -1
			break
		case LOCALIZATION_ZEISSPALM:
			xCol = 2
			yCol = 3
			zCol = -1
			break
		case LOCALIZATION_ELLIPSOIDAL2DGAUSS:
			xCol = 4
			yCol = 5
			zCol = -1
			break
		case LOCALIZATION_MLEWG:
			xCol = 3
			yCol = 4
			zCol = -1
			break
		case LOCALIZATION_ELLIPSGAUSS_ASTIG:
			xCol = 4
			yCol = 5
			zCol = -1
			break
		case LOCALIZATION_ASTIG_3D:
			xCol = 4
			yCol = 5
			zCol = 6
			break
		default:
			Print "Unable to deduce the kind of localized positions (check the wavenote)"
			return -1
			break
		endswitch
End

Function ReturnPointNearestFitPositions(positionsWave, xLoc, yLoc, nearestX, nearestY, index, [zLoc, nearestZ])
	wave positionsWave
	variable xLoc, yLoc
	variable& nearestX
	variable& nearestY
	variable& index
	variable zLoc
	variable &nearestZ
	
	// given the coordinates of a point xLoc, yLoc, the function returns (by reference) the coordinates of the closest matching point in positionsWave and its index in positionsWave
	
	Assert(WaveExists(positionsWave))
	
	variable haveZ = !ParamIsDefault(zLoc)
	if (haveZ && ParamIsDefault(nearestZ))
		Abort "Need nearestZ for returnPointNearestFitPositions()"
	endif
	
	variable closest_index = -1
	variable nPos = DimSize(positionsWave, 0), xCol, yCol, zCol
	
	if (nPos == 0)
		nearestX = NaN	// return NaN if no positions are found
		nearestY = NaN
		if (haveZ)
			nearestZ = NaN
		endif
		return 1
	endif
	
	// get the indices of the columns containing the x and y coordinates
	getColumnsForEmitterPositions(positionsWave, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1))
		Abort "The positions passed to returnPointNearestFitPositions() do not appear to contain any (x,y) information"
	endif
	if (haveZ && (zCol == -1))
		Abort "Requested z info from returnPointNearestFitPositions() but positions do not have any z data"
	endif
	
	if (!haveZ)
		MatrixOP /O/FREE W_SqDistances = magSqr(col(positionsWave, xCol) - xLoc) + magSqr(col(positionsWave, yCol) - yLoc)
	else
		MatrixOP /O/FREE W_SqDistances = magSqr(col(positionsWave, xCol) - xLoc) + magSqr(col(positionsWave, yCol) - yLoc) + magSqr(col(positionsWave, zCol) - zLoc)
	endif
	WaveStats /Q/M=1 W_SqDistances
	closest_index = V_minLoc
	
	nearestX = positionsWave[closest_index][xCol]
	nearestY = positionsWave[closest_index][yCol]
	if (haveZ)
		nearestZ = positionsWave[closest_index][zCol]
	endif
	index = closest_index
	
	return 0
End

Function /WAVE GetNPointsNearestPosition(positionsWave, xLoc, yLoc, nNearestPoints)
	wave positionsWave
	variable xLoc, yLoc, nNearestPoints
	
	Assert(WaveExists(positionsWave))
	Assert((xLoc >= 0) && (yLoc >= 0) && (nNearestPoints > 0))
	
	// given the coordinates of a point x, y, the function returns  the coordinates of the nNearestPoints closest matching points in positionsWave.
	// the coordinates are returned in an nPoints x 3 wave, with the format index / nearestX / nearestY.
	// If number of point returned will be min(nNearestPoints, nPositionsInWave).
	
	variable nPos = DimSize(positionsWave, 0)
	if (nPos == 0)
		Make /N=0 /FREE /D M_NearestPositions
		return M_NearestPositions
	endif
	
	variable nPositionsToInclude = min(nNearestPoints, nPos)
	
	// get the indices of the columns containing the x and y coordinates
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positionsWave, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1))
		Abort "The positions passed to returnPointNearestFitPositions() do not appear to contain any (x,y) information"
	endif
	
	MatrixOP /O/FREE W_SqDistances = magSqr(col(positionsWave, xCol) - xLoc) + magSqr(col(positionsWave, yCol) - yLoc)
	Make /FREE/N=(nPos)/I/U W_Indices = p
	Sort W_SqDistances, W_SqDistances, W_Indices
	
	Make /FREE /N=(nPositionsToInclude, 3) /D M_NearestPositions
	M_NearestPositions[][0] = W_Indices[p]
	M_NearestPositions[][1] = positionsWave[W_Indices[p]][xCol]
	M_NearestPositions[][2] = positionsWave[W_Indices[p]][yCol]
	return M_NearestPositions
End

Function GetTrackNearestPosition(tracksWave, frameNumber, xLoc, yLoc, nearestX, nearestY, trackIndex, posIndexInTrack)
	wave /wave tracksWave
	variable frameNumber, xLoc, yLoc
	variable &nearestX, &nearestY, &trackIndex, &posIndexInTrack
	
	// Returns the indices of the track and positions within the track that is closest to (nearestX, nearestY). If frameNumber > 0 then 
	// only tracks that are active in the provided frame are considered. Otherwise all tracks are sampled.
	
	trackIndex = -1
	posIndexInTrack = -1
	
	if (frameNumber >= 0)
		wave /WAVE W_ExtractionResult = ExtractTracksInFrame(tracksWave, frameNumber)
		wave /WAVE W_ExtractedTracks = W_ExtractionResult[0]
		wave W_ExtractedTracksIndices = W_ExtractionResult[1]
	else
		wave /wave W_ExtractedTracks = tracksWave
		Make /FREE/N=(DimSize(tracksWave, 0)) W_ExtractedTracksIndices = p
	endif
	
	variable nTracksInFrame = DimSize(W_ExtractedTracks, 0)
	variable nearestTrackIndex = -1, nearestPositionIndex = -1
	variable nearestDistance = inf, thisNearestX, thisNearestY
	variable i, posIndex, distance
	for (i = 0; i < nTracksInFrame; i+=1)
		wave thesePositions = W_ExtractedTracks[i]
		returnPointNearestFitPositions(thesePositions, xLoc, yLoc, thisNearestX, thisNearestY, posIndex)
		distance = sqrt((xLoc - thisNearestX)^2 + (yLoc - thisNearestY)^2)
		if (distance < nearestDistance)
			nearestDistance = distance
			nearestX = thisNearestX
			nearestY = thisNearestY
			trackIndex = W_ExtractedTracksIndices[i]
			posIndexInTrack = posIndex
		endif
	endfor
	
	return 0
End

Function CompressReadableFilesInFolder([doFullVerify])
	variable doFullVerify
	
	if (ParamIsDefault(doFullVerify))
		doFullVerify = 1
	endif
	
	NewPath /M="Selected folder containing the files to compress" /O/Q CompressPath
	if (V_flag)
		return 0
	endif
	PathInfo CompressPath
	string folderPath = S_path
	string dataFiles
	dataFiles = IndexedFile(CompressPath, -1, ".tif")
	dataFiles += IndexedFile(CompressPath, -1, ".spe")
	dataFiles += IndexedFile(CompressPath, -1, ".his")
	dataFiles += IndexedFile(CompressPath, -1, ".sif")
	dataFiles += IndexedFile(CompressPath, -1, ".btf")
	variable nFiles = ItemsInList(dataFiles)
	KillPath CompressPath
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct progressStruct
	progressStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
	
	variable i, refNum, j
	variable outputFormat = IMAGE_OUTPUT_TYPE_COMPR_TIFF
	for (i = 0; i < nFiles; i+=1)
		LocalizerProgressWindow("file", i+1, nFiles)
		string thisPath = folderPath + StringFromList(i, dataFiles)
		string fileName = ParseFilePath(3, thisPath, ":", 0, 1)
		string tempOutputPath = folderPath + fileName + "_compressed.tif"
		string finalOutputPath = folderPath + fileName + ".tif"
		// make sure output file doesn't exist
		for (j = 0; ; j+=1)
			Open /Z/R refNum as tempOutputPath
			if (V_flag == 0)
				Close refNum
				tempOutputPath = folderPath + fileName + "_compressed" + num2istr(j) + ".tif"
			else
				break
			endif
		endfor
		ProcessCCDImages /O /M=(PROCESSCCDIMAGE_CONVERTFORMAT) /OUT=(outputFormat) /PROG=progressStruct thisPath, tempOutputPath
		
		// verify if needed
		if (doFullVerify)
			ReadCCDImages /Q/H thisPath
			variable nImages = V_numberOfImages
			ReadCCDImages /Q/H tempOutputPath
			if (nImages != V_numberOfImages)
				Abort "Information lost for " + thisPath + ", aborting"
			endif
			for (j = 0; j < nImages; j+=1)
				LocalizerProgressWindow("verify", j+1, nImages)
				ReadCCDImages /O/Q/S=(j)/C=1 /DEST=M_LocalizerOriginal thisPath
				ReadCCDImages /O/Q/S=(j)/C=1 /DEST=M_LocalizerNew tempOutputPath
				if (!EqualWaves(M_LocalizerOriginal, M_LocalizerNew, 1 + 2 + 4))
					Abort "Information lost for " + thisPath + ", aborting"
				endif
			endfor
			KillWaves /Z M_LocalizerOriginal, M_LocalizerNew
		endif
		
		// how much space did we save?
		variable oldRefNum, newRefNum, oldBytes, newBytes
		do
			Open /R/Z oldRefNum as thisPath
			if (V_flag != 0)
				break
			endif
			Open /R/Z newRefNum as tempOutputPath
			if (V_flag != 0)
				Close oldRefNum
				break
			endif
			FStatus oldRefNum
			oldBytes = V_logEOF
			FStatus newRefNum
			newBytes = V_logEOF
			Close oldRefNum
			Close newRefNum
			Printf "%s: old %g MB, new %g MB (%g%%)\r", fileName, oldBytes / 1e6, newBytes / 1e6, newBytes / oldBytes * 100
		while (0)
		
		// delete the original, and replace it with the new version
		DeleteFile thisPath
		MoveFile tempOutputPath as finalOutputPath
	endfor
End

Function GetColumnsForLocalizationError(pos, xDeviation, yDeviation, zDeviation)
	wave pos
	variable &xDeviation, &yDeviation, &zDeviation
	
	Assert(WaveExists(pos))
	
	// The different types of localization routines in the software all return their own types of parameters for every positions
	// depending on what kind of information is offered by the routine
	// This means that the column containing the deviations on the x, y, and possibly z positions could be located in different columns of the positions wave
	// this function will inspect the type of positions and return the column indices containing that information by reference
	// if a particular type does not support one of these then -1 will be returned
	
	// the kind of localization used in generating the positions must be set in the wave note
	
	// special case: if LOCALIZATION METHOD is not set, then assume that the positions originate from a symmetric 2D Gaussian
	// this might be the case if the user is using positions fitted with a previous version of the package
	if (NumType(NumberByKey("LOCALIZATION METHOD", note(pos))) == 2)	// NaN
		xDeviation = 8
		yDeviation = 9
		zDeviation = -1
		// append the information for future use
		Note /NOCR pos, "LOCALIZATION METHOD:" + num2str(LOCALIZATION_GAUSS_FITTING) + ";"
	endif
	
	switch (NumberByKey("LOCALIZATION METHOD", note(pos)))
		case LOCALIZATION_GAUSS_FITTING:
			xDeviation = 8
			yDeviation = 9
			zDeviation = -1
			break
		case LOCALIZATION_GAUSS_FITTING_FIX:
			xDeviation = 6
			yDeviation = 7
			zDeviation = -1
			break
		case LOCALIZATION_CENTROID:
			xDeviation = -1
			yDeviation = -1
			zDeviation = -1
			break
		case LOCALIZATION_MULTIPLICATION:
			xDeviation = -1
			yDeviation = -1
			zDeviation = -1
			break
		case LOCALIZATION_ZEISSPALM:
			xDeviation = 4
			yDeviation = 4
			zDeviation = -1
			break
		case LOCALIZATION_ELLIPSOIDAL2DGAUSS:
			xDeviation = 11
			yDeviation = 12
			zDeviation = -1
			break
		case LOCALIZATION_MLEWG:
			xDeviation = 6
			yDeviation = 6
			zDeviation = -1
			break
		case LOCALIZATION_ELLIPSGAUSS_ASTIG:
			xDeviation = 11
			yDeviation = 12
			zDeviation = -1
			break
		case LOCALIZATION_ASTIG_3D:
			xDeviation = 11
			yDeviation = 12
			zDeviation = 13
			break
		default:
			Abort "Unable to deduce the kind of localized positions (check the wavenote)"
			break
		endswitch
End

ThreadSafe Function GetColumnForNFramesPresent(pos, nFramesPresentCol)
	wave pos
	variable &nFramesPresentCol
	
	// The different types of localization routines in the software all return their own types of parameters for every positions
	// depending on what kind of information is offered by the routine
	// This means that the column containing the number of frames it is present could be located in a different column of the positions wave
	// this function will inspect the type of positions and return the column indices containing that information by reference
	// if a particular type does not support one of these then -1 will be returned
	
	// the kind of localization used in generating the positions must be set in the wave note
	
	// special case: if LOCALIZATION METHOD is not set, then assume that the positions originate from a symmetric 2D Gaussian
	// this might be the case if the user is using positions fitted with a previous version of the package
	if (NumType(NumberByKey("LOCALIZATION METHOD", note(pos))) == 2)	// NaN
		nFramesPresentCol = 11
		// append the information for future use
		Note /NOCR pos, "LOCALIZATION METHOD:" + num2str(LOCALIZATION_GAUSS_FITTING) + ";"
	endif
	
	switch (NumberByKey("LOCALIZATION METHOD", note(pos)))
		case LOCALIZATION_GAUSS_FITTING:
			nFramesPresentCol = 11
			break
		case LOCALIZATION_GAUSS_FITTING_FIX:
			nFramesPresentCol = 9
			break
		case LOCALIZATION_CENTROID:
			nFramesPresentCol = 3
			break
		case LOCALIZATION_MULTIPLICATION:
			nFramesPresentCol = 4
			break
		case LOCALIZATION_ZEISSPALM:
			nFramesPresentCol = 5
			break
		case LOCALIZATION_ELLIPSOIDAL2DGAUSS:
			nFramesPresentCol = 15
			break
		case LOCALIZATION_MLEWG:
			nFramesPresentCol = 7
			break
		case LOCALIZATION_ELLIPSGAUSS_ASTIG:
			nFramesPresentCol = 15
			break
		default:
			Print "Unable to deduce the kind of localized positions (check the wavenote)"
			break
		endswitch
End

Function GetColumnForIntegratedIntensity(pos, integral)
	wave pos
	variable &integral
	
	Assert(WaveExists(pos))
	
	// The different types of localization routines in the software all return their own types of parameters for every positions
	// depending on what kind of information is offered by the routine
	// This means that the column containing integrated intensity could be located in different columns of the positions wave
	// this function will inspect the type of positions and return the column indices containing that information by reference
	// if a particular type does not support one of these then -1 will be returned
	
	// the kind of localization used in generating the positions must be set in the wave note
	
	// special case: if LOCALIZATION METHOD is not set, then assume that the positions originate from a symmetric 2D Gaussian
	// this might be the case if the user is using positions fitted with a previous version of the package
	if (NumType(NumberByKey("LOCALIZATION METHOD", note(pos))) == 2)	// NaN
		integral = 1
		// append the information for future use
		Note /NOCR pos, "LOCALIZATION METHOD:" + num2str(LOCALIZATION_GAUSS_FITTING) + ";"
	endif
	
	switch (NumberByKey("LOCALIZATION METHOD", note(pos)))
		case LOCALIZATION_GAUSS_FITTING:
		case LOCALIZATION_GAUSS_FITTING_FIX:
			integral = 1
			break
		case LOCALIZATION_CENTROID:
			integral = -1
			break
		case LOCALIZATION_MULTIPLICATION:
			integral = -1
			break
		case LOCALIZATION_ZEISSPALM:
			integral = 1
			break
		case LOCALIZATION_ELLIPSOIDAL2DGAUSS:
			integral = 1
			break
		case LOCALIZATION_MLEWG:
			integral = 1
			break
		case LOCALIZATION_ELLIPSGAUSS_ASTIG:
			integral = 1
			break
		case LOCALIZATION_ASTIG_3D:
			integral = 1
			break
		default:
			Abort "Unable to deduce the kind of localized positions (check the wavenote)"
			break
		endswitch
End

Function GetColumnsForFittedWidth(pos, xWidthColumn, yWidthColumn)
	wave pos
	variable &xWidthColumn, &yWidthColumn
	
	Assert(WaveExists(pos))
	
	// The different types of localization routines in the software all return their own types of parameters for every positions
	// depending on what kind of information is offered by the routine
	// This means that the columns containing the fitted widths could be located in different columns of the positions wave
	// this function will inspect the type of positions and return the column indices containing that information by reference
	// if a particular type does not support one of these then -1 will be returned
	
	// the kind of localization used in generating the positions must be set in the wave note
	
	// special case: if LOCALIZATION METHOD is not set, then assume that the positions originate from a symmetric 2D Gaussian
	// this might be the case if the user is using positions fitted with a previous version of the package
	if (NumType(NumberByKey("LOCALIZATION METHOD", note(pos))) == 2)	// NaN
		xWidthColumn = 2
		yWidthColumn = 2
		// append the information for future use
		Note /NOCR pos, "LOCALIZATION METHOD:" + num2str(LOCALIZATION_GAUSS_FITTING) + ";"
	endif
	
	switch (NumberByKey("LOCALIZATION METHOD", note(pos)))
		case LOCALIZATION_GAUSS_FITTING:
			xWidthColumn = 2
			yWidthColumn = -1
			break
		case LOCALIZATION_GAUSS_FITTING_FIX:
			xWidthColumn = -1
			yWidthColumn = -1
			break
		case LOCALIZATION_CENTROID:
			xWidthColumn = -1
			yWidthColumn = -1
			break
		case LOCALIZATION_MULTIPLICATION:
			xWidthColumn = -1
			yWidthColumn = -1
			break
		case LOCALIZATION_ZEISSPALM:
			xWidthColumn = -1
			yWidthColumn = -1
			break
		case LOCALIZATION_ELLIPSOIDAL2DGAUSS:
			xWidthColumn = 2
			yWidthColumn = 3
			break
		case LOCALIZATION_MLEWG:
			xWidthColumn = 2
			yWidthColumn = -1
			break
		case LOCALIZATION_ELLIPSGAUSS_ASTIG:
			xWidthColumn = 2
			yWidthColumn = 3
			break
		case LOCALIZATION_ASTIG_3D:
			xWidthColumn = 2
			yWidthColumn = 3
			break
		default:
			Abort "Unable to deduce the kind of localized positions (check the wavenote)"
			break
		endswitch
End

Function GetColumnForRotation(pos, rotationCol)
	wave pos
	variable &rotationCol
	
	Assert(WaveExists(pos))
	
	switch (NumberByKey("LOCALIZATION METHOD", note(pos)))
		case LOCALIZATION_ELLIPSGAUSS_ASTIG:
			rotationCol = 6
			break
		default:
			rotationCol = -1
			break
	endswitch
End

Function /S RecursiveWaveList(startDF, nameFilter, optionsFilter)
	DFREF startDF
	string nameFilter, optionsFilter
	
	// get a list of the waves in the requested data folder
	DFREF savDF = GetDataFolderDFR()
	SetDataFolder startDF
	
	string matchingWaves = WaveList(nameFilter, ";", optionsFilter)
	variable nMatches = ItemsInList(matchingWaves)
	variable i
	string matchingWavesFullPath = ""
	
	for (i = 0; i < nMatches; i+=1)
		wave matchingWaveRef = startDF:$StringFromList(i, matchingWaves)
		matchingWavesFullPath += GetWavesDataFolder(matchingWaveRef, 2 ) + ";"
	endfor
	
	// for every data folder in this one, recurse and append the result
	string DFName
	for (i = 0; ; i+=1)
		DFName = GetIndexedObjNameDFR(startDF, 4, i)
		if (strlen(DFName) == 0)
			break
		endif
		
		DFRef subFolder = startDF:$DFName
		
		matchingWavesFullPath += RecursiveWaveList(subFolder, nameFilter, optionsFilter)
		
	endfor
	
	SetDataFolder savDF
	return matchingWavesFullPath
End

Function Assert(condition, [errMessage])
	variable condition
	string errMessage
	
	string abortMessage
	
	if (!condition)
		if (ParamIsDefault(errMessage))
			abortMessage = "An assert failed"
		else
			abortMessage = "Assert failed: " + errMessage
		endif
		Debugger
		Abort abortMessage
	endif
End

Function ZapDataInFolderTree(path)
	String path

	String savDF= GetDataFolder(1)
	SetDataFolder path

	KillWaves/A/Z
	KillVariables/A/Z
	KillStrings/A/Z

	Variable i
	Variable numDataFolders = CountObjects(":",4)
	for(i=0; i<numDataFolders; i+=1)
		String nextPath = GetIndexedObjName(":",4,i)
		ZapDataInFolderTree(nextPath)
	endfor

	SetDataFolder savDF
End

Static Function unoise(top)
	Variable top
	Variable temp = (enoise(0.5)+ 0.5) * top
	return temp
End

Function /S PrettyPrintTimeElapsed(secondsElapsed)
	variable secondsElapsed
	
	Assert(secondsElapsed >= 0)
	
	variable days = 0, hours = 0, minutes = 0, seconds = 0
	variable remainder = secondsElapsed
	string prettyString = ""
	
	if (remainder >= 24 * 3600)
		days = floor(remainder / (24 * 3600))
		if (days > 1)
			prettyString += num2str(days) + " days, "
		else
			prettyString += num2str(days) + " day, "
		endif
		remainder -= days * 24 * 3600
	endif
	
	if (remainder >= 3600)
		hours = floor(remainder / 3600)
		if (hours > 1)
			prettyString += num2str(hours) + " hours, "
		else
			prettyString += num2str(hours) + " hour, "
		endif
		remainder -= hours * 3600
	endif
	
	if (remainder >= 60)
		minutes = floor(remainder / 60)
		if (minutes > 1)
			prettyString += num2str(minutes) + " minutes, "
		else
			prettyString += num2str(minutes) + " minute, "
		endif
		remainder -= minutes * 60
	endif
	
	// always show the seconds
	if (remainder != 1)
		prettyString += num2str(remainder) + " seconds"
	else
		prettyString += num2str(remainder) + " second"
	endif
	
	return prettyString
End

ThreadSafe Function Clip(a, b, c)
	variable a, b, c
	
	return limit(a, b, c)
End

ThreadSafe Function Within(a, b, c)
	variable a, b, c
	
	return (limit(a, b, c) == a)
End

Function ExtractInstanceNumber(instanceName)
	string instanceName
	
	variable hashIndex = StrSearch(instanceName, "#", inf, 1)
	if (hashIndex == -1)
		return -1
	endif
	
	variable instanceNumber = str2num(instanceName[hashIndex + 1, strlen(instanceName) - 1])
	return instanceNumber
End


// *********************************************
// ************** Progress Viewer **************
// *********************************************

Function ProgressFunctionPrototype(num, denom)
	variable num, denom
End

Function ProgressReporterFunc(num, denom)
	variable num, denom
	return LocalizerProgressWindow("Calculating", num, denom, allowAbort = 0)
End

// a hierarchic progress window based on that posted to IgorExchange by RGerkin,
// but with entirely separate code and improved functionality

constant kLocalizerProgressBarWidth = 250
constant kLocalizerProgressBarHeight = 20
constant kLocalizerProcessLabelWidth = 100
constant kLocalizerProcessTextBoxWidth = 100
constant kLocalizerProgressMargin = 10
constant kLocalizerProgressButtonHeight = 20
constant kLocalizerProgressBarMaxDepth = 10
constant kLocalizerAbortButtonHeight = 20
constant kLocalizerAbortButtonWidth = 75

Static Function ToPixels(val)
	variable val
	// convert a dimension in points to pixels
	return val * (ScreenResolution/72)
End

Static Function ToPoints(val)
	variable val
	// convert a dimension in points to pixels
	return val * (72/ScreenResolution)
End

Function LocalizerProgressWindow(processName, num, denom, [allowAbort])
	string processName
	variable num, denom, allowAbort
	
	// if allowAbort is set then execute "Abort" when the user click the abort button
	// if not then return kUserAbort. Default is to abort.
	variable doAbort
	if (ParamIsDefault(allowAbort))
		doAbort = 1
	else
		doAbort = allowAbort
	endif
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:ProgressViewer
	
	String /G root:Packages:ProgressViewer:S_DisplayedProcesses
	SVAR displayedProcesses = root:Packages:ProgressViewer:S_DisplayedProcesses
	
	Variable /G root:Packages:ProgressViewer:V_abortRequested
	NVAR abortRequested = root:Packages:ProgressViewer:V_abortRequested
	
	wave /Z W_PreviousTicks = root:Packages:ProgressViewer:W_PreviousTicks
	if (WaveExists(W_PreviousTicks) == 0)
		Make /D/N=0 root:Packages:ProgressViewer:W_PreviousTicks
		wave W_PreviousTicks = root:Packages:ProgressViewer:W_PreviousTicks
	endif
	
	// the index in the ticks wave is the same as that of processName in displayedProcesses
	variable index
	index = WhichListItem(processName, displayedProcesses)
	if (index == -1)
		// this is a new process, there is no tick count at which it was last called
		// make an entry in the wave
		Redimension /N=(DimSize(W_PreviousTicks, 0) + 1) W_PreviousTicks
		// and store the current tick count before continuing on
		W_PreviousTicks[DimSize(W_PreviousTicks, 0) - 1] = ticks
	else
		if ((ticks - W_PreviousTicks[index] < 10) && (num != denom) && (num != denom - 1))
			return 0
		else
			W_PreviousTicks[index] = ticks
		endif
	endif
	
	string controlName
	
	// if processName contains a ':' then that will lead to errors
	if (StringMatch(processName, "*:*"))
		Abort "Error: processName cannot contain a ':')"
	endif
	if ((StringMatch(CleanupName(processName, 0), processName) != 1) || (strlen(processName) >= 30))
		Abort "Error: process names must be valid Igor names and be shorter than 30 characters)"
	endif
	
	if (abortRequested != 0)
		if (doAbort != 0)
			Abort
		else
			return kUserAbort
		endif
	endif
	
	DoWindow /F ProgressViewer
	if (V_flag != 1)
		// window does not exist
		NewPanel /N=ProgressViewer /W=(300,300, 300 + kLocalizerProcessLabelWidth + kLocalizerProgressBarWidth + kLocalizerProcessTextBoxWidth + 4 * kLocalizerProgressMargin, 330) /K=1 as "Calculating..."
		Button BTAbort, win=ProgressViewer, pos={3 * kLocalizerProgressMargin + kLocalizerProcessLabelWidth + kLocalizerProgressBarWidth, kLocalizerProgressMargin}, size={kLocalizerAbortButtonWidth, kLocalizerAbortButtonHeight}, title="Abort", proc=BTLocalizerAbortProgressViewer
		Execute /P/Q "KillDataFolder root:Packages:ProgressViewer; DoWindow /K ProgressViewer" // Automatic cleanup of data folder and progress window at the end of function execution.  
	endif
	
	// check if the current process already exists
	// otherwise set up the necessary controls
	variable controlTop, additionalWindowHeight, abortButtonTop
	if (WhichListItem(processName, displayedProcesses) == -1)
		// new process, make room
		if (ItemsInList(displayedProcesses) > kLocalizerProgressBarMaxDepth)
			Abort "Error: too many processes"
		endif
		if (ItemsInList(displayedProcesses) == 0)
			// when the window is created, Igor gives it a height of at least 30 pixels
			// so that needs to be adapted for
			additionalWindowHeight = ToPoints(kLocalizerProgressBarHeight + 3 * kLocalizerProgressMargin + kLocalizerAbortButtonHeight - 30)
			controlTop = kLocalizerProgressMargin
			abortButtonTop = controlTop + kLocalizerProgressBarHeight + kLocalizerProgressMargin
		else
			GetWindow ProgressViewer, wsize
			additionalWindowHeight = ToPoints(kLocalizerProgressBarHeight + kLocalizerProgressMargin)
			controlTop = ToPixels(V_Bottom - V_top) - (kLocalizerProgressMargin + 20)
			abortButtonTop = controlTop + kLocalizerProgressMargin + kLocalizerProgressBarHeight
		endif
		GetWindow ProgressViewer, wsize
		MoveWindow /W=ProgressViewer V_Left, V_Top, V_Right, V_Bottom + ToPoints(additionalWindowHeight)
		controlName = "TBT" +  processName
		TitleBox $controlName, win=ProgressViewer, pos={kLocalizerProgressMargin, controlTop}, fixedSize=1, size={kLocalizerProcessLabelWidth, 20}, title=processName
		controlName = "PB" + processName
		ValDisplay $controlName, win=ProgressViewer, pos={2 * kLocalizerProgressMargin + kLocalizerProcessLabelWidth, controlTop},size={kLocalizerProgressBarWidth,kLocalizerProgressBarHeight},limits={0,1,0},barmisc={0,0}, mode=3
		controlName = "TBV" + processName
		TitleBox $controlName, win=ProgressViewer, pos={3 * kLocalizerProgressMargin + kLocalizerProcessLabelWidth + kLocalizerProgressBarWidth, controlTop}, fixedSize=1, size={kLocalizerProcessTextBoxWidth, 20}, title="hello"
		Button BTAbort, win=ProgressViewer, pos={3 * kLocalizerProgressMargin + kLocalizerProcessLabelWidth + kLocalizerProgressBarWidth, abortButtonTop}
		
		displayedProcesses += processName + ";"
	endif
	
	// update the requested process
	string progressText = num2str(num) + "/" + num2str(denom)
	controlName = "PB" + processName
	ValDisplay $controlName, win=ProgressViewer, value = _NUM:(num / denom)
	controlName = "TBV" + processName
	TitleBox $controlName, win=ProgressViewer, title=progressText
	
	DoUpdate /E=1 /W=ProgressViewer
	return 0
End

Function BTLocalizerAbortProgressViewer(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
				NVAR abortRequested = root:Packages:ProgressViewer:V_abortRequested
				abortRequested = 1
			break
	endswitch

	return 0
End

// *********************************************
// **************** Color Scale Sliders ***************
// *********************************************

Function AppendColorRangeSliders()
	
	string topGraphName = WinName(0, 1, 1)
	if (strlen(topGraphName) == 0)
		Abort "No graph window to append to"
	endif
	
	string imagesInGraph = ImageNameList(topGraphName, ";")
	if (ItemsInList(imagesInGraph) == 0)
		Abort "No image is shown in the top graph"
	endif
	if (ItemsInList(imagesInGraph) > 1)
		Abort "More than one image is present; not yet supported"
	endif
	
	// if the panel is already shown, do nothing
	DoUpdate /W=topGraphName#ColorScaleSliderPanel
	if (V_flag != 0)
		return 0
	endif
	
	string imageName = StringFromList(0, imagesInGraph)
	wave imageWave = ImageNameToWaveRef(topGraphName, imageName)
	
	// initialize the panel
	NewPanel /N=ColorScaleSliderPanel/W=(0,0,360,200) /EXT=0 /HOST=$topGraphName
	PopupMenu PMSelectColorTable,win=$topGraphName#ColorScaleSliderPanel,pos={15,8},size={205,20},bodyWidth=150,title="Color Table:"
	PopupMenu PMSelectColorTable,win=$topGraphName#ColorScaleSliderPanel,mode=1,value=CtabList(),proc=PMColorTableChanged
	Slider SLUpperLimit,win=$topGraphName#ColorScaleSliderPanel,pos={18,81},size={316,13},proc=SlidersModified
	Slider SLUpperLimit,win=$topGraphName#ColorScaleSliderPanel,limits={0,2,1},value= 0,side= 0,vert= 0
	CheckBox CBReverseColors,win=$topGraphName#ColorScaleSliderPanel,pos={243,10},size={98,14},title="Reverse colortable"
	CheckBox CBReverseColors,win=$topGraphName#ColorScaleSliderPanel,value= 0,proc=CBReverseColorTable
	GroupBox GBUpperLimit,win=$topGraphName#ColorScaleSliderPanel,pos={12,34},size={335,73},title="Upper limit"
	SetVariable SVUpperLimitSliderMin,win=$topGraphName#ColorScaleSliderPanel,pos={23,58},size={75,15},proc=SVSetVariablesModified
	SetVariable SVUpperLimitValue,win=$topGraphName#ColorScaleSliderPanel,pos={140,57},size={75,15},proc=SVSetVariablesModified
	SetVariable SVUpperLimitSliderMax,win=$topGraphName#ColorScaleSliderPanel,pos={263,57},size={75,15},proc=SVSetVariablesModified
	Slider SLLowerLimit,win=$topGraphName#ColorScaleSliderPanel,pos={19,160},size={316,13},proc=SlidersModified
	Slider SLLowerLimit,win=$topGraphName#ColorScaleSliderPanel,limits={0,2,1},value= 0,side= 0,vert= 0
	GroupBox GBLowerLimit,win=$topGraphName#ColorScaleSliderPanel,pos={13,113},size={335,73},title="Lower limit"
	SetVariable SVLowerLimitSliderMin,win=$topGraphName#ColorScaleSliderPanel,pos={24,137},size={75,15},proc=SVSetVariablesModified
	SetVariable SVLowerLimitValue,win=$topGraphName#ColorScaleSliderPanel,pos={141,136},size={75,15},proc=SVSetVariablesModified
	SetVariable SVLowerLimitSliderMax,win=$topGraphName#ColorScaleSliderPanel,pos={264,136},size={75,15},proc=SVSetVariablesModified
	
	// update the limits of the sliders to the min and max of the image
	variable minLimit = WaveMin(imageWave)
	variable maxLimit = WaveMax(imageWave)
	
	Slider SLUpperLimit,win=$topGraphName#ColorScaleSliderPanel,limits={minLimit, maxLimit, 0}
	Slider SLLowerLimit,win=$topGraphName#ColorScaleSliderPanel,limits={minLimit, maxLimit, 0}
	SetVariable SVUpperLimitSliderMax,win=$topGraphName#ColorScaleSliderPanel,value=_NUM:maxLimit
	SetVariable SVUpperLimitSliderMin,win=$topGraphName#ColorScaleSliderPanel,value=_NUM:minLimit
	SetVariable SVLowerLimitSliderMax,win=$topGraphName#ColorScaleSliderPanel,value=_NUM:maxLimit
	SetVariable SVLowerLimitSliderMin,win=$topGraphName#ColorScaleSliderPanel,value=_NUM:minLimit
		
	
	// and update the controls to the actual settings
	UpdateControlsToActualSettings(topGraphName)
End

Static Function UpdateControlsToActualSettings(windowName)
	string windowName
	
	string baseWindowName = StringFromList(0, windowName, "#")
	
	string imagesInGraph = ImageNameList(baseWindowName, ";")
	if (ItemsInList(imagesInGraph) == 0)
		Abort "No image is shown in the top graph"
	endif
	if (ItemsInList(imagesInGraph) > 1)
		Abort "More than one image is present; not yet supported"
	endif
	
	string imageName = StringFromList(0, imagesInGraph)
	wave imageWave = ImageNameToWaveRef(baseWindowName, imageName)
	
	// get the color scale that is currently selected
	string recreation = StringByKey("RECREATION", ImageInfo(baseWindowName, imageName, 0))
	string ctabStr = StringByKey("ctab", recreation, "=")
	
	string lowerLimitStr, upperLimitStr, ctabName, isReverseStr
	variable lowerLimit, upperLimit, isReverse
	
	SplitString  /E="\\{([^,]*),([^,]*),([^,]*),([0-9]+)\\}" ctabStr, lowerLimitStr, upperLimitStr, ctabName, isReverseStr
	
	if (GrepString(lowerLimitStr, "\\*"))
		lowerLimit = WaveMin(imageWave)
	else
		lowerLimit = str2num(lowerLimitStr)
	endif
	
	if (GrepString(upperLimitStr, "\\*"))
		upperLimit = WaveMax(imageWave)
	else
		upperLimit = str2num(upperLimitStr)
	endif
	
	isReverse = str2num(isReverseStr)
	
	// and update the controls to the actual settings
	PopupMenu PMSelectColorTable,win=$baseWindowName#ColorScaleSliderPanel,popmatch=ctabName
	CheckBox CBReverseColors,win=$baseWindowName#ColorScaleSliderPanel,value=isReverse
	Slider SLUpperLimit,win=$baseWindowName#ColorScaleSliderPanel,value=upperLimit
	Slider SLLowerLimit,win=$baseWindowName#ColorScaleSliderPanel,value=lowerLimit
	
	SetVariable SVUpperLimitValue,win=$baseWindowName#ColorScaleSliderPanel,value=_NUM:upperLimit
	SetVariable SVLowerLimitValue,win=$baseWindowName#ColorScaleSliderPanel,value=_NUM:lowerLimit

End

Function SlidersModified(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				string baseWindowName = StringFromList(0, sa.win, "#")
				
				// update both SetVariables
				ControlInfo /W=$sa.win SLUpperLimit
				SetVariable SVUpperLimitValue,win=$sa.win,value=_NUM:V_Value
				
				ControlInfo /W=$sa.win SLLowerLimit
				SetVariable SVLowerLimitValue,win=$sa.win,value=_NUM:V_Value
				
				UpdateColorScaling(sa.win)
			endif
			break
	endswitch

	return 0
End

Function SVSetVariablesModified(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			
			ControlInfo /W=$sva.win SVUpperLimitSliderMax
			variable upperLimitMax = V_Value
			ControlInfo /W=$sva.win SVUpperLimitSliderMin
			variable upperLimitMin = V_Value
			ControlInfo /W=$sva.win SVUpperLimitValue
			variable upperLimit = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMax
			variable lowerLimitMax = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMin
			variable lowerLimitMin = V_Value
			ControlInfo /W=$sva.win SVLowerLimitValue
			variable lowerLimit = V_Value
			
			Slider SLUpperLimit,win=$sva.win,value=upperLimit,limits={upperLimitMin, upperLimitMax, 0}
			Slider SLLowerLimit,win=$sva.win,value=lowerLimit,limits={lowerLimitMin, lowerLimitMax, 0}
			
			UpdateColorScaling(sva.win)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function PMColorTableChanged(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			UpdateColorScaling(pa.win)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CBReverseColorTable(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			
			UpdateColorScaling(cba.win)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Static Function UpdateColorScaling(windowName)
	string windowName
	
	string baseWindowName = StringFromList(0, windowName, "#")
	
	string imagesInGraph = ImageNameList(baseWindowName, ";")
	if (ItemsInList(imagesInGraph) == 0)
		Abort "No image is shown in the top graph"
	endif
	if (ItemsInList(imagesInGraph) > 1)
		Abort "More than one image is present; not yet supported"
	endif
	
	string imageName = StringFromList(0, imagesInGraph)
	wave imageWave = ImageNameToWaveRef(baseWindowName, imageName)
	
	// get the color scale to use
	string ctabName
	ControlInfo /W=$baseWindowName#ColorScaleSliderPanel PMSelectColorTable
	ctabName = S_Value
	
	variable minVal, maxVal
	ControlInfo /W=$baseWindowName#ColorScaleSliderPanel SLUpperLimit
	maxVal = V_Value
	
	ControlInfo /W=$baseWindowName#ColorScaleSliderPanel SLLowerLimit
	minVal = V_Value
	
	variable isReverse
	ControlInfo /W=$baseWindowName#ColorScaleSliderPanel CBReverseColors
	isReverse = V_Value
	
	ModifyImage /W=$baseWindowName $imageName, ctab={minVal, maxVal, $ctabName, isReverse}
End

// *********************************************
// ****************** RGB Merging *****************
// *********************************************

Function RunRGBMerge()
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:ColorScaleSliders
	
	DFREF savDF = GetDataFolderDFR()
	
	// the problem with image waves is that we need to consider both mxn and mxnx1 waves
	// do some hackery with that
	string candidateWaves = RecursiveWaveList(savDF, "*", "CMPLX:0,TEXT:0")
	string eligibleWaves = ""
	variable i
	for (i = 0; i < ItemsInList(candidateWaves); i+=1)
		wave thisWave = $StringFromList(i, candidateWaves)
		if ((DimSize(thisWave, 0) != 0) && (DimSize(thisWave, 1) != 0) && ((DimSize(thisWave, 2) == 0) || (DimSize(thisWave, 2) == 1)))
			eligibleWaves += StringFromList(i, candidateWaves) + ";"
		endif
	endfor
	
	string promptString = "_none_;" + eligibleWaves
	
	string redWaveName, greenWaveName, blueWaveName
	Prompt redWaveName, "Red channel:", popup, promptString
	Prompt greenWaveName, "Green channel:", popup, promptString
	Prompt blueWaveName, "Blue channel:", popup, promptString
	
	DoPrompt "Select the channels", redWaveName, greenWaveName, blueWaveName
	if (V_flag == 1)
		return 0
	endif
	
	if (StringMatch(redWaveName, "_none_") && StringMatch(greenWaveName, "_none_") && StringMatch(redWaveName, "_none_"))
		return 0
	endif
	
	variable useRed = (StringMatch(redWaveName, "_none_")) ? 0 : 1
	variable useGreen = (StringMatch(greenWaveName, "_none_")) ? 0 : 1
	variable useBlue = (StringMatch(blueWaveName, "_none_")) ? 0 : 1
	
	// get a unique name for the graph window
	string displayWindowName = UniqueName("RGBMerge", 6, 0)
	
	// now that we have a unique window name, set up the data folder
	// associated with this window
	NewDataFolder /O root:Packages:ColorScaleSliders:$displayWindowName
	DFREF windowDataFolder = root:Packages:ColorScaleSliders:$displayWindowName
	
	// the datafolder may already exist if another window was created and closed previously
	// so be sure to empty it first
	DFREF savDF = GetDataFolderDFR()
	SetDataFolder windowDataFolder
	KillWaves /A/Z
	SetDataFolder savDF
	
	// global strings that will contain the names of the waves to use
	string /G windowDataFolder:S_RedWaveName
	string /G windowDataFolder:S_GreenWaveName
	string /G windowDataFolder:S_BlueWaveName
	SVAR globalRedWaveName = windowDataFolder:S_RedWaveName
	SVAR globalGreenWaveName = windowDataFolder:S_GreenWaveName
	SVAR globalBlueWaveName = windowDataFolder:S_BlueWaveName
	
	globalRedWaveName = redWaveName
	globalGreenWaveName = greenWaveName
	globalBlueWaveName = blueWaveName
	wave /Z redWave = $redWaveName
	wave /Z greenWave = $greenWaveName
	wave /Z blueWave = $blueWaveName
	
	// Make a merged RGB wave
	if (useRed)
		Duplicate /O redWave, windowDataFolder:M_Merged_RGB
	elseif (useGreen)
		Duplicate /O greenWave, windowDataFolder:M_Merged_RGB
	elseif (useBlue)
		Duplicate /O blueWave, windowDataFolder:M_Merged_RGB
	endif
	Redimension /W/U/N=(-1,-1,3) windowDataFolder:M_Merged_RGB
	
	// display the waves in a single image
	Display /N=$displayWindowName
	AppendImage /W=$displayWindowName windowDataFolder:M_Merged_RGB
	ModifyGraph /W=$displayWindowName margin(left)=28,margin(bottom)=28,margin(top)=14,margin(right)=14,width={Aspect,1}, mirror=2, nticks=10, minor=1, fSize=8, standoff=0, tkLblRot(left)=90, btLen=3, tlOffset=-2
	ModifyGraph /W=$displayWindowName width={Aspect, DimSize(windowDataFolder:M_Merged_RGB, 0) / DimSize(windowDataFolder:M_Merged_RGB, 1)}
	
	// add the exterior subwindow with the controls
	string redTitle = "Red Channel: " + NameOfWave(redWave)
	string greenTitle = "Green Channel: " + NameOfWave(greenWave)
	string blueTitle = "Blue Channel: " + NameOfWave(blueWave)
	
	NewPanel /N=ColorScaleSliderPanel/W=(0,0,360,350) /EXT=0 /HOST=$displayWindowName
	GroupBox GBRed,pos={12,4},size={335,109},title=redTitle
	Slider SLUpperLimitRed,pos={18,50},size={316,13},disable=(!useRed * 2)
	Slider SLUpperLimitRed,limits={0,2,1},value= 0,side= 0,vert= 0,proc=SLRGBMergeProc
	SetVariable SVUpperLimitSliderMinRed,pos={23,26},size={75,15},disable=(!useRed * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVUpperLimitValueRed,pos={140,26},size={75,15},disable=(!useRed * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVUpperLimitSliderMaxRed,pos={263,26},size={75,15},disable=(!useRed * 2),proc=SVRGBMergeSetVariablesModified
	Slider SLLowerLimitRed,pos={19,94},size={316,13},disable=(!useRed * 2)
	Slider SLLowerLimitRed,limits={0,2,1},value= 0,side= 0,vert= 0,proc=SLRGBMergeProc
	SetVariable SVLowerLimitSliderMinRed,pos={24,70},size={75,15},disable=(!useRed * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVLowerLimitValueRed,pos={141,70},size={75,15},disable=(!useRed * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVLowerLimitSliderMaxRed,pos={264,70},size={75,15},disable=(!useRed * 2),proc=SVRGBMergeSetVariablesModified
	
	GroupBox GBGreen,pos={12,120},size={335,109},title=greenTitle
	Slider SLUpperLimitGreen,pos={18,166},size={316,13},disable=(!useGreen * 2)
	Slider SLUpperLimitGreen,limits={0,2,1},value= 0,side= 0,vert= 0,proc=SLRGBMergeProc
	SetVariable SVUpperLimitSliderMinGreen,pos={23,142},size={75,15},disable=(!useGreen * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVUpperLimitValueGreen,pos={140,142},size={75,15},disable=(!useGreen * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVUpperLimitSliderMaxGreen,pos={263,142},size={75,15},disable=(!useGreen * 2),proc=SVRGBMergeSetVariablesModified
	Slider SLLowerLimitGreen,pos={19,210},size={316,206},disable=(!useGreen * 2)
	Slider SLLowerLimitGreen,limits={0,2,1},value= 0,side= 0,vert= 0,proc=SLRGBMergeProc
	SetVariable SVLowerLimitSliderMinGreen,pos={24,186},size={75,15},disable=(!useGreen * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVLowerLimitValueGreen,pos={141,186},size={75,15},disable=(!useGreen * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVLowerLimitSliderMaxGreen,pos={264,186},size={75,15},disable=(!useGreen * 2),proc=SVRGBMergeSetVariablesModified
	
	GroupBox GBBlue,pos={12,236},size={335,109},title=blueTitle
	Slider SLUpperLimitBlue,pos={18,282},size={316,13},disable=(!useBlue * 2)
	Slider SLUpperLimitBlue,limits={0,2,1},value= 0,side= 0,vert= 0,proc=SLRGBMergeProc
	SetVariable SVUpperLimitSliderMinBlue,pos={23,258},size={75,15},disable=(!useBlue * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVUpperLimitValueBlue,pos={140,258},size={75,15},disable=(!useBlue * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVUpperLimitSliderMaxBlue,pos={263,258},size={75,15},disable=(!useBlue * 2),proc=SVRGBMergeSetVariablesModified
	Slider SLLowerLimitBlue,pos={19,326},size={316,206},disable=(!useBlue * 2)
	Slider SLLowerLimitBlue,limits={0,2,1},value= 0,side= 0,vert= 0,proc=SLRGBMergeProc
	SetVariable SVLowerLimitSliderMinBlue,pos={24,302},size={75,15},disable=(!useBlue * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVLowerLimitValueBlue,pos={141,302},size={75,15},disable=(!useBlue * 2),proc=SVRGBMergeSetVariablesModified
	SetVariable SVLowerLimitSliderMaxBlue,pos={264,302},size={75,15},disable=(!useBlue * 2),proc=SVRGBMergeSetVariablesModified
	
	// provide some default values
	variable redMin, redMax
	variable greenMin, greenMax
	variable blueMin, blueMax
	if (useRed)
		redMin = WaveMin(redWave)
		redMax = WaveMax(redWave)
		Slider SLUpperLimitRed,limits={redMin,redMax,0},value= redMax
		Slider SLLowerLimitRed,limits={redMin,redMax,0},value= redMin
		SetVariable SVUpperLimitSliderMinRed,value=_NUM:redMin
		SetVariable SVUpperLimitValueRed,value=_NUM:redMax
		SetVariable SVUpperLimitSliderMaxRed,value=_NUM:redMax
		SetVariable SVLowerLimitSliderMinRed,value=_NUM:redMin
		SetVariable SVLowerLimitValueRed,value=_NUM:redMin
		SetVariable SVLowerLimitSliderMaxRed,value=_NUM:redMax
	endif
	if (useGreen)
		greenMin = WaveMin(greenWave)
		greenMax = WaveMax(greenWave)
		Slider SLUpperLimitGreen,limits={greenMin,greenMax,0},value= greenMax
		Slider SLLowerLimitGreen,limits={greenMin,greenMax,0},value= greenMin
		SetVariable SVUpperLimitSliderMinGreen,value=_NUM:greenMin
		SetVariable SVUpperLimitValueGreen,value=_NUM:greenMax
		SetVariable SVUpperLimitSliderMaxGreen,value=_NUM:greenMax
		SetVariable SVLowerLimitSliderMinGreen,value=_NUM:greenMin
		SetVariable SVLowerLimitValueGreen,value=_NUM:greenMin
		SetVariable SVLowerLimitSliderMaxGreen,value=_NUM:greenMax
	endif
	if (useBlue)
		blueMin = WaveMin(blueWave)
		blueMax = WaveMax(blueWave)
		Slider SLUpperLimitBlue,limits={blueMin,blueMax,0},value= blueMax
		Slider SLLowerLimitBlue,limits={blueMin,blueMax,0},value= blueMin
		SetVariable SVUpperLimitSliderMinBlue,value=_NUM:blueMin
		SetVariable SVUpperLimitValueBlue,value=_NUM:blueMax
		SetVariable SVUpperLimitSliderMaxBlue,value=_NUM:blueMax
		SetVariable SVLowerLimitSliderMinBlue,value=_NUM:blueMin
		SetVariable SVLowerLimitValueBlue,value=_NUM:blueMin
		SetVariable SVLowerLimitSliderMaxBlue,value=_NUM:blueMax
	endif
	
	RGBMergeUpdate(displayWindowName)
End

Static Function RGBMergeUpdate(windowName)
	string windowName
	
	DFREF windowDataFolder = root:Packages:ColorScaleSliders:$windowName
	
	// update all of the RGB waves, and calculate a new merged image
	// first get the actual settings
	wave mergedImage = windowDataFolder:M_Merged_RGB
	
	variable minValue, maxValue
	string colors = "Red;Green;Blue;"
	variable i
	for (i = 0; i < ItemsInList(colors); i+=1)
		SVAR globalWaveName = windowDataFolder:$("S_" + StringFromList(i, colors) + "WaveName")
		if ((strlen(globalWaveName) != 0) && (StringMatch(globalWaveName, "_none_") != 1))
			ControlInfo /W=$windowName#ColorScaleSliderPanel $("SLUpperLimit" + StringFromList(i, colors))
			maxValue = V_Value
			ControlInfo /W=$windowName#ColorScaleSliderPanel $("SLLowerLimit" + StringFromList(i, colors))
			minValue = V_Value
			
			wave dataWave = $globalWaveName
			wave displayWave = windowDataFolder:M_Merged_RGB
			MatrixOP /O/FREE M_ThisImagePlane = (clip(dataWave, minValue, maxValue) - minValue) / (maxValue - minValue) * 65535
			ImageTransform /P=(i) /D=M_ThisImagePlane setPlane, mergedImage
			
		endif
	endfor
End

Function SLRGBMergeProc(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				
				// update the setvariable to the correct value
				string controlName = sa.ctrlName
				string baseWindowName = StringFromList(0, sa.win, "#")
				string setVariableName
				if (StringMatch(controlName, "*UpperLimit*"))
					setVariableName = ReplaceString("SLUpperLimit", controlName, "SVUpperLimitValue")
				else
					setVariableName = ReplaceString("SLLowerLimit", controlName, "SVLowerLimitValue")
				endif
				
				SetVariable $setVariableName,win=$baseWindowName#ColorScaleSliderPanel,value=_NUM:curval
				
				RGBMergeUpdate(baseWindowName)
			endif
			break
	endswitch

	return 0
End

Function SVRGBMergeSetVariablesModified(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			string baseWindowName = StringFromList(0, sva.win, "#")
			string controlName = sva.ctrlName
			// adjust the sliders to the values of the SetVariables
			
			variable upperLimitMax, upperLimitMin, upperLimit, lowerLimitMax, lowerLimitMin, lowerLimit
			
			//Red
			ControlInfo /W=$sva.win SVUpperLimitSliderMaxRed
			upperLimitMax = V_Value
			ControlInfo /W=$sva.win SVUpperLimitSliderMinRed
			upperLimitMin = V_Value
			ControlInfo /W=$sva.win SVUpperLimitValueRed
			upperLimit = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMaxRed
			lowerLimitMax = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMinRed
			lowerLimitMin = V_Value
			ControlInfo /W=$sva.win SVLowerLimitValueRed
			lowerLimit = V_Value
			
			Slider SLUpperLimitRed,win=$sva.win,value=upperLimit,limits={upperLimitMin, upperLimitMax, 0}
			Slider SLLowerLimitRed,win=$sva.win,value=lowerLimit,limits={lowerLimitMin, lowerLimitMax, 0}
			
			//Green
			ControlInfo /W=$sva.win SVUpperLimitSliderMaxGreen
			upperLimitMax = V_Value
			ControlInfo /W=$sva.win SVUpperLimitSliderMinGreen
			upperLimitMin = V_Value
			ControlInfo /W=$sva.win SVUpperLimitValueGreen
			upperLimit = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMaxGreen
			lowerLimitMax = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMinGreen
			lowerLimitMin = V_Value
			ControlInfo /W=$sva.win SVLowerLimitValueGreen
			lowerLimit = V_Value
			
			Slider SLUpperLimitGreen,win=$sva.win,value=upperLimit,limits={upperLimitMin, upperLimitMax, 0}
			Slider SLLowerLimitGreen,win=$sva.win,value=lowerLimit,limits={lowerLimitMin, lowerLimitMax, 0}
			
			//Blue
			ControlInfo /W=$sva.win SVUpperLimitSliderMaxBlue
			upperLimitMax = V_Value
			ControlInfo /W=$sva.win SVUpperLimitSliderMinBlue
			upperLimitMin = V_Value
			ControlInfo /W=$sva.win SVUpperLimitValueBlue
			upperLimit = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMaxBlue
			lowerLimitMax = V_Value
			ControlInfo /W=$sva.win SVLowerLimitSliderMinBlue
			lowerLimitMin = V_Value
			ControlInfo /W=$sva.win SVLowerLimitValueBlue
			lowerLimit = V_Value
			
			Slider SLUpperLimitBlue,win=$sva.win,value=upperLimit,limits={upperLimitMin, upperLimitMax, 0}
			Slider SLLowerLimitBlue,win=$sva.win,value=lowerLimit,limits={lowerLimitMin, lowerLimitMax, 0}
			
			RGBMergeUpdate(baseWindowName)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
