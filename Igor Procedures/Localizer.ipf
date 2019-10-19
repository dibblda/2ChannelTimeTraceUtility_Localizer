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
#pragma IndependentModule = Localizer

#include <Resize Controls>
#include <WaveSelectorWidget>
#include <SaveRestoreWindowCoords>
#include <ImageSlider>
#include "LocalizerUtilities"

Menu "Localizer"
	SubMenu "Read CCD Data"
		"Read CCD data from disk...", /Q, SetupImageViewer_loadFromDisk()
		"Read multi-file TIFF data...", /Q, SetupImageViewer_multiFileTIFF()
		"Read CCD data from Igor wave...", /Q, SetupImageViewer_IgorWaves()
		SubMenu "Special"
			"Read entire CCD file into memory...", /Q, ReadCCDFramesIntoMemory()
			"Read average image from...", /Q, OpenAverageImages()
		End
		"-"
		"Load localized positions from file...", /Q, LoadPositionsFromTextFile()
	End
	SubMenu "Manipulate localized positions"
		"Merge positions...", /Q, MergePositions()
		"Consolidate identical emitters...", /Q, ConsolidateSameEmitters_menu()
		"Remove outlier positions...", /Q, RemoveOutlierPositions_Menu()
		"Drift correction...", /Q, DoDriftCorrection()
		SubMenu "Two-Color Imaging"
			"Create registration map...", /Q, CreateRegistrationMap_Menu()
			"Apply registration map...", /Q, ApplyRegistrationMap_Menu()
		End
		SubMenu "3D Imaging"
			"Create 3D calibration...", /Q, Create3DCalibration_Menu()
			"Apply 3D calibration...", /Q, Apply3DCalibrationMap_Menu()
		End
		"-"
		"Save positions to text file...", /Q, SavePositionsToTextFile()
	End
	SubMenu "Analyze Localized Positions"
		"L-Clustering analysis...", /Q, DoClusteringAnalysis(1, 0)
		"Pairwise Correlation...", /Q, DoClusteringAnalysis(0, 0)
		SubMenu "Two-Color Imaging"
			"L-Clustering analysis...", /Q, DoClusteringAnalysis(1, 1)
			"Pairwise Correlation...", /Q, DoClusteringAnalysis(0, 1)
		End
	End
	SubMenu "Manipulate Particle Tracks"
		"Spatially filter tracks...", /Q, ManuallyFilterTracks()
		"-"
		"Save tracks to text file...", /Q, SaveTracksToTextFile()
	End
	SubMenu "Analyze Particle Tracks"
		"Calculate MSDs...", /Q, CalculateMSDs_menu()
		"Calculate CDFs...", /Q, CalculateCDFs_menu()
		"Displacement histograms...", /Q, DisplacementHistograms_menu()
	End
	SubMenu "Make Localization Images"
		"Generate scatter plot...", /Q, GenerateScatterPlot_Menu()
		"Generate accumulated image...", /Q, MakeAccumulatedImage_Menu()
		"Generate bitmap image...", /Q, MakePALMImage_Menu()
		SubMenu "3D Localization"
			"Generate Scatter plot...", /Q, Generate3DScatterPlot_Menu()
		End
	End
	SubMenu "Make Tracking Images"
		"Plot all tracks...", /Q, MakeParticleTrackingScatterPlot()
	End







// D. Dibble created functions-----------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
	SubMenu "Export to Single Particle Analysis"
		"Export Particle Data", /Q, ExportScatterDataMenu()
		//SubtractBackground()
		StrVarOrDefault("root:ExportedStr","(Subtract Background"), /Q, SubtractBackground()
		"Compare Background Channel", /Q, CompareBackgroundChannelImport()
		"Decimalize Imported Data", /Q, Decimalize()
	End




	SubMenu "Export Track Information"
		"Export Track Histogram", /Q, ExtractTrackHistogram()
		"Export Velocity Histograms", /Q, ExtractVelocityHistograms()
		"Export Mean Square Displacement Histograms", /Q, ExtractMSDHistograms()
	End



//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------







	"Save topmost movie as...", /Q, SaveCCDImagesAs()
	"Run batch localization...", /Q, SetupBatchLocalizationPanel()
End



Menu "Macros"
	"Append colorscale sliders", /Q, AppendColorRangeSliders()
	"RGB Merge.Print..", /Q, RunRGBMerge()
End

Menu "GraphMarquee"
	SubMenu "Localizer"
		"Average Intensity Trajectory", /Q, MakeAverageTraceFromMarquee()
		"Crop Images", /Q, MakeCroppedStackFromMarquee()
		SubMenu "Clustering"
			"Direct", /Q, Clustering_Marquee("Direct")
			"Fit Gaussian", /Q, Clustering_Marquee("MLEGauss")
			"L function", /Q, Clustering_Marquee("LClustering")
			"Pairwise correlation", /Q, Clustering_Marquee("PairwiseCorrelation")
		End
		SubMenu "Particle Tracking"
			"Calculate MSDs", /Q, CalculateMSDsCDFs_Marquee("msd")
			"Calculate CFDs", /Q, CalculateMSDsCDFs_Marquee("cdf")
		End
	End
End

Menu "TracePopup" 
	"-"
	SubMenu "Localizer"
		"Delete this localization", /Q, DeleteThisPos_TracePopup()
		SubMenu "Particle tracking"
			"Combine two tracks...", /Q, CombineTracks_TracePopup()
			"Split track here", /Q, SplitTracks_TracePopup()
			"-"
			"Delete this track position", /Q, DeleteTrackPosition_Popup()
			"-"
			"Delete this track", /Q, DeleteTrack_TracePopup()
		End
	End
End


Function SetupImageViewer_loadFromDisk()
	variable refNum
	string filePaths, currentFile






// D. Dibble created functions----------------------------------------------
//--------------------------------------------------------------------------
	string/g imageFileLocation
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------



	Open /D/R/T="????" /MULT=1 RefNum
	if (strlen(S_Filename) == 0)	// The user cancelled the dialog
		return -1
	endif
	filePaths = S_Filename
	variable nFiles = ItemsInList(filePaths, "\r")
	variable i
	
	for (i = 0; i < nFiles; i+=1)
		currentFile = StringFromList(i, filePaths, "\r")
		NewInterfacePanel(currentFile)
	endfor


// D. Dibble created functions----------------------------------------------
// **the assumption** is that we are only analyzing a single file-----------
	imageFileLocation = StringFromList(0, filePaths, "\r")
//--------------------------------------------------------------------------


End

Function SetupImageViewer_multiFileTIFF()
	variable refNum
	string filePath
	
	string fileFilters = "TIFF Files (*.tif,*.tiff,*.btf,*.tf8):.tif,.tiff,.btf,.tf8;"
	Open /D/R/F=fileFilters /M="Open one of the TIFF files in the sequence" RefNum
	if (strlen(S_Filename) == 0)	// The user cancelled the dialog
		return -1
	endif
	filePath = S_Filename
	
	NewInterfacePanel(filePath, fileType=CAMERA_TYPE_MULTIFILE_TIFF)
End

Function SetupImageViewer_IgorWaves()
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	
	// makeSetupImageViewer_loadFromDisk() a listbox window to choose from
	DoWindow /F candidateWavesViewer
	if (V_flag == 0)
		NewPanel /K=1/W=(101,286,401,486) /N=candidateWavesViewer as "Choose a matrix wave"
		ListBox ListBoxCandidateWaves,pos={1,2},size={298,183}
		MakeListIntoWaveSelector("candidateWavesViewer", "ListBoxCandidateWaves", content = WMWS_Waves, selectionMode = WMWS_SelectionNonContiguous,  nameFilterProc="WaveSelectorFilter2DAnd3DWaves")
		WS_SetNotificationProc("candidateWavesViewer", "ListBoxCandidateWaves", "wavesViewer_NotificationProc")
		ListBox ListBoxCandidateWaves,userdata(ResizeControlsInfo)= A"!!,<7!!#7a!!#BO!!#AFz!!#](Aon\"Qzzzzzzzzzzzzzz!!#o2B4uAezz"
		ListBox ListBoxCandidateWaves,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#u:Du]k<zzzzzzzzzzz"
		ListBox ListBoxCandidateWaves,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		SetWindow kwTopWin,hook(ResizeControls)=ResizeControls#ResizeControlsHook
		SetWindow kwTopWin,userdata(ResizeControlsInfo)= A"!!*'\"z!!#BP!!#AWzzzzzzzzzzzzzzzzzzzzz"
		SetWindow kwTopWin,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzzzzzzzzzzzzzzz"
		SetWindow kwTopWin,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzzzzzzzzz!!!"
		
		// if the window was opened previously then restore its size and position
		WC_WindowCoordinatesRestore("candidateWavesViewer")
		SetWindow candidateWavesViewer hook(WindowCoordsHook)=WC_WindowCoordinatesNamedHook
	endif
End

Function WaveSelectorFilter2DAnd3DWaves(aName, contents)
	String aName
	Variable contents
	
	if ((contents != WMWS_Waves) && (contents != WMWS_DataFolders))	// only accept waves or DFs
		return 0
	endif
	
	wave /Z proposedWave = $aName
	if (WaveExists(proposedWave))
		// only accept if the wave is 2D or 3D
		if ((DimSize(proposedWave, 1) == 0) || (DimSize(proposedWave, 3) != 0))
			return 0
		endif
	endif
	return 1
end

Function wavesViewer_NotificationProc(SelectedItem, EventCode)
	String selectedItem
	Variable eventCode
	
	switch (eventCode)
		case WMWS_DoubleClick:
			string wavePaths = WS_SelectedObjectsList("candidateWavesViewer", "ListBoxCandidateWaves")
			DoWindow /K candidateWavesViewer
			variable nWaves = ItemsInList(wavePaths)
			variable i
			
			for (i = 0; i < nWaves; i+=1)
				wave /Z dataWave = $StringFromList(i, wavePaths)
				if (WaveExists(dataWave))
					NewInterfacePanel(StringFromList(i, wavePaths))
				else
					Abort "Unable to find the selected wave!"
				endif
			endfor
			break
		default:
			return 0
	endswitch
End

Function OpenAverageImages()

	Open /D/R/T="????" /MULT=1 RefNum
	if (strlen(S_Filename) == 0)	// The user cancelled the dialog
		return -1
	endif
	string filePaths = S_Filename
	variable nFiles = ItemsInList(filePaths, "\r")
	
	// averages will be stored in this folder
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	NewDataFolder /O root:Packages:Localizer:AverageImages
	DFREF avgFolder = root:Packages:Localizer:AverageImages
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct progressStruct
	progressStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
	
	// calculate all required average images
	variable i
	string  thisFilePath, fileName, outputWaveName, macFilePath
	for (i = 0; i < nFiles; i+=1)
		if (nFiles > 1)
			LocalizerProgressWindow("File", i + 1, nFiles)
		endif
		
		thisFilePath = StringFromList(i, filePaths, "\r")
		macFilePath = ParseFilePath(5, thisFilePath, ":", 0, 1)
		fileName = ParseFilePath(3, macFilePath, ":", 0, 1)
		
		outputWaveName = GetValidWaveNameFromFileName(fileName, "avg")
		
		AnalyzeCCDImages /M=(ANALYZECCD_AVERAGE_IMAGE) /PROG=progressStruct /DEST=avgFolder:$outputWaveName thisFilePath
		wave M_AverageImage = avgFolder:$outputWaveName
		NewInterfacePanel("root:Packages:Localizer:AverageImages:" + outputWaveName)
	endfor
	
End

Function NewInterfacePanel(filePath, [fileType])
	String filePath
	variable fileType
	
	String windowName
	
	// check if the user already has this file open
	// in another file. If so then do nothing
	if (CheckIfPanelFromSameFileExists(filePath) != 0)
		return 0
	endif
	
	string fileName = ParseFilePath(5, filePath, ":", 0, 1)
	fileName = ParseFilePath(3, fileName, ":", 0, 1)
	
	// get the coordinates for the panel
	variable panelLeft, panelTop, panelRight, panelBottom
	GetNewInterfacePanelCoordinates(panelLeft, panelTop, panelRight, panelBottom)
	
	// get the name for the panel
	// make sure that no panel ever shares the same name
	// by keeping an incrementing counter
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	variable /G root:Packages:Localizer:V_nextWindowIndex
	NVAR nextWindowIndex = root:Packages:Localizer:V_nextWindowIndex
	
	// if the user is mixing windows created with an older version
	// then the name is not necessarily unique
	// so do some trickery for that
	do
		sprintf windowName, "LocalizerViewer%d", nextWindowIndex
		nextWindowIndex += 1
		DFREF testDF = root:Packages:Localizer:$windowName
	while (DataFolderRefStatus(testDF) != 0)
	
	// now make the panel itself
	NewPanel /K=1 /N=$windowName /W=(panelLeft,panelTop,panelRight,panelBottom) as fileName
	Assert(StringMatch(windowName, S_name) == 1)
	windowName = S_name
	
	// now that we know the name of the panel, set up all of the
	// other options we will use
	NewDataFolder /O root:Packages:Localizer:$windowName
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	DFREF packageDataFolder = root:Packages:Localizer
	
	// check if we need to force a particular file type (exists due to multifile tiff)
	if (ParamIsDefault(fileType))
		fileType = -1
	endif
	
	// set up a global that will keep track of the filetype
	variable /G windowDataFolder:V_fileType
	NVAR gFileType = windowDataFolder:V_fileType
	
	// now read the actual image
	ReadCCDImages /Y=(fileType) /C=1 /S=0 /DEST=windowDataFolder:M_CCDImage filePath
	if (V_flag != 0)
		DoWindow /K $windowName
		KillDataFolder /Z windowDataFolder
		DoAlert 0, "unable to read from " + filePath
		return 1
	endif
	gFileType = V_fileType
	
	// convert the image from mxnx1 (which can cause trouble) into mxn
	Redimension /N=(-1, -1) windowDataFolder:M_CCDImage
	
	string /G windowDataFolder:S_filePath = filePath
	variable /G windowDataFolder:V_nImages
	variable /G windowDataFolder:V_currentImage
	variable /G windowDataFolder:V_RightMouseH
	variable /G windowDataFolder:V_RightMouseV
	
	variable /G windowDataFolder:V_CCDPixelSize
	variable /G windowDataFolder:V_currentThreshold
	variable /G windowDataFolder:V_PFA
	variable /G windowDataFolder:V_smoothSigmaFactor
	variable /G windowDataFolder:V_gaussianWidth
	variable /G windowDataFolder:V_particleVerifierOverlap
	variable /G windowDataFolder:V_particleVerifierSymm
	variable /G windowDataFolder:V_particleVerifierEllipse
	variable /G windowDataFolder:V_particleVerifierEllipseAstig
	variable /G windowDataFolder:V_firstFrameInLocalization
	variable /G windowDataFolder:V_lastFrameInLocalization
	variable /G windowDataFolder:V_CCDPixelSize
	
	variable /G windowDataFolder:V_trackingMaxJump
	variable /G windowDataFolder:V_trackingMaxBlinking
	variable /G windowDataFolder:V_trackingMaxLength
	variable /G windowDataFolder:V_trackingMinLength
	variable /G windowDataFolder:V_trackingTimeDelta
	
	variable /G windowDataFolder:V_SOFInFramesToSkip
	variable /G windowDataFolder:V_SOFILastFrameInAnalysis
	variable /G windowDataFolder:V_SOFInFramesInMovieFrame
	variable /G windowDataFolder:V_SOFLagTime0
	variable /G windowDataFolder:V_SOFLagTime1
	variable /G windowDataFolder:V_SOFInDeconvolutionIterations
	
	variable /G windowDataFolder:V_bleachingMaxJump
	variable /G windowDataFolder:V_nFramesAverageSubtraction
	variable /G windowDataFolder:V_cameraMultiplication
	variable /G windowDataFolder:V_cameraOffset
	
	NVAR currentImage = windowDataFolder:V_currentImage
	NVAR nImages = windowDataFolder:V_nImages
	nImages = V_numberOfImages
	
	NVAR smoothSigmaFactor = windowDataFolder:V_smoothSigmaFactor
	NVAR PFA = windowDataFolder:V_PFA
	NVAR currentThreshold = windowDataFolder:V_currentThreshold
	NVAR gaussianWidth = windowDataFolder:V_gaussianWidth
	NVAR pixelSize = windowDataFolder:V_CCDPixelSize
	NVAR particleVerifierOverlap = windowDataFolder:V_particleVerifierOverlap
	NVAR particleVerifierSymm = windowDataFolder:V_particleVerifierSymm
	NVAR particleVerifierEllipse = windowDataFolder:V_particleVerifierEllipse
	NVAR particleVerifierEllipseAstig = windowDataFolder:V_particleVerifierEllipseAstig
	NVAR firstFrameInLocalization = windowDataFolder:V_firstFrameInLocalization
	NVAR lastFrameInLocalization = windowDataFolder:V_lastFrameInLocalization
	
	NVAR trackingMaxJump = windowDataFolder:V_trackingMaxJump
	NVAR trackingMaxBlinking = windowDataFolder:V_trackingMaxBlinking
	NVAR trackingMaxLength =  windowDataFolder:V_trackingMaxLength
	NVAR trackingMinLength =  windowDataFolder:V_trackingMinLength
	NVAR trackingTimeDelta = windowDataFolder:V_trackingTimeDelta
	
	NVAR SOFInFramesToSkip = windowDataFolder:V_SOFInFramesToSkip
	NVAR SOFIlastFrameInAnalysis = windowDataFolder:V_SOFILastFrameInAnalysis
	NVAR SOFInFramesInMovieFrame = windowDataFolder:V_SOFInFramesInMovieFrame
	NVAR SOFLagTime0 = windowDataFolder:V_SOFLagTime0
	NVAR SOFLagTime1 = windowDataFolder:V_SOFLagTime1
	NVAR SOFInDeconvolutionIterations = windowDataFolder:V_SOFInDeconvolutionIterations
	
	NVAR bleachingMaxJump = windowDataFolder:V_bleachingMaxJump
	NVAR nFramesAverageSubtraction = windowDataFolder:V_nFramesAverageSubtraction
	NVAR cameraMultiplication = windowDataFolder:V_cameraMultiplication
	NVAR cameraOffset = windowDataFolder:V_cameraOffset
	
	//  set some defaults
	pfa = 25
	smoothSigmaFactor = 5
	gaussianWidth = 1.6
	particleVerifierOverlap = 1
	lastFrameInLocalization = nImages - 1
	trackingMaxJump = 5
	trackingMaxLength = inf
	trackingMinLength = 2
	trackingTimeDelta = 0
	SOFInDeconvolutionIterations = 2
	bleachingMaxJump = 1
	cameraMultiplication = 1
	SOFInFramesToSkip = 0
	SOFIlastFrameInAnalysis = nImages - 1
	SOFInFramesInMovieFrame = 50
	SOFLagTime0 = 0
	SOFLagTime1 = 0
	
	string listParticleVerifiersFunc = GetIndependentModuleName() + "#PMParticleVerificationListItems()"
	string listPositionsFunc = GetIndependentModuleName() + "#GetPossiblePositionsWavesForPM()"
	string listTracksFunc = GetIndependentModuleName() + "#GetPossibleTrackingWaves()"
	string listCombinationWeightsFunc = GetIndependentModuleName() + "#GetPossibleSOFICombWeights()"
	
	// add the viewer graph
	DefineGuide /W=$windowName GViewerTop={FT,kViewerGraph_EdgeMargin},GViewerBottom={FB,-(kViewerControls_Margin + kViewerControls_Height + kViewerGraph_EdgeMargin)},GViewerLeft={FL,kViewerGraph_EdgeMargin}
	DefineGuide /W=$windowName GViewerRight={FR,-kViewerGraph_EdgeMargin}
	Display /N=CCDViewer /W=(182,133,546,401)/FG=(GViewerLeft,GViewerTop,GViewerRight,GViewerBottom)/HOST=$windowName
	
	// add the bottom controls in a subwindow
	DefineGuide /W=$windowName GLowerPanelBottom={FB,-kViewerControls_Margin}, GLowerPanelTop={GLowerPanelBottom, -kViewerControls_Height}
	DefineGuide /W=$windowName GLowerPanelLeft={GViewerLeft, 0}, GlowerPanelRight={GViewerRight, 0}
	
	NewPanel /N=CCDViewerControls /HOST=$windowName /FG=(GLowerPanelLeft, GLowerPanelTop, GLowerPanelRight, GLowerPanelBottom)
	SetDrawLayer UserBack
	Slider SLFrameSlider,win=$windowName#CCDViewerControls,pos={9,15},size={298,13},proc=SLViewerSliderProc
	Slider SLFrameSlider,win=$windowName#CCDViewerControls,limits={0,nImages - 1,1},value= 0,side= 0,vert= 0
	ValDisplay VDMaxFrames,win=$windowName#CCDViewerControls,pos={385,14},size={39,13},limits={0,0,0},barmisc={0,1000}
	ValDisplay VDMaxFrames,win=$windowName#CCDViewerControls,value= #("root:Packages:Localizer:" + windowName + ":V_nImages - 1")
	Button BTShowAnalysisTools,win=$windowName#CCDViewerControls,pos={429,10},size={30,20},proc=BTToggleAnalysisTools,title=">>"
	SetVariable SVCurrentFrame,win=$windowName#CCDViewerControls,pos={316,13},size={66,15}
	SetVariable SVCurrentFrame,win=$windowName#CCDViewerControls,limits={0,nImages -1 ,0},value=currentImage, proc=SVSetImage,title=" "
	
	// Guides for the panel containing the controls
	// the 'ControlFrame' guides are those that should be used by the actual controls to become available to the user (when the corresponding tab is selected)
	// they will be redefined appropriately by the tab procedure
	DefineGuide /W=$windowName GControlFrameLeft={FR,kAnalysisControls_Margin}, GControlFrameRight={FR,kAnalysisControls_Margin + kAnalysisControls_Width}
	DefineGuide /W=$windowName GControlFrameBottom={FB,-kAnalysisControls_Margin}, GControlFrameTop={GControlFrameBottom,-kAnalysisControls_Height}
	// the 'DisabledControl' guides are those that should be used to hide controls (when another tab is selected)
	DefineGuide /W=$windowName GDisabledControlsLeft={FR,kAnalysisControls_Margin}, GDisabledControlsRight={FR,kAnalysisControls_Margin + kAnalysisControls_Width}
	DefineGuide /W=$windowName GDisabledControlsBottom={FB,-kAnalysisControls_Margin}, GDisabledControlsTop={GControlFrameBottom,-kAnalysisControls_Height}
	DefineGuide /W=$windowName GTabBoxTop={GControlFrameBottom,-kAnalysisControls_Height - kTabHeight}
	NewPanel /N=TabBox /HOST=$windowName /FG=(GControlFrameLeft, GTabBoxTop, GControlFrameRight, GControlFrameBottom)
	TabControl TBSelectControls,win=$windowName#TabBox,pos={0,0},size={283,327},tabLabel(0)="Localization",tabLabel(1)="Tracking",tablabel(2)="SOFI",tabLabel(3)="Other",proc=TBSelectControlsProc
	
	// add the analysis controls (initially hidden)
	DefineGuide /W=$windowName GLocalizationControlsLeft={GControlFrameLeft,0}, GLocalizationControlsRight={GControlFrameRight,0}
	NewPanel /N=AnalysisControls /HOST=$windowName /FG=(GLocalizationControlsLeft, GControlFrameTop, GLocalizationControlsRight,GControlFrameBottom)
	SetDrawLayer UserBack
	SetDrawEnv fsize= 9
	DrawText 54,136,"Particle verification:"
	CheckBox CBThresholdCheck,win=$windowName#AnalysisControls,pos={9,11},size={62,14},proc=CBShowThresholdProc
	CheckBox CBThresholdCheck,win=$windowName#AnalysisControls,value= 0,title=""
	Button BTAnalyzePALM,win=$windowName#AnalysisControls,pos={199,275},size={70,20},proc=BTAnalyzePALMProc,title="Analyze!"
	CheckBox CBDisplayFittedPositions,win=$windowName#AnalysisControls,pos={9,173},size={120,14},proc=CBDisplayFittedPositionsProc,title="Display Fitted Positions"
	CheckBox CBDisplayFittedPositions,win=$windowName#AnalysisControls,value= 0
	PopupMenu PMSelectPositionsWave,win=$windowName#AnalysisControls,pos={142,171},size={125,20},bodyWidth=125,proc=PMUpdateViewerProc
	PopupMenu PMSelectPositionsWave,win=$windowName#AnalysisControls,mode=1,value=#listPositionsFunc
	SetVariable SVSetPixelSize,win=$windowName#AnalysisControls,pos={9,197},size={195,15},title="CCD pixel size (nm, optional):"
	SetVariable SVSetPixelSize,win=$windowName#AnalysisControls,value= pixelSize
	PopupMenu PMThresholdMethod,win=$windowName#AnalysisControls,pos={142,8},size={125,20},bodyWidth=125,proc=PMUpdateThresholdMethodProc
	PopupMenu PMThresholdMethod,win=$windowName#AnalysisControls,mode=1,value= #"\"GLRT;SmoothSigma;Absolute\"",title="Segmentation algorithm:"
	SetVariable SVSetGaussianWidthLocalization,win=$windowName#AnalysisControls,pos={9,214},size={260,15},title="Standard deviation of the PSF (pixels):"
	SetVariable SVSetGaussianWidthLocalization,win=$windowName#AnalysisControls,limits={0.1,inf,0.1},value= gaussianWidth,proc=SVUpdateThreshold
	SetVariable SVSetPFA,win=$windowName#AnalysisControls,pos={10,233},size={175,15},disable=0,title="GLRT Insensitivity:"
	SetVariable SVSetPFA,win=$windowName#AnalysisControls,limits={1,inf,1},value= pfa,proc=SVUpdateThreshold
	SetVariable SVSetSmoothSigmaFactor,win=$windowName#AnalysisControls,pos={10,233},size={175,15},disable=1,title="SmoothSigma factor:"
	SetVariable SVSetSmoothSigmaFactor,win=$windowName#AnalysisControls,limits={1e-3,inf,0.1},value= smoothSigmaFactor,proc=SVUpdateThreshold
	SetVariable SVSetThreshold,win=$windowName#AnalysisControls,pos={10,233},size={125,15},title="Threshold:"
	SetVariable SVSetThreshold,win=$windowName#AnalysisControls,value= currentThreshold,disable=1,proc=SVUpdateThreshold
	PopupMenu PMThresholdPreprocessing,win=$windowName#AnalysisControls,pos={27,49},size={240,20},bodyWidth=125,proc=PMUpdateViewerProc,title="Threshold preprocessing: "
	PopupMenu PMThresholdPreprocessing,win=$windowName#AnalysisControls,mode=1,popvalue="None",value= #"\"None;3x3 Median Filter;5x5 Median Filter;1x1 Gaussian Smoothing;2x2 Gaussian Smoothing;3x3 Mean Filter;5x5 Mean Filter\""
	PopupMenu PMThresholdPostProcessing,win=$windowName#AnalysisControls,pos={22,72},size={246,20},bodyWidth=125,proc=PMUpdateViewerProc,title="Threshold postprocessing: "
	PopupMenu PMThresholdPostProcessing,win=$windowName#AnalysisControls,mode=1,popvalue="None",value= #"\"None;Remove Isolated Pixels\""
	CheckBox CBShowParticles,win=$windowName#AnalysisControls,pos={27,31},size={128,14},proc=CBCheckShowParticlesProc,title="Show estimated particles"
	CheckBox CBShowParticles,win=$windowName#AnalysisControls,value= 0,disable=0
	PopupMenu PMParticleFinder,win=$windowName#AnalysisControls,pos={73,95},size={195,20},bodyWidth=125,proc=PMUpdateViewerProc,title="Particle finding:"
	PopupMenu PMParticleFinder,win=$windowName#AnalysisControls,mode=1,value= #"\"8-way adjacency;4-way adjacency\""
	PopupMenu PMParticleVerification,win=$windowName#AnalysisControls,pos={53,119},size={215,20},bodyWidth=125,proc=PMParticleVerificationProc,title="Choose..."
	PopupMenu PMParticleVerification,win=$windowName#AnalysisControls,mode=0,value=#listParticleVerifiersFunc
	SetVariable SVNFramesToSkip,win=$windowName#AnalysisControls,pos={15,254},size={130,15},bodyWidth=60,title="Frames between",value=firstFrameInLocalization
	SetVariable SVNFramesToSkip,win=$windowName#AnalysisControls,limits={0, nImages - 1, 1}
	SetVariable SVLastFrameToInclude,win=$windowName#AnalysisControls,pos={185,254},size={80,15},bodyWidth=60,title="and"
	SetVariable SVLastFrameToInclude,win=$windowName#AnalysisControls,value=lastFrameInLocalization,limits={0, nImages - 1, 1}
	PopupMenu PMLocalizationMethod,win=$windowName#AnalysisControls,pos={41,143},size={227,20},bodyWidth=125,title="Localization algorithm:",proc=PMUpdateViewerProc
	PopupMenu PMLocalizationMethod,win=$windowName#AnalysisControls,mode=1,popvalue="Gaussian Fitting",value= #"\"Gaussian Fitting;Gaussian Fitting (Fixed Width);Ellipsoidal Gaussian Fitting;Ellipsoidal Gaussian Fitting (astigmatism);Iterative Multiplication;Center-of-Mass;MLEwG\""
	Button BTCopyPositions,win=$windowName#AnalysisControls,pos={11,275},size={100,20},proc=BTCopyPositionsProc,title="Copy Positions"
	
	// add the particle tracking controls (initially inactive)
	DefineGuide /W=$windowName GTrackingControlsLeft={GDisabledControlsLeft,0}, GTrackingControlsRight={GDisabledControlsRight,0}
	NewPanel /N=TrackingControls /HOST=$windowName /FG=(GTrackingControlsLeft, GControlFrameTop, GTrackingControlsRight, GControlFrameBottom)
	GroupBox GBMakeTrack,win=$windowName#TrackingControls,pos={10,12},size={262,154},title="Particle tracking"
	GroupBox GBMakeTrack,win=$windowName#TrackingControls,frame=0
	SetVariable SVTrackingMaxJump,win=$windowName#TrackingControls,pos={20,57},size={178,15},bodyWidth=72,title="Max jump distance (px)"
	SetVariable SVTrackingMaxJump,win=$windowName#TrackingControls,limits={0,inf,0.1},value=trackingMaxJump
	SetVariable SVTrackingTimeDelta,win=$windowName#TrackingControls,pos={20,142},size={178,15},bodyWidth=72,title="Time between acq's (s)"
	SetVariable SVTrackingTimeDelta,win=$windowName#TrackingControls,limits={0,inf,0.1},value=trackingTimeDelta
	Button BTDoTracking,win=$windowName#TrackingControls,pos={209,140},size={55,20},proc=BTRunParticleTrackingProc,title="Do it"
	SetVariable SVTrackingMaxBlinking,win=$windowName#TrackingControls,pos={29,80},size={169,15},bodyWidth=72,title="Max blinking (frames)"
	SetVariable SVTrackingMaxBlinking,win=$windowName#TrackingControls,limits={0,inf,1},value=trackingMaxBlinking
	PopupMenu PMSelectPosWaveForTracking,win=$windowName#TrackingControls,pos={28,31},size={235,20},bodyWidth=150,title="Positions to track: "
	PopupMenu PMSelectPosWaveForTracking,win=$windowName#TrackingControls,mode=1,value= #listPositionsFunc
	SetVariable SVTrackMinTrackLength,win=$windowName#TrackingControls,pos={52,101},size={146,15},bodyWidth=72,title="Min track length"
	SetVariable SVTrackMinTrackLength,win=$windowName#TrackingControls,limits={2,inf,1},value=trackingMinLength
	SetVariable SVTrackMaxTrackLength,win=$windowName#TrackingControls,pos={50,122},size={148,15},bodyWidth=72,title="Max track length"
	SetVariable SVTrackMaxTrackLength,win=$windowName#TrackingControls,limits={0,inf,1},value=trackingMaxLength
	CheckBox CBDisplayTracks,win=$windowName#TrackingControls,pos={19,178},size={109,14},proc=CBDisplayFittedPositionsProc,title="Display Fitted Tracks"
	CheckBox CBDisplayTracks,win=$windowName#TrackingControls,value= 0
	CheckBox CBDisplayTrackIdentifier,win=$windowName#TrackingControls,pos={37,198},size={115,14},proc=CBDisplayFittedPositionsProc,title="Display track identifier"
	CheckBox CBDisplayTrackIdentifier,win=$windowName#TrackingControls,value=0,proc=CBDisplayFittedPositionsProc
	PopupMenu PMSelectTrackWave,win=$windowName#TrackingControls,pos={142,176},size={125,20},bodyWidth=125,proc=PMUpdateViewerProc
	PopupMenu PMSelectTrackWave,win=$windowName#TrackingControls,mode=1,value=#listTracksFunc
	PopupMenu PMDisplayTracksMode,win=$windowName#TrackingControls,pos={19,217},size={126,14},title="Mode", proc=PMUpdateViewerProc, value="Display tracks in entirety;Tracks are constructed with playback;Show all calculated tracks;"
	PopupMenu PMTracksColorMode,win=$windowName#TrackingControls,pos={19,242},size={210,14},title="Colors", proc=PMUpdateViewerProc, value="Arbitrary Track Colors;Color tracks by average displacement;Color each step individually by distance;"
	
	// add the SOFI controls (initially inactive)
	DefineGuide /W=$windowName GSOFIControlsLeft={GDisabledControlsLeft,0}, GSOFIControlsRight={GDisabledControlsRight,0}
	NewPanel /N=SOFIControls /HOST=$windowName /FG=(GSOFIControlsLeft, GControlFrameTop, GSOFIControlsRight, GControlFrameBottom)
	GroupBox GBSOFISettings,win=$windowName#SOFIControls,pos={10,12},size={262,195},title="SOFI Calculation"
	GroupBox GBSOFISettings,win=$windowName#SOFIControls,frame=0
	PopupMenu PMSOFIOrder,win=$windowName#SOFIControls,pos={21,31},size={120,20},title="Order: ",proc=PMSOFIOrderProc
	PopupMenu PMSOFIOrder,win=$windowName#SOFIControls,mode=1,value="2;3;4;5;6;"
	PopupMenu PMCalculationQuality,win=$windowName#SOFIControls,pos={109,32},size={89,20},title="Pixel Combinations:"
	PopupMenu PMCalculationQuality,win=$windowName#SOFIControls,mode=1,value= #"\"Minimum;Few;More;All;\""
	CheckBox RBAutoSOFI,win=$windowName#SOFIControls,pos={169,33},size={86,14},title="Autocumulant",disable=1
	CheckBox RBAutoSOFI,win=$windowName#SOFIControls,value= 0,mode=1,proc=RBSOFICorrelationModeProc
	CheckBox RBCrossSOFI,win=$windowName#SOFIControls,pos={169,55},size={88,14},title="Crosscumulant",disable=1
	CheckBox RBCrossSOFI,win=$windowName#SOFIControls,value= 1,mode=1,proc=RBSOFICorrelationModeProc
	CheckBox CBSOFIDoAverageImage,win=$windowName#SOFIControls,pos={21,56},size={103,14},title="Also average image"
	CheckBox CBSOFIDoAverageImage,win=$windowName#SOFIControls,value= 1
	SetVariable SVSOFINFramesToSkip,win=$windowName#SOFIControls,pos={19,78},size={136,15},bodyWidth=60,title="Frames between"
	SetVariable SVSOFINFramesToSkip,win=$windowName#SOFIControls,value=SOFInFramesToSkip,limits={0, nImages - 1, 1}
	SetVariable SVSOFILastFrameToInclude,win=$windowName#SOFIControls,pos={185,78},size={80,15},bodyWidth=60,title="and"
	SetVariable SVSOFILastFrameToInclude,win=$windowName#SOFIControls,value=SOFIlastFrameInAnalysis,limits={0, nImages - 1, 1}
	CheckBox CBSOFIDropSaturatedFrames,win=$windowName#SOFIControls,pos={19,101},size={203,14},size={203,14},title="Drop frames with saturated pixels (spikes)"
	CheckBox CBSOFIDropSaturatedFrames,win=$windowName#SOFIControls,value=0, disable=1
	CheckBox CBSOFIMakeMovie,win=$windowName#SOFIControls,pos={19,123},size={72,14},title="Make movie:"
	CheckBox CBSOFIMakeMovie,win=$windowName#SOFIControls,value= 0
	CheckBox CBUseCombinationWeights,win=$windowName#SOFIControls,pos={19.00,145.00},size={128.00,16.00},title="Combination weights"
	CheckBox CBUseCombinationWeights,win=$windowName#SOFIControls,value= 0
	PopupMenu PMCombinationWeights,win=$windowName#SOFIControls,pos={153.00,143.00},size={61.00,23.00}
	PopupMenu PMCombinationWeights,win=$windowName#SOFIControls,mode=1,value= #listCombinationWeightsFunc
	Button BTCalculateCombinationWeights,pos={36.00,160.00},size={100.00,20.00},proc=BTCalcCombinationWeights,title="Calc weights"
	SetVariable SVSOFInFramesInMovie,win=$windowName#SOFIControls,pos={97,123},size={168,15},bodyWidth=60,title=" frames for every image"
	SetVariable SVSOFInFramesInMovie,win=$windowName#SOFIControls,limits={25, nImages,1},value=SOFInFramesInMovieFrame
	SetVariable SVSOFILagTime0,win=$windowName#SOFIControls,pos={21,143},size={125,15},title="Lag time 0: ",value=SOFLagTime0
	SetVariable SVSOFILagTime0,win=$windowName#SOFIControls,limits={0,100,1},disable=1
	SetVariable SVSOFILagTime1,win=$windowName#SOFIControls,pos={21,163},size={125,15},title="Lag time 1: ",value=SOFLagTime1
	SetVariable SVSOFILagTime1,win=$windowName#SOFIControls,limits={0,100,1},disable=1
	Button BTDoSOFI,win=$windowName#SOFIControls,pos={209,184},size={55,20},proc=BTDoSOFIProc,title="Do it"
	GroupBox GBSOFIDeconvolution,win=$windowName#SOFIControls,pos={10,215},size={262,90},title="Richardson-Lucy Deconvolution"
	GroupBox GBSOFIDeconvolution,win=$windowName#SOFIControls,frame=0
	SetVariable SVSetGaussianWidthSOFI,win=$windowName#SOFIControls,pos={16,233},size={230,15},title="Standard deviation of the PSF (pixels):"
	SetVariable SVSetGaussianWidthSOFI,win=$windowName#SOFIControls,limits={0.1,inf,0.1},value=gaussianWidth
	SetVariable SVDeconvolutionIterations,win=$windowName#SOFIControls,pos={16,254},size={230,15},title="Number of iterations:"
	SetVariable SVDeconvolutionIterations,win=$windowName#SOFIControls,limits={1,inf,1},value=SOFInDeconvolutionIterations
	Button BTDoSOFIDeconvolution,win=$windowName#SOFIControls,pos={209,278},size={55,20},proc=BTDoSOFIDeconvolutionProc,title="Do it"
	
	// add the 'other' analysis controls (initially inactive)
	DefineGuide /W=$windowName GOtherControlsLeft={GDisabledControlsLeft,0}, GOtherControlsRight={GDisabledControlsRight,0}
	NewPanel /N=OtherControls /HOST=$windowName /FG=(GOtherControlsLeft, GControlFrameTop, GOtherControlsRight, GControlFrameBottom)
	SetDrawLayer UserBack
	SetDrawEnv fsize= 9,fstyle= 2
	DrawText 20,132,"Other settings will be taken from 'Localization' tab"
	GroupBox GBBasicAnalysis,win=$windowName#OtherControls, pos={10,12},size={262,54},title="Basic analysis",frame=0
	PopupMenu PMChooseCCDAnalysisMethod,win=$windowName#OtherControls,pos={20,36},size={175,20},bodyWidth=175
	PopupMenu PMChooseCCDAnalysisMethod,win=$windowName#OtherControls,mode=1,value= #"\"Average Intensity Trace;Summed Intensity Trace;Average Image;Variance Image;\""
	Button BTDoBasicAnalysis,win=$windowName#OtherControls,pos={207,36},size={55,20},title="Do it",proc=BTBasicAnalysisProc
	GroupBox GBBleaching,win=$windowName#OtherControls,pos={10,72},size={262,72},title="Bleaching analysis",frame=0
	SetVariable SVBleachingMaxJump,win=$windowName#OtherControls,pos={20,96},size={178,15},limits={0,100,0.1},bodyWidth=72,title="Max jump distance (px)",value=bleachingMaxJump
	Button BTDoBleaching,win=$windowName#OtherControls,pos={207,93},size={55,20},title="Do it",proc=BTBleachingAnalysisProc
	GroupBox GBDriftCorrection,win=$windowName#OtherControls,pos={10,152},size={262,72},title="Drift correction",frame=0
	PopupMenu PMDriftCorrectionMode,win=$windowName#OtherControls,pos={19,173},size={118,20},title="Mode:",mode=1,value= #"\"subimages;ad-hoc fiducials;manual fiducials;average position;\""
	PopupMenu PMDriftCorrSelectPositionsWave,win=$windowName#OtherControls,pos={19,199},size={244,20},bodyWidth=150,title="Positions to correct: "
	PopupMenu PMDriftCorrSelectPositionsWave,win=$windowName#OtherControls,mode=1,value= #listPositionsFunc
	Button BTDoDriftCorrection,win=$windowName#OtherControls,pos={207,173},size={55,20},proc=BTCorrectDriftProc,title="Do it"
	GroupBox GBMovieProcessing,win=$windowName#OtherControls, pos={10,232},size={262,72},title="Process movie",frame=0
	PopupMenu PMChooseMovieProcessingMethod,win=$windowName#OtherControls,pos={20,254},size={175,20},bodyWidth=175,proc=PMChooseMovieProcessingProc
	PopupMenu PMChooseMovieProcessingMethod,win=$windowName#OtherControls,mode=1,value= #"\"Subtract Pixelwise Average;Difference Image;Convert File Format;Convert to Photons\""
	SetVariable SVFramesAverageSubtraction,win=$windowName#OtherControls,pos={19,281},size={179,15},bodyWidth=72,title="Frames to average over"
	SetVariable SVFramesAverageSubtraction,win=$windowName#OtherControls,limits={0,inf,1},value=nFramesAverageSubtraction
	SetVariable SVCameraMultiplication,win=$windowName#OtherControls,pos={146,285},size={121,15},bodyWidth=60,title="Multiplication"
	SetVariable SVCameraMultiplication,win=$windowName#OtherControls,limits={0,100,0.1},value=cameraMultiplication,disable=1
	SetVariable SVCameraOffset,win=$windowName#OtherControls,pos={20,285},size={91,15},bodyWidth=60,title="Offset"
	SetVariable SVCameraOffset,win=$windowName#OtherControls,limits={0,100,0.1},value=cameraOffset,disable=1
	Button BTDoMovieProcessing,win=$windowName#OtherControls,pos={207,254},size={55,20},title="Do it",proc=BTProcessMovieProc
	
	// add a panel for the LUT controls
	DefineGuide /W=$windowName GLUTControlsLeft={GControlFrameLeft, 0}, GLUTControlsRight={GControlFrameRight, 0}
	DefineGuide /W=$windowName GLUTControlsBottom={GTabBoxTop, -kHistogram_Margin}, GLUTControlsTop={GLUTControlsBottom, -kLUTControls_Height}
	NewPanel /N=LUTControls /HOST=$windowName /FG=(GLUTControlsLeft, GLUTControlsTop, GLUTControlsRight, GLUTControlsBottom)
	PopupMenu PMLUTColorTable,win=$windowName#LUTControls,pos={12,13},size={181,20},bodyWidth=125,title="Color table: "
	PopupMenu PMLUTColorTable,win=$windowName#LUTControls,mode=1,value= #"\"*COLORTABLEPOPNONAMES*\"",proc=PMLUTChangeColorTableProc
	CheckBox CBReverseColorTable,win=$windowName#LUTControls,pos={208,15},size={52,14},title="Reverse",value= 0,proc=CBLUTProc
	GroupBox GBMinColorScale,win=$windowName#LUTControls,pos={11,39},size={135,61},title="LUT min"
	GroupBox GBMaxColorScale,win=$windowName#LUTControls,pos={149,39},size={126,61},title="LUT max"
	CheckBox CBAutoMax,win=$windowName#LUTControls,pos={153,58},size={39,14},title="Auto",value= 1,mode=1,proc=CBLUTProc
	CheckBox CBAutoMin,win=$windowName#LUTControls,pos={15,57},size={39,14},title="Auto",value= 1,mode=1,proc=CBLUTProc
	CheckBox CBManMin,win=$windowName#LUTControls,pos={15,77},size={50,14},title="Manual:",value= 0,mode=1,proc=CBLUTProc
	CheckBox CBManMax,win=$windowName#LUTControls,pos={153,77},size={50,14},title="Manual:",value= 0,mode=1,proc=CBLUTProc
	SetVariable SVLUTMinVal,win=$windowName#LUTControls,pos={71,77},size={50,15},proc=SVLUTProc
	SetVariable SVLUTMinVal,win=$windowName#LUTControls,limits={-inf,inf,0},value= _NUM:0
	SetVariable SVLUTMaxVal,win=$windowName#LUTControls,pos={207,77},size={50,15},proc=SVLUTProc
	SetVariable SVLUTMaxVal,win=$windowName#LUTControls,limits={-inf,inf,0},value= _NUM:0
	
	// add a graph for showing the histogram
	DefineGuide /W=$windowName GHistogramLeft={GControlFrameLeft, 0}, GHistogramRight={GControlFrameRight, 0}
	DefineGuide /W=$windowName GHistogramTop={GViewerTop, 0}, GHistogramBottom={GLUTControlsTop, -kHistogram_Margin}
	Display /HOST=$windowName /N=HistogramViewer /FG=(GHistogramLeft, GHistogramTop, GHistogramRight, GHistogramBottom)
	
	// append the CCD image to the viewer
	AppendImage /W=$windowName#CCDViewer windowDataFolder:M_CCDImage
	AdjustGraphForImageDisplay(windowName + "#CCDViewer")
	
	// also set up the histogram
	wave hist = MakeImageHistogram(windowDataFolder:M_CCDImage, windowName)
	AppendToGraph /W=$windowName#HistogramViewer hist
	Label /W=$windowName#HistogramViewer Bottom, "Intensity"
	Label /W=$windowName#HistogramViewer Left, "Occurrence"
	
	SetWindow $windowName,hook(keycontrol)=ViewerKeyHook
	SetWindow $windowName,hook(LUTControls)=LUTWindowHook
	SetWindow $windowName,hook(RightClickPos)=SaveRightClickCoordinatesHook
End

Function AdjustGraphForImageDisplay(graphNameStr)
	String graphNameStr
	
	Assert(strlen(graphNameStr) != 0)
	
	ModifyGraph /W=$graphNameStr margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14
	ModifyGraph /W=$graphNameStr mirror=2
	ModifyGraph /W=$graphNameStr nticks=4
	ModifyGraph /W=$graphNameStr minor=1
	ModifyGraph /W=$graphNameStr fSize=9
	ModifyGraph /W=$graphNameStr standoff=0
	ModifyGraph /W=$graphNameStr tkLblRot(left)=90
	ModifyGraph /W=$graphNameStr btLen=3
	ModifyGraph /W=$graphNameStr tlOffset=-2
	ModifyGraph /W=$graphNameStr tickUnit=1
End

Function CheckIfPanelFromSameFileExists(filePath)
	string filePath
	
	// called when a new panel is about to be created
	// check if the file requested by the user is not
	// already open in a different window
	
	// return 0 if there is no such window
	// otherwise bring it to the front and return 1
	variable i
	string windowName
	DFREF windowDataFolder
	for (i = 0; ; i+=1)
		windowName = WinName(i, 64, 1)
		if (strlen(windowName) == 0)
			return 0
		endif
		
		if (GrepString(windowName, "LocalizerViewer[0-9]*") == 1)
			// this is a window created by this package
			// check if the filePath is the same as the one
			// the user requests
			windowDataFolder = GetWindowDataFolder(windowName)
			SVAR windowFilePath = windowDataFolder:S_filePath
			
			if (StringMatch(filePath, windowFilePath) == 1)
				DoWindow /F $windowName
				return 1
			endif
		endif
	endfor
End

Function GetNewInterfacePanelCoordinates(left, top, right, bottom)
	variable &left, &top, &right, &bottom
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	
	variable /G root:Packages:Localizer:V_ViewerLeftCoordinate
	variable /G root:Packages:Localizer:V_ViewerTopCoordinate
	
	NVAR previousLeft = root:Packages:Localizer:V_ViewerLeftCoordinate
	NVAR previousTop = root:Packages:Localizer:V_ViewerTopCoordinate
	
	variable panelWidth = 595
	variable panelHeight = 590
	variable dx = 30, dy = 30	// offset new window
	
	string screenInfo = StringByKey("SCREEN1", IgorInfo(0))
	variable dummy, displayWidth, displayHeight
	sscanf screenInfo, "DEPTH=%g,RECT=0,0,%g,%g", dummy, displayWidth, displayHeight
	
	if (previousLeft == 0)
		left = 150
		top = 50
	else
		left = previousLeft + dx
		top = previousTop + dy
		if (left + panelWidth > displayWidth)
			// wrap around to the left side of the screen
			left = 150
			top = 50
		endif
		top = (top + panelHeight > displayHeight) ? 50 : top
	endif
	
	previousLeft = left
	previousTop = top
	right = left + panelWidth
	bottom = top + panelHeight
End

Function /DF GetWindowDataFolder(windowName)
	string windowName
	
	// window names can sometimes consist of paths to subwindows
	// e.g. when returned in a control procedure
	// this function takes care of all of that
	string baseWindowName = StringFromList(0, windowName, "#")
	return root:Packages:Localizer:$baseWindowName
End

Function /S GetBaseWindowName(windowName)
	string windowName
	
	// sometimes the window names returned can include subwindows
	// however, I take these into account explicitly
	// so I really want the base name
	return StringFromList(0, windowName, "#")
End

Function /S GetNameOfTopInterfacePanel()
	
	string windowName
	variable i
	
	for (i = 0; ; i+=1)
		windowName = WinName(i, 64, 1)
		if (strlen(windowName) == 0)
			return ""
		endif
		
		if (GrepString(windowName, "LocalizerViewer[0-9]*") == 1)
			return windowName
		endif
	endfor
End

Function ReadCCDFramesIntoMemory()
	
	variable RefNum, nFilesToLoad, i
	variable skipThisFile
	string filePath, filePaths, MacintoshFilePath, fileName
	string wName
	string errorMessage
	
	Open /D/R/MULT=1/T="????" RefNum

	if (strlen(S_Filename) == 0)	// The user cancelled the dialog
		return -1
	endif
	
	filePaths = S_FileName
	nFilesToLoad = ItemsInList(filePaths, "\r")
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct progressStruct
	progressStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
	
	for (i = 0; i < nFilesToLoad; i+=1)
		skipThisFile = 0
		filePath = StringFromList(i, filePaths, "\r")
		MacintoshFilePath = ParseFilePath(5, FilePath, ":", 0, 0)
		
		// get the number of images in the file
		ReadCCDImages /Z /H MacintoshFilePath
		if (V_flag != 0)		// did we fail to read the data?
			errorMessage = "Unable to read from  " + MacintoshFilePath
			Abort errorMessage
		endif
		
		wName = ParseFilePath(3, MacintoshFilePath, ":", 0, 0)
		wName = CleanupName(wName, 0)
		for ( ; ; )
			if (WaveExists(root:'Localizer Movies':$wName))
				wName = GetNewOutputWaveName(wName, root:'Localizer Movies', "Name already exists, please provide a different one")
				if (strlen(wName) == 0)
					skipThisFile = 1
					nFilesToLoad -= 1
				endif
				break
			else
				break
			endif
		endfor
		
		if (skipThisFile == 1)
			continue
		endif
		
		NewDataFolder /O root:'Localizer Movies'
		ReadCCDImages /O /PROG=progressStruct /DEST=root:'Localizer Movies':$wName MacintoshFilePath
		
		NewInterfacePanel("root:'Localizer Movies':" + wName)
	endfor
End

Function ChangeCurrentImage(n, windowName)
	variable n
	string windowName
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	NVAR nImages = windowDataFolder:V_nImages
	SVAR filePath = windowDataFolder:S_filePath
	NVAR currentImage = windowDataFolder:V_currentImage
	
	variable fileType
	NVAR /Z gFileType = windowDataFolder:V_fileType
	if (!NVAR_Exists(gFileType))
		fileType = -1
	else
		fileType = gFileType
	endif
	
	string baseWindowName = GetBaseWindowName(windowName)
	
	if (n < 0)
		n = 0
	endif
	
	if (n >= nImages)
		n = nImages - 1
	endif
	
	ReadCCDImages /Y=(fileType) /O /C=1 /S=(n) /DEST=windowDataFolder:M_CCDImage filePath
	currentImage = n
	
	// the XOP returns a mxnx1 wave, which causes errors when some expect mxn
	// so to handle that, do a redimension here
	Redimension /N=(-1, -1) windowDataFolder:M_CCDImage
	
	GenerateGUIThreshold(baseWindowName)
	
	// update the slider and the histogram
	Slider SLFrameSlider win=$baseWindowName#CCDViewerControls, value=currentImage
	MakeImageHistogram(windowDataFolder:M_CCDImage, windowName)
End

Function ChangeToNextImageAndLoop(windowName)
	string windowName
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	NVAR nImages = windowDataFolder:V_nImages
	NVAR currentImage = windowDataFolder:V_currentImage
	
	string baseWindowName = GetBaseWindowName(windowName)
	
	variable nextImage = (currentImage + 1 >= nImages) ? 0 : currentImage + 1
	ChangeCurrentImage(nextImage, baseWindowName)
End

Function BTViewerProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	string controlName = ba.ctrlName
	string windowName = ba.win
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	NVAR nImages = windowDataFolder:V_nImages
	NVAR currentImage = windowDataFolder:V_currentImage
	
	switch (ba.eventCode)
		case 1: 	// mouse down	
			strswitch(controlName)
			case "BTNext":
				ChangecurrentImage(currentImage + 1, windowName)
				break
			case "BTPrevious":
				ChangecurrentImage(currentImage - 1, windowName)
				break
			case "BTFirst":
				ChangecurrentImage(0, windowName)
				break
				
			case "BTLast":
				ChangecurrentImage(nImages - 1, windowName)
				break
			endswitch
		break
	EndSwitch
End

Function SLViewerSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	
	variable curval
	string windowName = sa.win
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	NVAR currentImage = windowDataFolder:V_currentImage

	switch( sa.eventCode )
		case -1: // kill
			break
		default:
			if( sa.eventCode & 1 ) // value set
				curval = sa.curval
				ChangecurrentImage(curval, windowName)
			endif
			break
	endswitch

	return 0
End

Function ViewerKeyHook(s)
	STRUCT WMWinHookStruct &s
	
	string baseWindowName = GetBaseWindowName(s.winName)
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	NVAR currentImage = windowDataFolder:V_currentImage
	
	variable didSomething = 0
	switch (s.eventcode)
		case 11:	// keyboard event
			switch (s.keyCode)
				case 29:	// right arrow
					ChangeCurrentImage(currentImage + 1, baseWindowName)
					didSomething = 1
					break
				case 28:	// left arrow
					ChangeCurrentImage(currentImage -1, baseWindowName)
					didSomething = 1
					break
				case 32:	// space
					ToggleMoviePlayback(baseWindowName)
					didSomething = 1
					break
			EndSwitch
			break
		case 22:	// scrollwheel event
			variable framesToScroll = - s.wheelDy
			ChangeCurrentImage(currentImage + framesToScroll, baseWindowName)
			didSomething = 1
			break
	EndSwitch
	
	return didSomething
End

Function LUTWindowHook(s)
	STRUCT WMWinHookStruct &s
	
	string baseWindowName = GetBaseWindowName(s.winName)
	
	GetWindow $s.winName activeSW
	String activeSubwindow = S_value
	if (StringMatch(activeSubwindow, "*CCDViewer") != 1)
		return 0
	endif
	
	Switch (s.eventcode)
		case 5:	// mouseup
			//Debugger
			break
		case 8:	// modified
			UpdateLUTSetttings(baseWindowName, 1)
			break
	EndSwitch
	
	return 0
End

Function SaveRightClickCoordinatesHook(s)
	STRUCT WMWinHookStruct &s
	
	string baseWindowName = GetBaseWindowName(s.winName)
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	NVAR currentImage = windowDataFolder:V_currentImage
	
	variable didSomething = 0
	switch (s.eventcode)
		case 3:	// mouse down
			variable isRightMouse = (s.eventMod & 2^4) != 0
			if (!isRightMouse)
				return 0
			endif
			NVAR rightMouseH = windowDataFolder:V_RightMouseH
			NVAR rightMouseV = windowDataFolder:V_RightMouseV
			rightMouseH = s.mouseLoc.h
			rightMouseV = s.mouseLoc.v
			break
	EndSwitch
	
	return 0
End

Function BTToggleAnalysisTools(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	string windowName = ba.win
	string controlName = ba.ctrlName
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	
	Variable /G windowDataFolder:V_ControlsAreShown
	NVAR controlsAreShown = windowDataFolder:V_ControlsAreShown
	
	switch (ba.eventCode)
		case 1: 	// mouse down
			if (!controlsAreShown)
				ShowAnalysisControls(windowName)
				Button $controlName, win=$windowName, title="<<"
				controlsAreShown = 1
			else
				HideAnalysisControls(windowName)
				Button $controlName, win=$windowName, title=">>"
				controlsAreShown = 0
			endif
			break
	EndSwitch
	
End

Function ShowAnalysisControls(wName)
	string wName
	
	string windowName = StringFromList(0, wName, "#")
	
	GetWindow $windowName, wsize
	
	variable widthOfControls = kAnalysisControls_Width + kAnalysisControls_Margin
	variable newWidth = (V_Right - V_Left + widthOfControls* 72 / ScreenResolution)
	
	// now redefine the guides
	// due a bug in the way Igor saves recreation macros, GViewerRight must be defined in terms of 
	// FR. The expression below was more clear, but triggers an error upon experiment load
	//DefineGuide /W=$windowName GViewerRight={GControlFrameLeft, -kViewerGraph_EdgeMargin}
	DefineGuide /W=$windowName GViewerRight={FR, - (kAnalysisControls_Margin + kAnalysisControls_Width + kViewerGraph_EdgeMargin)}
	DefineGuide /W=$windowName GControlFrameLeft={FR,- (kAnalysisControls_Margin + kAnalysisControls_Width)}, GControlFrameRight={FR,-kAnalysisControls_Margin}
	
	
	MoveWindow /W=$windowName V_Left, V_Top, V_Left + newWidth, V_Bottom
End

Function HideAnalysisControls(wName)
	string wName
	
	string windowName = StringFromList(0, wName, "#")
	
	GetWindow $windowName, wsize
	
	variable widthOfControls = kAnalysisControls_Width + kAnalysisControls_Margin
	variable newWidth = (V_Right - V_Left - widthOfControls * 72 / ScreenResolution)
	
	// now redefine the guides
	DefineGuide /W=$windowName GViewerRight={FR,-kViewerGraph_EdgeMargin}
	DefineGuide /W=$windowName GControlFrameLeft={FR,kAnalysisControls_Margin}, GControlFrameRight={FR,kAnalysisControls_Margin + kAnalysisControls_Width}
	
	MoveWindow /W=$windowName V_Left, V_Top, V_Left + newWidth, V_Bottom
End

Function TBSelectControlsProc(tca) : TabControl
	STRUCT WMTabControlAction &tca
	
	// handle the user selection of different controls using the tabbox
	// depending on the selection, move the unwanten controls outside the panel
	// and the wanted ones inside

	switch( tca.eventCode )
		case 2: // mouse up
			Variable tab = tca.tab
			string baseWindowName = GetBaseWindowName(tca.win)
			// get the text for the currently selected tab item
			ControlInfo /W=$baseWindowName#TabBox TBSelectControls
			string tabText = S_Value
			
			// start by deactivating all panels, then show the one that is requested
			DefineGuide /W=$baseWindowName GLocalizationControlsLeft={GDisabledControlsLeft,0}, GLocalizationControlsRight={GDisabledControlsRight,0}
			DefineGuide /W=$baseWindowName GTrackingControlsLeft={GDisabledControlsLeft,0}, GTrackingControlsRight={GDisabledControlsRight,0}
			DefineGuide /W=$baseWindowName GSOFIControlsLeft={GDisabledControlsLeft,0}, GSOFIControlsRight={GDisabledControlsRight,0}
			DefineGuide /W=$baseWindowName GOtherControlsLeft={GDisabledControlsLeft,0}, GOtherControlsRight={GDisabledControlsRight,0}
			
			StrSwitch(tabText)
				case "Localization":
					DefineGuide /W=$baseWindowName GLocalizationControlsLeft={GControlFrameLeft,0}, GLocalizationControlsRight={GControlFrameRight,0}
					break
				case "Tracking":
					DefineGuide /W=$baseWindowName GTrackingControlsLeft={GControlFrameLeft,0}, GTrackingControlsRight={GControlFrameRight,0}
					break
				case "SOFI":
					DefineGuide /W=$baseWindowName GSOFIControlsLeft={GControlFrameLeft,0}, GSOFIControlsRight={GControlFrameRight,0}
					break
				case "Other":
					DefineGuide /W=$baseWindowName GOtherControlsLeft={GControlFrameLeft,0}, GOtherControlsRight={GControlFrameRight,0}
					break
				default:
					Abort "Missing handler in TBSelectControlsProc, got " + tabText
					break
			EndSwitch
			break
	endswitch

	return 0
End

Function ToggleMoviePlayback(windowName)
	string windowName
	// see comments below
	
	string BaseWindowName = GetBaseWindowName(windowName)
	String /G root:Packages:Localizer:S_moviesToPlay
	SVAR moviesToPlay = root:Packages:Localizer:S_moviesToPlay
	
	if (FindListItem(baseWindowName, moviesToPlay) == -1)
		// start the movie
		StartPlayingMovie(baseWindowName)
	else
		//stop the movie
		StopPlayingMovie(baseWindowName)
	endif
End

Function StartPlayingMovie(windowName)
	string windowName
	
	// start playing the movie contained in windowName in a loop
	// the strategy used is to add a background task that updates the
	// frame each time it is run
	// a limitation is that the function executed by the task can not take
	// any parameters, so we have kludge our way around that so the
	// task knows which window needs to be updated
	// for now this is done by maintaining a global string variable
	// in the root:Packages:localizer folder, that contains a list
	// of the windows to be run in a movie
	// if this list is empty then the background task will terminate,
	// though manual termination is preferred
	
	string baseWindowName = GetBaseWindowName(windowName)
	// add this windowName to the string variable
	// if it already is in there then that is an error
	String /G root:Packages:Localizer:S_moviesToPlay
	SVAR moviesToPlay = root:Packages:Localizer:S_moviesToPlay
	if (FindListItem(baseWindowName, moviesToPlay) != -1)
		Abort "StopPlayingMovie error: the requested window is already present in S_moviesToPlay"
	endif
	
	moviesToPlay += baseWindowName + ";"
	
	// if there was a window already in the string then the task will be running
	// if not then start it
	if (ItemsInList(moviesToPlay) == 1)
		StartBackGroundMovieTask()
	endif
End

Function StopPlayingMovie(windowName)
	string windowName
	
	SVAR moviesToPlay = root:Packages:Localizer:S_moviesToPlay
	string baseWindowName = GetBaseWindowName(windowName)
	
	// windowName must be in the list, otherwise it's an error
	if (FindListItem(baseWindowName, moviesToPlay) == -1)
		Abort "StopPlayingMovie error: the requested window is not present in S_moviesToPlay"
	endif
	
	moviesToPlay = RemoveFromList(baseWindowName, moviesToPlay)
	
	// if there are no more windows left then stop the task
	if (ItemsInList(moviesToPlay) == 0)
		StopBackGroundMovieTask()
	endif
End

Function StartBackGroundMovieTask()
	SVAR moviesToPlay = root:Packages:Localizer:S_moviesToPlay
	
	Assert(ItemsInList(moviesToPlay) > 0)
	
	variable numTicks = 4 // period in ticks, or 1/60th of second
	CtrlNamedBackground MoviePlayer, period=numTicks, proc=RunBackgroundMovie
	CtrlNamedBackground MoviePlayer, start
End

Function StopBackGroundMovieTask()
	CtrlNamedBackground MoviePlayer, stop
End

Function RunBackgroundMovie(s)
	// this function is called as a background task
	STRUCT WMBackgroundStruct &s
	
	SVAR moviesToPlay = root:Packages:Localizer:S_moviesToPlay
	
	variable nMoviesToPlay = ItemsInList(moviesToPlay)
	if (nMoviesToPlay == 0)
		// no movies to play, stop the task
		return 1
	endif
	
	variable i
	string currentBaseWindowName
	DFREF windowDataFolder
	for (i = 0; i < nMoviesToPlay; i+=1)
		currentBaseWindowName = GetBaseWindowName(StringFromList(i, moviesToPlay))
		try
			ChangeToNextImageAndLoop(currentBaseWindowName)
		catch
			return 1
		endtry
	endfor
	
	return 0
End

Function PMLUTChangeColorTableProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			string baseWindowName = GetBaseWindowName(pa.win)
			UpdateLUTSetttings(baseWindowName, 0)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CBLUTProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string baseWindowName = GetBaseWindowName(cba.win)
			
			if (!StringMatch(cba.ctrlName, "CBReverseColorTable"))		
				StrSwitch (cba.ctrlName)
					case "CBAutoMin":
						CheckBox CBManMin,win=$cba.win,value=0
						break
					case "CBManMin":
						CheckBox CBAutoMin, win=$cba.win,value=0
						break
					case "CBAutoMax":
						CheckBox CBManMax,win=$cba.win,value=0
						break
					case "CBManMax":
						CheckBox CBAutoMax,win=$cba.win,value=0
						break
					default:
						Abort "Unknown control in CBLUTProc()"
						break
				EndSwitch
			endif
			
			UpdateLUTSetttings(baseWindowName, 0)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function SVLUTProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			string baseWindowName = GetBaseWindowName(sva.win)
			
			StrSwitch(sva.ctrlName)
				case "SVLUTMinVal":
					CheckBox CBManMin,win=$sva.win,value=1
					CheckBox CBAutoMin,win=$sva.win,value=0
					break
				case "SVLUTMaxVal":
					CheckBox CBManMax,win=$sva.win,value=1
					CheckBox CBAutoMax,win=$sva.win,value=0
					break
				default:
					Abort "Unknown control in SVLUTProc()"
			EndSwitch
			
			UpdateLUTSetttings(baseWindowName, 0)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function UpdateLUTSetttings(baseWindowName, synchronizeToCurrentSettings)
	string baseWindowName
	variable synchronizeToCurrentSettings	// if true, update the controls to match the current settings. Otherwise set the colortable to match the controls.
	
	if (synchronizeToCurrentSettings == 0)
		ControlInfo /W=$baseWindowName#LUTControls PMLUTColorTable
		string colorTable = S_value
		ControlInfo /W=$baseWindowName#LUTControls CBReverseColorTable
		variable reverseCTAB = V_value
		
		ControlInfo /W=$baseWindowName#LUTControls CBAutoMin
		variable autoMin = V_value
		ControlInfo /W=$baseWindowName#LUTControls CBAutoMax
		variable autoMax = V_value
		
		ModifyImage /W=$baseWindowName#CCDViewer M_CCDImage, ctab={, , $colorTable, reverseCTAB}
		
		if (!autoMin)
			ControlInfo /W=$baseWindowName#LUTControls SVLUTMinVal
			variable minLut = V_value
			ModifyImage /W=$baseWindowName#CCDViewer M_CCDImage, ctab={minLut, , , reverseCTAB}
		else
			ModifyImage /W=$baseWindowName#CCDViewer M_CCDImage, ctab={*, , , reverseCTAB}
		endif
		if (!autoMax)
			ControlInfo /W=$baseWindowName#LUTControls SVLUTMaxVal
			variable maxLut = V_value
			ModifyImage /W=$baseWindowName#CCDViewer M_CCDImage, ctab={, maxLUT, , reverseCTAB}
		else
			ModifyImage /W=$baseWindowName#CCDViewer M_CCDImage, ctab={, *, , reverseCTAB}
		endif
	else
		// fetch the current table settings
		string recreationMacro = WinRecreation(baseWindowName, 0)
		// find the line describing the color macro
		string ctabMacro="", thisLine
		variable i
		for (i = 0; ; i+=1)
			thisLine = StringFromList(i, recreationMacro, "\r")
			if (strlen(thisLine) == 0)
				break
			endif
			if (StringMatch(thisLine, "*ModifyImage M_CCDImage ctab=*") == 1)
				ctabMacro = thisLine
				break
			endif
		endfor
		
		if (strlen(ctabMacro) == 0)
			return 0	// no info on color scale available, so nothing we can do
		endif
		
		// parse the ctab descriptor
		variable openBracket = StrSearch(ctabMacro, "{", 0)
		variable closeBracket = StrSearch(ctabMacro, "}", 0)
		ctabMacro = ctabMacro[openBracket + 1, closeBracket - 1]
		
		string lowerLimitStr = StringFromList(0, ctabMacro, ",")
		string upperLimitStr = StringFromList(1, ctabMacro, ",")
		string colorTableName = StringFromList(3, ctabMacro, ",")
		string doReverseStr = StringFromList(4, ctabMacro, ",")
		
		if (GrepString(lowerLimitStr, "\\*"))	// auto lower limit
			CheckBox CBAutoMin, win=$baseWindowName#LUTControls,value=1
			CheckBox CBManMin, win=$baseWindowName#LUTControls,value=0
		else
			variable lowerLimit = str2num(lowerLimitStr)
			CheckBox CBAutoMin, win=$baseWindowName#LUTControls,value=0
			CheckBox CBManMin, win=$baseWindowName#LUTControls,value=1
			SetVariable SVLUTMinVal, win=$baseWindowName#LUTControls, value=_NUM:lowerLimit
		endif
		
		if (GrepString(upperLimitStr, "\\*"))	// auto upper limit
			CheckBox CBAutoMax, win=$baseWindowName#LUTControls,value=1
			CheckBox CBManMax, win=$baseWindowName#LUTControls,value=0
		else
			variable upperLimit = str2num(upperLimitStr)
			CheckBox CBAutoMax, win=$baseWindowName#LUTControls,value=0
			CheckBox CBManMax, win=$baseWindowName#LUTControls,value=1
			SetVariable SVLUTMaxVal, win=$baseWindowName#LUTControls, value=_NUM:upperLimit
		endif
		
		PopupMenu PMLUTColorTable,win=$baseWindowName#LUTControls,popMatch=colorTableName
		
		if (StringMatch(doReverseStr, "0"))
			CheckBox CBReverseColorTable,win=$baseWindowName#LUTControls,value=1
		else
			CheckBox CBReverseColorTable,win=$baseWindowName#LUTControls,value=0
		endif
	endif
	
End

Function CBShowThresholdProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			string baseWindowName = GetBaseWindowName(cba.win)
			GenerateGUIThreshold(baseWindowName)
			break
	endswitch

	return 0
End

Function CBCheckShowParticlesProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			string BaseWindowName = GetBaseWindowName(cba.win)
			GenerateGUIThreshold(BaseWindowName)
			break
	endswitch

	return 0
End

Function SVUpdateThreshold(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string baseWindowName = GetBaseWindowName(sva.win)
			
			GenerateGUIthreshold(baseWindowName)
			
			break
	endswitch

	return 0
End

Function PMParticleVerificationProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			string baseWindowName = GetBaseWindowName(pa.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			
			NVAR particleVerifierOverlap = windowDataFolder:V_particleVerifierOverlap
			NVAR particleVerifierSymm = windowDataFolder:V_particleVerifierSymm
			NVAR particleVerifierEllipse = windowDataFolder:V_particleVerifierEllipse
			NVAR particleVerifierEllipseAstig = windowDataFolder:V_particleVerifierEllipseAstig
			
			StrSwitch (popStr)
				case "Remove Overlapping Particles":
					particleVerifierOverlap = !particleVerifierOverlap
					break
				case "Gaussian Fitting":
					particleVerifierSymm = !particleVerifierSymm
					break
				case "Ellipsoidal Gaussian Fitting":
					particleVerifierEllipse = !particleVerifierEllipse
					break
				case "Ellipsoidal Gaussian Fitting (astigmatism)":
					particleVerifierEllipseAstig = !particleVerifierEllipseAstig
					break
				default:
					Abort "Unknown particle verifier"
			EndSwitch
			
			GenerateGUIThreshold(baseWindowName)
			
			break
	endswitch

	return 0
End

Function /S PMParticleVerificationListItems()
	// this function needs to know which particle verifiers are enabled
	// the problem is that that depends on the particular window we're looking at
	// but this function doesn't have a clean way of finding that out
	// so we will use a bit of a hack, looking at what the top window is
	string windowName = GetNameOfTopInterfacePanel()
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	
	NVAR particleVerifierOverlap = windowDataFolder:V_particleVerifierOverlap
	NVAR particleVerifierSymm = windowDataFolder:V_particleVerifierSymm
	NVAR particleVerifierEllipse = windowDataFolder:V_particleVerifierEllipse
	NVAR particleVerifierEllipseAstig = windowDataFolder:V_particleVerifierEllipseAstig
	
	string listString = ""
	string checkMark = "\\M0:!" + num2char(18)+":"
	
	if (particleVerifierOverlap != 0)
		listString += checkMark + "Remove Overlapping Particles" + ";"
	else
		listString += "Remove Overlapping Particles" + ";"
	endif
	
	if (particleVerifierSymm != 0)
		listString += checkMark + "Gaussian Fitting" + ";"
	else
		listString += "Gaussian Fitting" + ";"
	endif
	
	if (particleVerifierEllipse != 0)
		listString += checkMark + "Ellipsoidal Gaussian Fitting" + ";"
	else
		listString += "Ellipsoidal Gaussian Fitting" + ";"
	endif
	
	if (particleVerifierEllipseAstig != 0)
		listString += checkMark + "Ellipsoidal Gaussian Fitting (astigmatism)" + ";"
	else
		listString += "Ellipsoidal Gaussian Fitting (astigmatism)" + ";"
	endif
	
	return listString
		
End

Function PMUpdateThresholdMethodProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			string baseWindowName = GetBaseWindowName(pa.win)
			
			SetVariable SVSetThreshold, win=$baseWindowName#AnalysisControls, disable=1
			SetVariable SVSetPFA, win=$baseWindowName#AnalysisControls, disable=1
			SetVariable SVSetSmoothSigmaFactor, win=$baseWindowName#AnalysisControls, disable=1
			
			StrSwitch (pa.popStr)
				case "GLRT":
					SetVariable SVSetPFA, win=$baseWindowName#AnalysisControls, disable=0
					break
				case "Absolute":
					SetVariable SVSetThreshold, win=$baseWindowName#AnalysisControls, disable=0
					break
				case "SmoothSigma":
					SetVariable SVSetSmoothSigmaFactor, win=$baseWindowName#AnalysisControls, disable=0
					break
					break
				default:
					Abort "Unknown segmentation option in PMUpdateThresholdMethodProc"
					break
			EndSwitch
			
			GenerateGUIThreshold(baseWindowName)
			break
	endswitch
	return 0
End

Function PMUpdateViewerProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			string baseWindowName = GetBaseWindowName(pa.win)
			GenerateGUIThreshold(baseWindowName)
			break
	endswitch
	return 0
End

Function CBDisplayFittedPositionsProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			string baseWindowName = GetBaseWindowName(cba.win)
			GenerateGUIthreshold(baseWindowName)
				
			break
	endswitch

	return 0
End

Function RBParticleTrackDisplayProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string controlName = cba.ctrlName
			string baseWindowName = GetBaseWindowName(cba.win)
			
			CheckBox RBDisplayFullTracks, win=$baseWindowName#TrackingControls,value=0
			CheckBox RBDisplayUpdatingTracks, win=$baseWindowName#TrackingControls,value=0
			CheckBox RBDisplayAllKnownTracks, win=$baseWindowName#TrackingControls,value=0
			
			StrSwitch (controlName)
				case "RBDisplayFullTracks":
					CheckBox RBDisplayFullTracks, win=$baseWindowName#TrackingControls,value=1
					break
				case "RBDisplayUpdatingTracks":
					CheckBox RBDisplayUpdatingTracks, win=$baseWindowName#TrackingControls,value=1
					break
				case "RBDisplayAllKnownTracks":
					CheckBox RBDisplayAllKnownTracks, win=$baseWindowName#TrackingControls,value=1
					break
				default:
					Abort "Unknown control name in RBParticleTrackDisplayProc"
					break
			EndSwitch
			
			GenerateGUIThreshold(baseWindowName)
			break
	endswitch

	return 0
End

Function RBSOFICorrelationModeProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			string controlName = cba.ctrlName
			string baseWindowName = GetBaseWindowName(cba.win)
			
			CheckBox RBAutoSOFI, win=$baseWindowName#SOFIControls, value=0
			CheckBox RBCrossSOFI, win=$baseWindowName#SOFIControls, value=0
			
			StrSwitch (controlName)
				case "RBAutoSOFI":
					CheckBox RBAutoSOFI, win=$baseWindowName#SOFIControls, value=1
					break
				case "RBCrossSOFI":
					CheckBox RBCrossSOFI, win=$baseWindowName#SOFIControls, value=1
					break
				default:
					Abort "Unknown control name in RBSOFICorrelationModeProc"
					break
			EndSwitch
			break
	endswitch

	return 0
End

Function SVSetImage(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			string baseWindowName = GetBaseWindowName(sva.win)
			
			ChangeCurrentImage(dval, baseWindowName)
			
			break
	endswitch

	return 0
End

Function BTAnalyzePALMProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	variable msTimer, sec, err

	switch( ba.eventCode )
		case 2: // mouse up
				// click code here
				Struct FitData fitParams
				string baseWindowName = GetBaseWindowName(ba.win)
				
				// get all the fit options and parameters the user defined in the GUI
				GetGUIFitSettingsAndOptions(fitParams, baseWindowName)
				
				msTimer = startMSTimer
				
				string posName
				err = DoPALMFitting(fitParams, baseWindowName, posName, 1)
				if (err != 0)	// an error occured
					Abort "An error occurred during the analysis"
				endif
				
				sec = stopMSTimer(msTimer)
				sec /= 1e6
				Printf "Elapsed time is %s\r", PrettyPrintTimeElapsed(sec)
				
				ControlUpdate /W=$baseWindowName#AnalysisControls PMSelectPositionsWave
				ControlUpdate /W=$baseWindowName#TrackingControls PMSelectPosWaveForTracking
				ControlUpdate /W=$baseWindowName#OtherControls PMDriftCorrSelectPositionsWave
				PopupMenu PMSelectPositionsWave win=$baseWindowName#AnalysisControls, mode=1, popMatch=posName
				
				// update the displayed waves in case we were visualizing positions that have now been overwritten
				GenerateGUIThreshold(baseWindowName)
				break

	endswitch
	return 0
End

Function BTCopyPositionsProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			string baseWindowName = GetBaseWindowName(ba.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWIndowName)
			NewDataFolder /O root:'Localized Positions'
			DFREF outputFolder = root:'Localized Positions'
			string posName
			
			SVAR filePath = windowDataFolder:S_filePath
			
			ControlInfo /W=$baseWindowName#AnalysisControls PMSelectPositionsWave
			posName = S_Value
			
			if (StringMatch(posName, "* No Positions *") == 1)
				return 0
			endif
			
			wave /Z posWave = GetPositionsWaveReference(posName)
			
			if (WaveExists(posWave) == 0)
				SetDataFolder savDF
				Abort "The specified positions wave doesn't exist"
			endif
			
			// propose a name based on the filename
			string fileName = ParseFilePath(3, filePath, ":", 0, 1)
			string proposedName = CleanupName(fileName, 0)
			
			posName = GetNewPositionsWaveName(proposedName, outputFolder, "Enter the new name of the positions wave")
			if (strlen(posName) == 0)	// the user canceled the dialog
				return -1
			endif
			
			Duplicate /O posWave, outputFolder:$posName
			break
	endswitch

	return 0
End

Function BTRunParticleTrackingProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	variable msTimer, sec, err

	Switch( ba.eventCode )
		case 2: // mouse up
				// click code here
			string baseWindowName = GetBaseWindowName(ba.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			DFREF trackDataFolder = root:'Particle Tracks'
			
			NVAR trackingMaxJump = windowDataFolder:V_trackingMaxJump
			NVAR trackingMaxBlinking = windowDataFolder:V_trackingMaxBlinking
			NVAR trackingMaxLength =  windowDataFolder:V_trackingMaxLength
			NVAR trackingMinLength =  windowDataFolder:V_trackingMinLength
			NVAR trackingTimeDelta =  windowDataFolder:V_trackingTimeDelta
			
			// get the name of the selected positions wave
			ControlInfo /W=$baseWindowName#TrackingControls PMSelectPosWaveForTracking
			string posName = S_Value
			
			if (StringMatch(posName, "* No Positions *") == 1)
				return 0
			endif
			
			wave /Z pos = GetPositionsWaveReference(posName)
			if (WaveExists(pos) != 1)
				Abort "The selected positions wave was not found"
			endif
			
			string suggestedName = UniqueName(CleanupName(posName + "_track", 0), 1, 0)
			string outputWaveName = GetNewPositionsWaveName(suggestedName, trackDataFolder, "")
			if (strlen(outputWaveName) == 0)
				return 0
			endif
			
			wave /WAVE W_GroupedEmitters = GroupEmitters(pos, CombinedPosEstimator_LastPos, trackingMaxJump, trackingMaxBlinking)
			
			variable nGroups = DimSize(W_GroupedEmitters, 0)
			variable nGroupsRemaining = nGroups
			variable i
			
			// eliminate all entries that do not fit between the requested lengths
			for (i = 0; i < nGroups; i+=1)
				wave currentEmitter = W_GroupedEmitters[i]
				if ((DimSize(currentEmitter, 0) < trackingMinLength) || (DimSize(currentEmitter, 0) > trackingMaxLength))
					nGroupsRemaining -= 1
				endif
			endfor
			
			Make /O/N=(nGroupsRemaining) /WAVE root:'Particle Tracks':$outputWaveName
			wave /WAVE outputWave = root:'Particle Tracks':$outputWaveName
			
			variable offset = 0
			for (i = 0; i < nGroups; i+=1)
				wave currentEmitter = W_GroupedEmitters[i]
				if ((DimSize(currentEmitter, 0) < trackingMinLength) || (DimSize(currentEmitter, 0) > trackingMaxLength))
					continue
				endif
				outputWave[offset] = currentEmitter
				offset += 1
			endfor
			
			// add a wavenote to the output wave
			string waveNote = Note(pos) + "TRACKING MAX SHIFT:" + num2str(trackingMaxJump) + ";" + "TRACKING MAX BLINKING:" + num2str(trackingMaxBlinking) + ";"
			waveNote += "TRACKING MIN LENGTH:" + num2str(trackingMinLength) + ";" + "TRACKING MAX LENGTH:" + num2str(trackingMaxLength) + ";"
			if ((NumType(trackingTimeDelta) == 0) && (trackingTimeDelta > 0))
				waveNote += "TRACKING TIME DELTA:" + num2str(trackingTimeDelta) + ";"
			endif
			Note /K outputWave, waveNote
			
			// update the popup menu showing the tracks
			ControlUpdate /W=$baseWindowName#TrackingControls PMSelectTrackWave
			PopupMenu PMSelectTrackWave, win=$baseWindowName#TrackingControls, popMatch=outputWaveName
			
			// update the waves are displayed in case we are replacing a tracking wave that is currently displayed
			GenerateGUIThreshold(baseWindowName)
			break
	EndSwitch
End

Function PMSOFIOrderProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			// time lags are not supported with the NewSOFI operation
//			string baseWindowName = GetBaseWindowName(pa.win)
//			string controlName
//			variable orderToUse = str2num(popStr)
//			variable i
//			for (i = 0; i < 10; i+=1)
//				controlName = "SVSOFILagTime" + num2istr(i)
//				if (i < orderToUse - 1)
//					SetVariable /Z $controlName win=$baseWindowName#SOFIControls,disable=0
//				else
//					SetVariable /Z $controlName win=$baseWindowName#SOFIControls,disable=1
//				endif
//			endfor
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function GetSOFIQualityFactor(baseWindowName)
	String baseWindowName
	
	ControlInfo /W=$baseWindowName#SOFIControls PMCalculationQuality
	StrSwitch (S_value)
		case "Minimum":
			return 1
		case "Few":
			return 10
		case "More":
			return 50
		case "All":
			return 1e50
		default:
			Abort "Unknown SOFI quality factor"
	EndSwitch
End

Function BTCalcCombinationWeights(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	Switch( ba.eventCode )
		case 2: // mouse up
				// click code here
			string baseWindowName = GetBaseWindowName(ba.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			NewDataFolder /O root:'SOFI Combination Weights'
			DFREF combinationWeightsFolder = root:'SOFI Combination Weights'
			
			SVAR dataFilePath = windowDataFolder:S_filePath
			NVAR framesToSkip = windowDataFolder:V_SOFInFramesToSkip
			NVAR lastFrameInAnalysis = windowDataFolder:V_SOFILastFrameInAnalysis
			variable fileType
			NVAR /Z gFileType = windowDataFolder:V_fileType
			if (!NVAR_Exists(gFileType))
				fileType = -1
			else
				fileType = gFileType
			endif
			if (lastFrameInAnalysis <= framesToSkip)
				Abort "Invalid frame range"
			endif
			variable nImagesToInclude = lastFrameInAnalysis - framesToSkip + 1
			
			variable sofiOrder
			ControlInfo /W=$baseWindowName#SOFIControls PMSOFIOrder
			sofiOrder = str2num(S_value)
			variable nPoints
			switch (sofiOrder)
				case 2:
					nPoints = 250
					break
				case 3:
					nPoints = 1000
					break
				case 4:
					nPoints = 2000
					break
				default:
					nPoints = 5000
					break
			endswitch
			
			string fileName = ParseFilePath(3, dataFilePath, ":", 0, 1)
			string outputWaveName = GetNewPixelCombinationsWaveName(CleanupName("Wght_" + num2str(sofiOrder) + "_" + fileName, 0), combinationWeightsFolder, "Output wave:")
			if (strlen(outputWaveName) == 0)
				return 0
			endif
			
			ReadCCDImages /Y=(fileType)/O/Q/S=(framesToSkip) /C=(nImagesToInclude) /DEST=M_CCDFrames_Weights dataFilePath
			wave M_CCDFrames_Weights
			variable nImages = DimSize(M_CCDFrames_Weights, 2)
			MatrixOP /FREE M_Average_Weights = sumbeams(M_CCDFrames_Weights) / nImages 
			wave points = SampleProbability(M_Average_Weights, nPoints)
			variable comb = inf
			wave /Z W_CombinationWeights = PrepareWeights(comb, sofiOrder, points, M_CCDFrames_Weights)
			KillWaves /Z M_CCDFrames_Weights, points, M_Average_Weights
			if (!WaveExists(W_CombinationWeights))
				return 0
			endif
			
			Duplicate /O W_CombinationWeights, combinationWeightsFolder:$outputWaveName
			wave W_Combs = combinationWeightsFolder:$outputWaveName
			Note W_Combs, "TYPE:SOFICombinationWeights;"
			KillWaves /Z W_CombinationWeights
			break
		EndSwitch
End

Function BTDoSOFIProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	Switch( ba.eventCode )
		case 2: // mouse up
				// click code here
			string baseWindowName = GetBaseWindowName(ba.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			
			SVAR dataFilePath = windowDataFolder:S_filePath
			NVAR framesToSkip = windowDataFolder:V_SOFInFramesToSkip
			NVAR lastFrameInAnalysis = windowDataFolder:V_SOFILastFrameInAnalysis
			NVAR SOFLagTime0 = windowDataFolder:V_SOFLagTime0
			NVAR SOFLagTime1 = windowDataFolder:V_SOFLagTime1
			variable fileType
			NVAR /Z gFileType = windowDataFolder:V_fileType
			if (!NVAR_Exists(gFileType))
				fileType = -1
			else
				fileType = gFileType
			endif
			
			if (lastFrameInAnalysis <= framesToSkip)
				Abort "Invalid frame range"
			endif
			variable nImagesToInclude = lastFrameInAnalysis - framesToSkip + 1
			
			variable sofiOrder
			ControlInfo /W=$baseWindowName#SOFIControls PMSOFIOrder
			sofiOrder = str2num(S_value)
			
			ControlInfo /W=$baseWindowName#SOFIControls RBCrossSOFI
			variable doCrossCorrelation = V_value
			ControlInfo /W=$baseWindowName#SOFIControls CBSOFIDoAverageImage
			variable doAverage = V_value
			ControlInfo /W=$baseWindowName#SOFIControls CBSOFIMakeMovie
			variable makeMovie = V_value
			ControlInfo /W=$baseWindowName#SOFIControls CBSOFIDropSaturatedFrames
			variable dropSaturation = V_value
			variable nFramesToGroup
			if (makeMovie != 0)
				ControlInfo /W=$baseWindowName#SOFIControls SVSOFInFramesInMovie
				nFramesToGroup = V_value
			else
				nFramesToGroup = 0
			endif
			variable sofiQualityFactor = GetSOFIQualityFactor(baseWindowName)
			ControlInfo /W=$baseWindowName#SOFIControls CBUseCombinationWeights
			variable useWeights = V_Value
			if (useWeights)
				ControlInfo /W=$baseWindowName#SOFIControls PMCombinationWeights
				wave /Z W_CombinationWeights = root:'SOFI Combination Weights':$S_value
				if (!WaveExists(W_CombinationWeights))
					Abort "no weights wave selected"
				endif
				sofiQualityFactor = 1e50
			endif
			
			// lag times
			Make /N=2 /D /FREE W_LagTimes = {SOFLagTime0, SOFLagTime1}
			
			// add a structure with a reference to the progress reporting function
			STRUCT LocalizerProgStruct progressStruct
			progressStruct.version = kLocalizerProgStructVersion
			FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
			
			//SOFIAnalysis /Y=(fileType)/XC=(doCrossCorrelation) /ORDR=(sofiOrder) /GRP=(nFramesToGroup) /SUB={framesToSkip, nImagesToInclude} /AVG=(doAverage) /NSAT=(dropSaturation) /LAGW=W_LagTimes /PROG=progressStruct /DEST=windowDataFolder:M_SOFI dataFilePath
			variable i
			if (nFramesToGroup == 0)
				if (!useWeights)
					NewSOFI /Y=(fileType) /ORDR=(sofiOrder) /COMB=(sofiQualityFactor) /SUB={framesToSkip, nImagesToInclude} /AVG=(doAverage) /NSAT=(dropSaturation) /PROG=progressStruct /DEST=windowDataFolder:M_SOFI dataFilePath
				else
					NewSOFI /Y=(fileType) /ORDR=(sofiOrder) /COMB=(sofiQualityFactor) /SUB={framesToSkip, nImagesToInclude} /WGHT=W_CombinationWeights /AVG=(doAverage) /NSAT=(dropSaturation) /PROG=progressStruct /DEST=windowDataFolder:M_SOFI dataFilePath
				endif
			else
				ReadCCDImages /H/Y=(fileType) dataFilePath
				variable nImages = V_numberOfImages
				variable localNFramesToSkip = max(localNFramesToSkip, 0)
				nImagesToInclude = min(nImagesToInclude, nImages - localNFramesToSkip)
				variable nGroups = floor(nImagesToInclude / nFramesToGroup)
				if (nGroups <= 1)
					Abort "not enough images to make the requested movie"
				endif
				for (i = 0; i < nGroups; i+=1)
					if (!useWeights)
						NewSOFI /Y=(fileType) /ORDR=(sofiOrder) /COMB=(sofiQualityFactor) /SUB={framesToSkip + i * nFramesToGroup, nFramesToGroup} /AVG=(doAverage) /NSAT=(dropSaturation) /PROG=progressStruct /DEST=windowDataFolder:M_MovieSOFI dataFilePath
					else
						NewSOFI /Y=(fileType) /ORDR=(sofiOrder) /COMB=(sofiQualityFactor) /SUB={framesToSkip + i * nFramesToGroup, nFramesToGroup} /WGHT=W_CombinationWeights /AVG=(doAverage) /NSAT=(dropSaturation) /PROG=progressStruct /DEST=windowDataFolder:M_MovieSOFI dataFilePath
					endif
					wave M_MovieSOFI = windowDataFolder:M_MovieSOFI
					wave /Z M_MovieSOFI_avg = windowDataFolder:M_MovieSOFI_Avg
					if (i == 0)
						Duplicate /O M_MovieSOFI, windowDataFolder:M_SOFI
						Redimension /N=(-1,-1,nGroups) windowDataFolder:M_SOFI
						if (WaveExists(M_MovieSOFI_avg))
							Duplicate /O M_MovieSOFI_avg, windowDataFolder:M_SOFI_Avg
							Redimension /N=(-1,-1,nGroups) windowDataFolder:M_SOFI_Avg
						endif
					else
						ImageTransform /D=M_MovieSOFI /P=(i) setPlane, windowDataFolder:M_SOFI
						if (WaveExists(M_MovieSOFI_avg))
							ImageTransform /D=M_MovieSOFI_avg /P=(i) setPlane, windowDataFolder:M_SOFI_Avg
						endif
					endif
				endfor
				KillWaves /Z M_MovieSOFI, M_MovieSOFI_Avg
			endif
			wave M_SOFI = windowDataFolder:M_SOFI
			
			// provide a wave note that encapsulates the relevant information
			string waveNote = "KIND:SOFI;ORDER:" + num2istr(sofiOrder) + ";CROSSCORRELATION:" + num2istr(doCrossCorrelation) + ";GROUPING:" + num2istr(nFramesToGroup) + ";SKIP:" + num2istr(framesToSkip) + ";DROPSATURATION:" + num2istr(dropSaturation) + ";ORIGINAL FILE PATH:" + dataFilePath + ";"
			Note /K M_SOFI, waveNote
			
			string SOFIWindowTitle = "SOFI image from " + ParseFilePath(3, dataFilePath, ":", 0, 1)
			string SOFIWindowName = ReplaceString("LocalizerViewer", baseWindowName, "SOFIViewer")
			
			DoWindow /F $SOFIWindowName
			if (V_flag == 0)
				Display /K=1 /N=$SOFIWindowName as SOFIWindowTitle
				AppendImage /G=1 /W=$SOFIWindowName M_SOFI
				AdjustGraphForImageDisplay(SOFIWindowName)
				
				ModifyGraph /W=$SOFIWindowName width={Plan,1,bottom,left}
				
				// if this is a movie (more than one image)
				// then show a slider automatically
				if (DimSize(M_SOFI, 2) > 1)
					WMAppend3DImageSlider()
				else
					// this is really a 2D image, and some procedures (e.g. ImageSave)
					// choke on mxmx1 matrices. So turn it into one explictly
					Redimension /N=(-1, -1) M_SOFI
				endif
			endif
			
			if (doAverage)
				string averageWindowTitle = "Average/SOFI image from " + ParseFilePath(3, dataFilePath, ":", 0, 1)
				string averageWindowName = ReplaceString("LocalizerViewer", baseWindowName, "AverageSOFIViewer")
			
				DoWindow /F $averageWindowName
				if (V_flag == 0)
					wave M_SOFI_avg = windowDataFolder:M_SOFI_avg
					
					Display /K=1 /N=$averageWindowName as averageWindowTitle
					AppendImage /G=1 /W=$averageWindowName M_SOFI_avg
					AdjustGraphForImageDisplay(averageWindowName)
					
					ModifyGraph /W=$averageWindowName width={Plan,1,bottom,left}
					
					// if this is a movie (more than one image)
					// then show a slider automatically
					if (DimSize(M_SOFI_avg, 2) > 1)
						WMAppend3DImageSlider()
					else
						// this is really a 2D image, and some procedures (e.g. ImageSave)
						// choke on mxmx1 matrices. So turn it into one explictly
						Redimension /N=(-1, -1) M_SOFI_avg
					endif
					
					AutoPositionWindow /M=0 /R=$SOFIWindowName $averageWindowName
				endif
			endif
			break
	EndSwitch
End

Function BTDoSOFIDeconvolutionProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	Switch( ba.eventCode )
		case 2:	// mouse up
			string baseWindowName = GetBaseWindowName(ba.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			
			NVAR gaussianWidth = windowDataFolder:V_gaussianWidth
			NVAR SOFInDeconvolutionIterations = windowDataFolder:V_SOFInDeconvolutionIterations
			
			SVAR dataFilePath = windowDataFolder:S_filePath
			
			wave /Z M_SOFI = windowDataFolder:M_SOFI
			if (WaveExists(M_SOFI) != 1)
				Abort "No SOFI image was found"
			endif
			
			// get the calculation settings from the wavenote
			string waveNote = Note(M_SOFI)
			variable order = NumberByKey("ORDER", waveNote)
			variable isCross = NumberByKey("CROSSCORRELATION", waveNote)
			
			variable effectivePSFWidth
			if (isCross != 0)
				effectivePSFWidth = gaussianWidth
			else
				// generate more pixels in XCSOFI, so the effective psfWidth is larger
				effectivePSFWidth = gaussianWidth * order
			endif
			
			// make an image of the PSF
			variable psfDimensions = ceil(2 * 4 * effectivePSFWidth)
			if (mod(psfDimensions, 2) != 1)
				psfDimensions += 1
			endif
			
			Make /D/FREE/N=(psfDimensions, psfDimensions) M_PSF
			SetScale /P x, - (psfDimensions - 1) / 2, 1, M_PSF
			SetScale /P y, - (psfDimensions - 1) / 2, 1, M_PSF
			
			M_PSF = StatsNormalPDF(sqrt(x^2 + y^2), 0, effectivePSFWidth)^order
			
			// if the SOFI image really is a movie then multiple frames need to be processed
			variable nFrames = DimSize(M_SOFI, 2)
			if (nFrames == 0)
				nFrames = 1
			endif
			
			Make /O/N=(DimSize(M_SOFI, 0), DimSize(M_SOFI, 1), nFrames) /D windowDataFolder:M_DeconvolvedSOFI
			wave M_DeconvolvedSOFI = windowDataFolder:M_DeconvolvedSOFI
			
			variable i
			for (i = 0; i < nFrames; i+=1)
				MatrixOP /O /FREE M_ExtractedPlane = M_SOFI[][][i]
				WaveStats /Q/M=1 M_ExtractedPlane
				if (V_min < 0)	// contains negative values, these appear to produce noise in the deconvolution
					M_ExtractedPlane = (M_ExtractedPlane[p][q] < 0) ? 0 : M_ExtractedPlane[p][q]
				endif
				ImageRestore /ITER=(SOFInDeconvolutionIterations) /DEST=windowDataFolder:M_DeconvolvedSOFI_Temp srcWave=M_ExtractedPlane, psfWave=M_PSF
				ImageTransform /D=windowDataFolder:M_DeconvolvedSOFI_Temp /P=(i) setPlane, M_DeconvolvedSOFI
			endfor
			
			KillWaves /Z windowDataFolder:M_DeconvolvedSOFI_Temp
			
			// copy the original wave scaling as well
			CopyScales M_SOFI, M_DeconvolvedSOFI
			
			string windowTitle = "Deconvolved SOFI image from " + ParseFilePath(3, dataFilePath, ":", 0, 1)
			string windowName = ReplaceString("LocalizerViewer", baseWindowName, "DeconvolvedSOFIViewer")
			
			DoWindow /F $windowName
			if (V_flag == 0)
				Display /K=1 /N=$windowName as windowTitle
				AppendImage /G=1 /W=$windowName M_DeconvolvedSOFI
				AdjustGraphForImageDisplay(windowName)
				ModifyGraph /W=$windowName width={Plan,1,bottom,left}
				
				// if this is a movie (more than one image)
				// then show a slider automatically
				if (DimSize(M_DeconvolvedSOFI, 2) > 1)
					WMAppend3DImageSlider()
				else
					// this is really a 2D image, and some procedures (e.g. ImageSave)
					// choke on mxmx1 matrices. So turn it into one explictly
					Redimension /N=(-1, -1) M_DeconvolvedSOFI
				endif
			endif
			break
	EndSwitch
End

Function BTBasicAnalysisProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string baseWindowName = GetBaseWindowName(ba.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			
			variable analysisMethod
			SVAR CCDFilePath = windowDataFolder:S_filePath
			variable fileType
			NVAR /Z gFileType = windowDataFolder:V_fileType
			if (!NVAR_Exists(gFileType))
				fileType = -1
			else
				fileType = gFileType
			endif
			
			string fileName = ParseFilePath(5, CCDFilePath, ":", 0, 0)
			fileName = ParseFilePath(3, fileName, ":", 0, 0)
			
			ControlInfo /W=$baseWindowName#OtherControls PMChooseCCDAnalysisMethod
			StrSwitch (S_Value)
				case "Average Intensity Trace":
					analysisMethod = ANALYZECCD_AVERAGETRACE
					break
				case "Summed Intensity Trace":
					analysisMethod = ANALYZECCD_SUMMEDINTENSITYTRACE
					break
				case "Average Image":
					analysisMethod = ANALYZECCD_AVERAGE_IMAGE
					break
				case "Variance Image":
					analysisMethod = ANALYZECCD_VARIANCE_IMAGE
					break
				default:
					Abort "Unknown CCD analysis method"
			EndSwitch
			
			// add a structure with a reference to the progress reporting function
			STRUCT LocalizerProgStruct progressStruct
			progressStruct.version = kLocalizerProgStructVersion
			FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
			
			string graphWindowName
			string graphTitle
			GetWindow $baseWindowName, title
			string baseWindowTitle = S_Value
			
			switch (analysisMethod)
				case ANALYZECCD_SUMMEDINTENSITYTRACE:
					AnalyzeCCDImages /Y=(fileType)/M=(analysisMethod) /PROG=progressStruct /DEST=windowDataFolder:W_SummedIntensityTrace CCDFilePath
					wave W_SummedIntensityTrace = windowDataFolder:W_SummedIntensityTrace
					graphWindowName = ReplaceString("LocalizerViewer", baseWindowname, "SummedIntensityTrace")
					graphTitle = "Summed intensity trace from " + baseWindowTitle
					DoWindow /F $graphWindowName
					if (V_flag != 1)
						Display /K=1 /N=$graphWindowName
						AppendToGraph /W=$S_name W_SummedIntensityTrace
						Label /W=$S_name Left, "Summed Intensity"
						Label /W=$S_name Bottom, "Frame Number"
					endif
					DoWindow /T $graphWindowName, graphTitle
					break
				case ANALYZECCD_AVERAGETRACE:
					AnalyzeCCDImages /Y=(fileType)/M=(analysisMethod) /PROG=progressStruct /DEST=windowDataFolder:W_AverageIntensityTrace CCDFilePath
					wave W_AverageIntensityTrace = windowDataFolder:W_AverageIntensityTrace
					graphWindowName = ReplaceString("LocalizerViewer", baseWindowname, "AverageIntensityTrace")
					graphTitle = "Average intensity trace from " + baseWindowTitle
					DoWindow /F $graphWindowName
					if (V_flag != 1)
						Display /K=1 /N=$graphWindowName
						AppendToGraph /W=$S_name W_AverageIntensityTrace
						Label /W=$S_name Left, "Average Intensity"
						Label /W=$S_name Bottom, "Frame Number"
					endif
					DoWindow /T $graphWindowName, graphTitle
					break
				case ANALYZECCD_AVERAGE_IMAGE:	// calculate a image containing an average of the entire trace
					AnalyzeCCDImages /Y=(fileType)/M=(analysisMethod) /PROG=progressStruct /DEST=windowDataFolder:M_AverageImage CCDFilePath
					wave M_AverageImage = windowDataFolder:M_AverageImage
					graphWindowName = ReplaceString("LocalizerViewer", baseWindowname, "AverageImage")
					graphTitle = "Average image from " + baseWindowTitle
					DoWindow /F $graphWindowName
					if (V_flag != 1)
						Display /K=1 /N=$graphWindowName
						AppendImage /W=$S_name M_AverageImage
						AdjustGraphForImageDisplay(S_name)
						ModifyGraph /W=$S_name width={Plan,1,bottom,left}
					endif
					DoWindow /T $graphWindowName, graphTitle
					break
				case ANALYZECCD_VARIANCE_IMAGE:	// calculate an image with the variance of each pixel
					AnalyzeCCDImages /Y=(fileType)/M=(analysisMethod) /PROG=progressStruct /DEST=windowDataFolder:M_Variance CCDFilePath
					wave M_Variance = windowDataFolder:M_Variance
					graphWindowName = ReplaceString("LocalizerViewer", baseWindowname, "VarianceImage")
					graphTitle = "Variance image from " + baseWindowTitle
					DoWindow /F $graphWindowName
					if (V_flag != 1)
						Display /K=1 /N=$graphWindowName
						AppendImage /W=$S_name M_Variance
						AdjustGraphForImageDisplay(S_name)
						ModifyGraph /W=$S_name width={Plan,1,bottom,left}
					endif
					break
				default:
					Abort "Unknown CCD analysis method"
					break
			endswitch
			
			break
	endswitch
	
	return 0
End

Function BTBleachingAnalysisProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string baseWindowName = GetBaseWindowName(ba.win)
			NVAR bleachingMaxJump = windowDataFolder:V_bleachingMaxJump
			
			RecoverEmittersByBleachingSteps(bleachingMaxJump, baseWindowName)
			break
	EndSwitch
End

Function PMChooseMovieProcessingProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			string baseWindowName = GetBaseWindowName(pa.win)
			string selection = pa.popStr
			
			StrSwitch(selection)
				case "Subtract Pixelwise Average":
					SetVariable SVFramesAverageSubtraction,win=$baseWindowName#OtherControls,disable=0
					SetVariable SVCameraMultiplication,win=$baseWindowName#OtherControls,disable=1
					SetVariable SVCameraOffset,win=$baseWindowName#OtherControls,disable=1
					break
				case "Convert to Photons":
					SetVariable SVFramesAverageSubtraction,win=$baseWindowName#OtherControls,disable=1
					SetVariable SVCameraMultiplication,win=$baseWindowName#OtherControls,disable=0
					SetVariable SVCameraOffset,win=$baseWindowName#OtherControls,disable=0
					break
				default:
					SetVariable SVFramesAverageSubtraction,win=$baseWindowName#OtherControls,disable=1
					SetVariable SVCameraMultiplication,win=$baseWindowName#OtherControls,disable=1
					SetVariable SVCameraOffset,win=$baseWindowName#OtherControls,disable=1
					break
			EndSwitch
			break
	endswitch
	return 0
End

Function BTProcessMovieProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string baseWindowName = GetBaseWindowName(ba.win)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			NewDataFolder /O root:'Localizer Movies'
			DFREF movieFolder = root:'Localizer Movies'
			
			NVAR nFramesAverageSubtraction = windowDataFolder:V_nFramesAverageSubtraction
			NVAR cameraMultiplication = windowDataFolder:V_cameraMultiplication
			NVAR cameraOffset = windowDataFolder:V_cameraOffset
			
			SVAR filePath = windowDataFolder:S_filePath
			
			variable fileType
			NVAR /Z gFileType = windowDataFolder:V_fileType
			if (!NVAR_Exists(gFileType))
				fileType = -1
			else
				fileType = gFileType
			endif
			
			// ask how the user wants to save the data  and get the output path
			variable outputFormat
			string outputPath
			GetOutputStorageTypeAndPath(filePath, outputFormat, outputPath)
			if (strlen(outputPath) == 0)
				return 0
			endif
			
			// determine the processing to be done
			variable processMethod
			ControlInfo /W=$baseWindowName#OtherControls PMChooseMovieProcessingMethod
			StrSwitch (S_Value)
				case "Subtract Pixelwise Average":
					processMethod = PROCESSCCDIMAGE_SUBAVG
					break
				case "Difference Image":
					processMethod = PROCESSCCDIMAGE_DIFFIMAGE
					break
				case "Convert File Format":
					processMethod = PROCESSCCDIMAGE_CONVERTFORMAT
					break
				case "Convert to photons":
					processMethod = PROCESSCCDIMAGE_CONVERTPHOTONS
					break
				default:
					Abort "Unknown CCD processing method in StrSwitch in BTProcessMovieProc"
			EndSwitch
			
			// add a structure with a reference to the progress reporting function
			STRUCT LocalizerProgStruct progressStruct
			progressStruct.version = kLocalizerProgStructVersion
			FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
			
			// do the actual work
			ProcessCCDImages /Y=(fileType)/M=(processMethod) /CAL={cameraOffset, cameraMultiplication} /AVG=(nFramesAverageSubtraction) /OUT=(outputFormat) /PROG=progressStruct /O filePath, outputPath
			
			// open the newly created movie
			NewInterfacePanel(outputPath)
			break
	EndSwitch
End

Function SavePositionsToTextFile()
	
	// get the positions we want to save
	String posName
	Prompt posName, "Positions wave: ", popup, GetPossiblePositionsWaves()
	DoPrompt "Select the positions to be saved", posName
	if (V_flag == 1)
		return 0
	endif
	
	wave /Z pos = GetPositionsWaveReference(posName)
	if (WaveExists(pos) != 1)
		Abort "Internal error: cannot find the requested positions"
	endif
	
	// get the output filepath
	variable refNum
	string outputFilePath
	Open /D /F="Text File (*.txt):.txt;All Files:.*;" refNum
	if (strlen(S_fileName) == 0)
		return 0
	endif
	outputFilePath = S_fileName
	
	// open the file for writing
	Open refNum as outputFilePath
	WritePositionsToFile(pos, refNum)
	Close refNum
End

Function SaveTracksToTextFile()
	// get the tracks we want to save
	String tracksName
	Prompt tracksName, "Tracks wave: ", popup, GetPossibleTrackingWaves()
	DoPrompt "Select the tracks to be saved", tracksName
	if (V_flag == 1)
		return 0
	endif
	
	wave /WAVE/Z W_Tracks = GetTrackingWaveReference(tracksName)
	if (WaveExists(W_Tracks) != 1)
		Abort "Internal error: cannot find the requested positions"
	endif
	
	// get the output filepath
	variable refNum
	string outputFilePath
	Open /D /F="Text File (*.txt):.txt;All Files:.*;" refNum
	if (strlen(S_fileName) == 0)
		return 0
	endif
	outputFilePath = S_fileName
	
	// open the file for writing
	Open refNum as outputFilePath
	WriteTracksToFile(W_Tracks, refNum)
	Close refNum
End

Function MakeAverageTraceFromMarquee()
	
	GetMarquee /Z left, bottom
	if (V_flag == 0)
		return 0
	endif
	
	string baseWindowName = GetBaseWindowName(S_MarqueeWin)
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	
	if (GrepString(baseWindowName, "LocalizerViewer[0-9]*") != 1)
		Abort "This procedure can only be called on viewer windows created with the Localizer package"
	endif
	
	wave image = windowDataFolder:M_CCDImage
	SVAR filePath = windowDataFolder:S_FilePath
	variable fileType
	NVAR /Z gFileType = windowDataFolder:V_fileType
	if (!NVAR_Exists(gFileType))
		fileType = -1
	else
		fileType = gFileType
	endif
	
	variable xSize = DimSize(image, 0)
	variable ySize = DimSize(image, 1)
	Variable bottom, left, right, top
	
	bottom = round(min(V_bottom, V_top))
	top = round(max(V_bottom, V_top))
	left = round(min(V_left, V_right))
	right = round(max(V_left, V_right))
	
	// crop the marquee to the dimensions of the image
	bottom = (bottom < 0) ? 0 : bottom
	left = (left < 0) ? 0 : left
	top = (top >= ySize) ? ySize - 1 : top
	right = (right >= xSize) ? xSize - 1 : right
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct progressStruct
	progressStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
	
	// run the operation that will generate the trace
	AnalyzeCCDImages /Y=(fileType)/M=(ANALYZECCD_AVERAGETRACE) /ROI={left, right, bottom, top} /PROG=progressStruct /DEST=windowDataFolder:W_AverageIntensityTrace filePath
	
	// now display the resulting image
	string graphWindowName = ReplaceString("LocalizerViewer", baseWindowName, "AverageIntensityTrace")
	
	GetWindow $baseWindowName, title
	string graphWindowTitle = "Average trace from " + S_value
	
	DoWindow /F $graphWindowName
	if (V_flag == 0)
		Display /K=1 /N=$graphWindowName windowDataFolder:W_AverageIntensityTrace as "Average intensity trace"
		Label /W=$S_name left, "Average intensity"
		Label /W=$S_name bottom, "Frame number"
	endif
	
	DoWindow /T $graphWindowName, graphWindowTitle
End

Function MakeCroppedStackFromMarquee()
	GetMarquee /Z left, bottom
	if (V_flag == 0)
		return 0
	endif
	
	string baseWindowName = GetBaseWindowName(S_MarqueeWin)
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	
	if (GrepString(baseWindowName, "LocalizerViewer[0-9]*") != 1)
		Abort "This procedure can only be called on viewer windows created with the Localizer package"
	endif
	
	wave image = windowDataFolder:M_CCDImage
	SVAR filePath = windowDataFolder:S_FilePath
	variable fileType
	NVAR /Z gFileType = windowDataFolder:V_fileType
	if (!NVAR_Exists(gFileType))
		fileType = -1
	else
		fileType = gFileType
	endif
	
	variable xSize = DimSize(image, 0)
	variable ySize = DimSize(image, 1)
	Variable bottom, left, right, top
	
	bottom = round(min(V_bottom, V_top))
	top = round(max(V_bottom, V_top))
	left = round(min(V_left, V_right))
	right = round(max(V_left, V_right))
	
	// crop the marquee to the dimensions of the image
	bottom = (bottom < 0) ? 0 : bottom
	left = (left < 0) ? 0 : left
	top = (top >= ySize) ? ySize - 1 : top
	right = (right >= xSize) ? xSize - 1 : right
	
	if ((top == bottom) || (left == right))	// it's not a real box
		return -1
	endif
	
	variable outputFormat
	string outputPath
	GetOutputStorageTypeAndPath(filePath, outputFormat, outputPath)
	if (strlen(outputPath) == 0)
		return 0
	endif
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct progressStruct
	progressStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
	
	// do the actual processing
	ProcessCCDImages /Y=(fileType)/O /ROI={left, right, bottom, top} /M=(PROCESSCCDIMAGE_CROP) /OUT=(outputFormat) /PROG=progressStruct filePath, outputPath
	
	// open the newly created movie
	NewInterfacePanel(outputPath)
End

Function Clustering_Marquee(procedureName)
	string procedureName	// either "LClustering", "PairwiseCorrelation", "Direct", or "MLEGauss"
	
	DFREF packageDF = root:Packages:Localizer
	
	variable xPixelSize, yPixelSize, xStart, yStart, xEnd, yEnd
	string marqueeWin
	
	wave /Z positions = GetCoordinatesFromMarquee("positions", marqueeWin, xPixelSize, yPixelSize, xStart, xEnd, yStart, yEnd)
	if (WaveExists(positions) == 0)
		// signals a user cancel
		return 0
	endif
	
	variable nPos = DimSize(positions, 0), xCol, yCol, zCol
	variable nPositionsWithinLimits = 0, i, offset
	
	// get the indices of the columns containing the x and y coordinates
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1))
		Abort "The positions passed to ExtractSubsetOfPositions_menu() do not appear to contain any (x,y) information"
	endif
	
	
	wave M_SelectedPositions = GetPositionsWithinLimits(positions, xStart, xEnd, yStart, yEnd)
	nPositionsWithinLimits = DimSize(M_SelectedPositions, 0)
	
	// if there are no positions within the box then don't do anything
	if (nPositionsWithinLimits == 0)
		return 0
	endif
	
	StrSwitch (procedureName)
		case "LClustering":
		case "PairwiseCorrelation":
			variable isLClustering = StringMatch(procedureName, "LClustering")
			NVAR /Z nMonteCarloSamples= root:Packages:Localizer:V_lClusteringNMonteCarloSamples
			variable nRandomSamples = 0
			if (NVAR_Exists(nMonteCarloSamples))
				nRandomSamples = nMonteCarloSamples
			endif
			variable calculationRange = 0.45 * sqrt((xEnd - xStart)^2 + (yEnd - yStart)^2)
			variable nBins = ceil(calculationRange / 0.1)
			
			 CalculateAndDisplayClustering(isLClustering, M_SelectedPositions, calculationRange, nBins, inf, nRandomSamples, xStart, xEnd, yStart, yEnd)
			
			break
		case "MLEGauss":
			if (nPositionsWithinLimits < 2)
				Abort "Need more positions to fit a Gaussian"
			endif
			wave W_GaussParams = FitMLEGaussToPositions(M_SelectedPositions)
			Make /N=(2,2) /D/FREE M_Covar
			M_Covar[0][0] = W_GaussParams[2]^2
			M_Covar[0][1] = W_GaussParams[4] * W_GaussParams[2] * W_GaussParams[3]
			M_Covar[1][0] = W_GaussParams[4] * W_GaussParams[2] * W_GaussParams[3]
			M_Covar[1][1] = W_GaussParams[3]^2
			MatrixEigenV /L M_Covar
			wave /C W_eigenvalues, M_L_eigenVectors
			variable stdDev1 = sqrt(real(W_eigenvalues[0]))
			variable stdDev2 = sqrt(real(W_eigenvalues[1]))
			variable theta = imag(r2polar(cmplx(M_L_eigenVectors[0][0], M_L_eigenVectors[1][0]))) * 360 / (2 * pi)
			KillWaves /Z W_eigenvalues, M_L_eigenVectors
			
			// display a contour showing the fitted Gaussian
			Make /O/N=(min(max(stdDev1, stdDev2) * 50, 100), min(max(stdDev1, stdDev2) * 50, 100)) /D packageDF:M_MLEContour
			wave M_MLEContour = packageDF:M_MLEContour
			SetScale /I x, W_GaussParams[0] - max(stdDev1, stdDev2) * 3, W_GaussParams[0] + max(stdDev1, stdDev2) * 3, M_MLEContour
			SetScale /I y, W_GaussParams[1] - max(stdDev1, stdDev2) * 3, W_GaussParams[1] + max(stdDev1, stdDev2) * 3, M_MLEContour
			M_MLEContour = MLEGaussFitFunc(x, y, W_GaussParams)
			variable strIndex = FindListItem("M_MLEContour", ContourNameList("", ";"))
			if (strIndex == -1)
				AppendMatrixContour M_MLEContour
				ModifyContour M_MLEContour labels=0
			endif
			
			Printf "In pixels: %u points, center is at (%g, %g), stdDev1 = %g, stdDev2 = %g, angle is %g degrees\r", nPositionsWithinLimits, W_GaussParams[0], W_GaussParams[1], stdDev1, stdDev2, theta
			Make /O/N=6 /D W_FittedGaussianParams = {nPositionsWithinLimits, W_GaussParams[0], W_GaussParams[1], stdDev1, stdDev2, theta}
			if (xPixelSize != 0)
				Printf "In nanometer: %u points, center is at (%g nm, %g nm), stdDev1 = %g nm, stdDev2 = %g nm, angle is %g degrees\r", nPositionsWithinLimits, W_GaussParams[0] * xPixelSize, W_GaussParams[1] * xPixelSize, stdDev1 * xPixelSize, stdDev2 * xPixelSize, theta
				Make /O/N=6 /D W_FittedGaussianParams_nm = {nPositionsWithinLimits, W_GaussParams[0] * xPixelSize, W_GaussParams[1] * xPixelSize, stdDev1 * xPixelSize, stdDev2 * xPixelSize, theta}
			else
				KillWaves /Z W_FittedGaussianParams_nm
			endif
			break
		case "Direct":
			// get the center of the cluster
			variable clusterCenterX = 0, clusterCenterY = 0
			for (i = 0; i < nPositionsWithinLimits; i+=1)
				clusterCenterX += M_SelectedPositions[i][xCol]
				clusterCenterY += M_SelectedPositions[i][yCol]
			endfor
			clusterCenterX /= nPositionsWithinLimits
			clusterCenterY /= nPositionsWithinLimits
			
			// calculate the distances of all the positions from the center
			Make /D/O/N=(nPositionsWithinLimits) /FREE W_CenterDistances
			for (i = 0; i < nPositionsWithinLimits; i+=1)
				W_CenterDistances[i] = sqrt((M_SelectedPositions[i][xCol] - clusterCenterX)^2 + (M_SelectedPositions[i][yCol] - clusterCenterY)^2)
			endfor
			
			WaveStats /Q W_CenterDistances
			variable stdDev = V_sdev
			variable maxDistance = V_max
			variable avgDistance = V_avg
			Printf "In pixels: %u points, center is at (%g, %g), average distance from the center is %g, standard deviation is %g, radius of enclosing circle is %g\r", nPositionsWithinLimits, clusterCenterX, clusterCenterY, avgDistance, stdDev, maxDistance
			Make /O/N=6 /D W_DirectClusterParams = {nPositionsWithinLimits, clusterCenterX, clusterCenterY, avgDistance, stdDev, maxDistance}
			// also print the information in pixels if available
			if (xPixelSize != 0)
				string formatString = "In nanometer: %u points, center is at (%g nm, %g nm), average distance from the center is %g nm, standard deviation is %g nm, radius of enclosing circle is %g nm\r"
				Printf formatString, nPositionsWithinLimits, clusterCenterX * xPixelSize, clusterCenterY * xPixelSize, avgDistance * xPixelSize, stdDev * xPixelSize, maxDistance * xPixelSize
				Make /O/N=6 /D W_DirectClusterParams_nm = {nPositionsWithinLimits, clusterCenterX * xPixelSize, clusterCenterY * xPixelSize, avgDistance * xPixelSize, stdDev * xPixelSize, maxDistance * xPixelSize}
			else
				KillWaves /Z W_DirectClusterParams_nm
			endif
			break
		default:
			Abort "Unknown analysis method in Clustering_Marquee"
			break
	EndSwitch
End

Function CalculateMSDsCDFs_Marquee(msdOrCDFStr)
	string msdOrCDFStr	// "msd" or "cdf"
	variable xPixelSize, yPixelSize, xStart, yStart, xEnd, yEnd
	string marqueeWin
	
	variable useMSD = 0
	StrSwitch (msdOrCDFStr)
		case "msd":
			useMSD = 1
			break
		case "cdf":
			useMSD = 0
			break
		default:
			Abort "Unknown type passed to CalculateMSDsCDFs_Marquee()"
			break
	EndSwitch
	
	wave /WAVE/Z tracksWave = GetCoordinatesFromMarquee("tracks", marqueeWin, xPixelSize, yPixelSize, xStart, xEnd, yStart, yEnd)
	if (WaveExists(tracksWave) == 0)
		// signals a user cancel
		return 0
	endif
	
	wave /WAVE selectedTracksWave = GetTracksWithinLimits(tracksWave, "allinbox", xStart, xEnd, yStart, yEnd)
	
	if (useMSD)
		CalculateMSDs_menu(tracksWave = selectedTracksWave, providedTracksWaveName = NameOfWave(tracksWave))
	else
		CalculateCDFs_menu(tracksWave = selectedTracksWave, providedTracksWaveName = NameOfWave(tracksWave))
	endif
End

Function /WAVE GetCoordinatesFromMarquee(positionsOrTracksString, marqueeWin, xPixelSize, yPixelSize, xStart, xEnd, yStart, yEnd)
	string positionsOrTracksString	// fetch pixel size from positions or from tracks?
	string &marqueeWin
	variable &xPixelSize, &yPixelSize, &xStart, &xEnd, &yStart, &yEnd
	
	// if positionsOrTracksString == "positions" then work on positions, if == "tracks" then work on tracks. Everything else is an error.
	// returns the coordinates of the marquee by reference in units of pixels, as well as the pixel sizes in nm
	// return value is that positions wave or tracks wave selected by the user
	
	variable workingOnPositions = 0
	StrSwitch (positionsOrTracksString)
		case "positions":
			workingOnPositions = 1
			break
		case "tracks":
			workingOnPositions = 0
			break
		default:
			Abort "Unknown type passed to GetCoordinatesFromMarquee()"
			break
	EndSwitch
	
	GetMarquee /Z left, bottom
	if (V_flag == 0)
		return $""
	endif
	
	marqueeWin = S_MarqueeWin
	variable bottom, top, left, right
	bottom = min(V_bottom, V_top)
	top = max(V_bottom, V_top)
	left = min(V_left, V_right)
	right = max(V_left, V_right)
	
	if ((top == bottom) || (left == right))	// it's not a real box
		return $""
	endif
	
	variable userCancelled
	if (workingOnPositions)
		// try to get the positions wave from the graph automatically
		string posName
		wave /Z pos = GetPositionsWaveFromGraph(marqueeWin, userCancelled)
		if (userCancelled)
			return $""
		endif
		if (WaveExists(pos) == 0)
			Abort "The selected positions wave doesn't seem to exist!"
		endif
	else
		// working on tracks
		wave /WAVE /Z tracksWave = GetTrackingWaveFromGraph(marqueeWin, userCancelled)
		if (userCancelled)
			return $""
		endif
		if (!WaveExists(tracksWave))
			Abort "The selected tracks wave doesn't seem to exist!"
		endif
	endif
	
	variable err
	string units
	string dataWaveNote
	if (workingOnPositions)
		dataWaveNote = note(pos)
	else
		dataWaveNote = note(tracksWave)
	endif
	
	// the main localizerviewer panel is an exception, since its coordinates
	// are always in units of pixels, so it doesn't bother with unit strings
	// and so forth. Handle this case separately
	if (GrepString(marqueeWin, "^LocalizerViewer[0-9]*#CCDViewer$"))
		xPixelSize = NumberByKey("X PIXEL SIZE", dataWaveNote)
		yPixelSize = NumberByKey("Y PIXEL SIZE", dataWaveNote)
		if ((NumType(xPixelSize) == 2) || (NumType(yPixelSize) == 2))
			xPixelSize = 0
			yPixelSize = 0
		endif
		xStart = left
		yStart = bottom
		xEnd = right
		yEnd = top
		if (workingOnPositions)
			return pos
		else
			return tracksWave
		endif
	endif
	
	// determine if the plot is in units of pixels or um
	string traceNames, imageNames, unitsString
	traceNames = TraceNameList("", ";", 1)
	imageNames = ImageNameList("", ";")
	if (ItemsInList(traceNames) > 0)	// the plot contains some traces, use these to get the scale
		wave traceWave = TraceNameToWaveRef("", StringFromList(0, traceNames))
		unitsString = WaveUnits(traceWave, 0)
	elseif (ItemsInList(imageNames) > 0) // the plot does not contain any traces, try using an image instead
		wave imageWave = ImageNameToWaveRef("", StringFromList(0, imageNames))
		unitsString = WaveUnits(imageWave, 0)
	else	// no traces or images on the plot
		Abort "The graph does not seem to contain any traces or images"
	endif
	
	// adjust the marquee coordinates depending on whether the graph
	// is in units of pixels or um
	strswitch (unitsString)
		case "pixel":
		case "":	// treat no units as pixels
			xPixelSize = 0
			yPixelSize = 0
			break
		case "um":
			// get the wave units
			xPixelSize = NumberByKey("X PIXEL SIZE", dataWaveNote)
			yPixelSize = NumberByKey("Y PIXEL SIZE", dataWaveNote)
			if ((NumType(xPixelSize) == 2) || (NumType(yPixelSize) == 2))	// NaN
				Abort "Unable to recognize the units of the plot."
			endif
			break
		default:
			Abort "Unable to recognize the units of the plot."
	endswitch
	
	if (xPixelSize != 0)
		// convert the box to pixels
		xStart = left / (xPixelSize / 1000)
		yStart = bottom / (yPixelSize / 1000)
		XEnd = right / (xPixelSize / 1000)
		yEnd = top / (yPixelSize / 1000)
	else
		// the box already is in units of pixels
		xStart = left
		yStart = bottom
		XEnd = right
		yEnd = top
	endif
	
	if (workingOnPositions)
		return pos
	else
		return tracksWave
	endif
End

Function ManuallyFilterTracks()
	DFREF packageFolder = root:Packages:Localizer

	string trackingWaveList = GetPossibleTrackingWaves()
	if (ItemsInList(trackingWaveList) == 0)
		Abort "No particle tracks seem to have been calculated"
	endif
	
	string trackingWaveName
	Prompt trackingWaveName, "Track wave:", popup, trackingWaveList
	DoPrompt "Select the tracks to filter", trackingWaveName
	if (V_flag != 0)
		return 0
	endif
	
	string /G packageFolder:S_manFilterTracksOutName = trackingWaveName + "_filt"
	
	wave /WAVE /Z trackWave = GetTrackingWaveReference(trackingWaveName)
	if (WaveExists(trackWave) == 0)
		Abort "Unable to find the requested track"
	endif
	
	Duplicate /O trackWave, packageFolder:W_TracksToFilter
	GenerateFilterScatterWaves(packageFolder:W_TracksToFilter)
	wave W_TracksFilterX = packageFolder:W_TracksFilterX
	wave W_TracksFilterY = packageFolder:W_TracksFilterY
	
	DoWindow /K ManuallyFilterTracksWindow
	NewPanel /N=ManuallyFilterTracksWindow /W=(232,58,1027,603)
	DoWindow /T ManuallyFilterTracksWindow, "Filter tracks from " + trackingWaveName
	DefineGuide GraphLeft={FL,0.02,FR}, GraphTop={FT,0.02,FB},GraphRight={FR,-200},GraphBottom={FT,0.98,FB}
	DefineGuide ControlLeft={GraphRight, 10}, ControlTop={GraphTop, 0}, ControlRight={FL,0.98,FR}, ControlBottom={GraphBottom,0}
	Display /N=FilteredTracksGraph /HOST=ManuallyFilterTracksWindow /FG=(GraphLeft,GraphTop,GraphRight,GraphBottom) W_TracksFilterY vs W_TracksFilterX
	Label /W=ManuallyFilterTracksWindow#FilteredTracksGraph left, "Distance (pixels)"
	Label /W=ManuallyFilterTracksWindow#FilteredTracksGraph bottom, "Distance (pixels)"
	NewPanel /N=FilterControlsPanel /HOST=ManuallyFilterTracksWindow /FG=(ControlLeft,ControlTop,ControlRight,ControlBottom)
	SetDrawLayer /W=ManuallyFilterTracksWindow#FilterControlsPanel UserBack
	DrawText /W=ManuallyFilterTracksWindow#FilterControlsPanel 11,55,"Position the mouse\rnear the track(s) to\rremove and type 'x'"
	DrawText /W=ManuallyFilterTracksWindow#FilterControlsPanel 13,111,"Brush size (pixels)"
	Button BTDoneFiltering,win=ManuallyFilterTracksWindow#FilterControlsPanel,pos={113,210},size={50,20},title="Done",proc=BTFinishManualTrackFilterProc
	Button BTCancelFiltering,win=ManuallyFilterTracksWindow#FilterControlsPanel,pos={14,210},size={50,20},title="Cancel",proc=BTCancelManualTrackFilterProc
	Slider SLFilterRadius,win=ManuallyFilterTracksWindow#FilterControlsPanel,pos={11,113},size={150,45},limits={1,50,1},value= 10,vert= 0
	SetWindow ManuallyFilterTracksWindow, hook(filterHook)=ManuallyFilterTracksHook
End

Function GenerateFilterScatterWaves(W_TracksToFilter)
	wave /WAVE W_TracksToFilter
	
	DFREF packageFolder = root:Packages:Localizer
	
	wave M_ScatterCoordinates = ParticleTrackToScatterPairs(W_TracksToFilter)
	Make /O/N=(DimSize(M_ScatterCoordinates, 0)) /D packageFolder:W_TracksFilterX, packageFolder:W_TracksFilterY
	wave W_TracksFilterX = packageFolder:W_TracksFilterX
	wave W_TracksFilterY = packageFolder:W_TracksFilterY
	W_TracksFilterX = M_ScatterCoordinates[p][0]
	W_TracksFilterY = M_ScatterCoordinates[p][1]
End

Function ManuallyFilterTracksHook(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	DFREF packageFolder = root:Packages:Localizer

	switch(s.eventCode)
		case 11:				// keyboard
			if ((s.keyCode == 120) || (s.keyCode == 88))	// x or X
				ControlInfo /W=ManuallyFilterTracksWindow#FilterControlsPanel SLFilterRadius
				variable radius = V_Value
				wave /WAVE W_TracksToFilter = packageFolder:W_TracksToFilter
				GetMouse /W=ManuallyFilterTracksWindow#FilteredTracksGraph
				wave W_TracksWave_Filtered = DeleteTracksAroundLocation(W_TracksToFilter, V_left, V_top, radius)
				Duplicate /O W_TracksWave_Filtered, packageFolder:W_TracksToFilter
				GenerateFilterScatterWaves(packageFolder:W_TracksToFilter)
				hookResult = 1
			endif
			break
	endswitch

	return hookResult		// 0 if nothing done, else 1
End

Function /WAVE DeleteTracksAroundLocation(W_TracksWave, pixelX, pixelY, radius)
	wave /WAVE W_TracksWave
	variable pixelX, pixelY, radius
	
	variable radius2 = radius^2
	variable xx = AxisValFromPixel("ManuallyFilterTracksWindow#FilteredTracksGraph", "Bottom", pixelX)
	variable yy = AxisValFromPixel("ManuallyFilterTracksWindow#FilteredTracksGraph", "Left", pixelY)
	
	variable nTracks = DimSize(W_TracksWave, 0)
	Make /FREE /N=(nTracks) /WAVE W_TracksWave_Filtered
	Note /K W_TracksWave_Filtered, note(W_TracksWave)
	variable i, j, offset = 0, xCol, yCol, zCol
	variable shouldSkip
	for (i = 0; i < nTracks; i+=1)
		wave thesePositions = W_TracksWave[i]
		variable nPositionsInTrack = DimSize(thesePositions, 0)
		GetColumnsForEmitterPositions(thesePositions, xCol, yCol, zCol)
		shouldSkip = 0
		for (j = 0; j < nPositionsInTrack; j+=1)
			if ((thesePositions[j][xCol] - xx)^2 + (thesePositions[j][yCol] - yy)^2 < radius2)
				shouldSkip = 1
				break
			endif
		endfor
		
		if (!shouldSkip)
			Duplicate /FREE thesePositions, thesePositions_Filtered
			W_TracksWave_Filtered[offset] = thesePositions_Filtered
			offset += 1
		endif
		
	endfor
	
	Redimension /N=(offset) W_TracksWave_Filtered
	return W_TracksWave_Filtered
End

Function BTCancelManualTrackFilterProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			DFREF packageFolder = root:Packages:Localizer
			DoWindow /K ManuallyFilterTracksWindow
			KillWaves /Z packageFolder:W_TracksFilterX, packageFolder:W_TracksFilterY, packageFolder:W_TracksToFilter
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTFinishManualTrackFilterProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			DFREF packageFolder = root:Packages:Localizer
			DFREF tracksFolder = root:'Particle Tracks'
			SVAR gManFilterTracksOutName = packageFolder:S_manFilterTracksOutName
			wave FilteredTracks = packageFolder:W_TracksToFilter
			string newTracksName = GetNewTrackingWaveName(gManFilterTracksOutName, tracksFolder, "Name of the filtered tracks wave:")
			if (strlen(newTracksName) == 0)
				KillWaves /Z packageFolder:W_TracksFilterX, packageFolder:W_TracksFilterY, packageFolder:W_TracksToFilter
				return 0
			endif
			DoWindow /K ManuallyFilterTracksWindow
			Duplicate /O FilteredTracks, tracksFolder:$newTracksName
			KillWaves /Z packageFolder:W_TracksFilterX, packageFolder:W_TracksFilterY, packageFolder:W_TracksToFilter
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function /S GetNewOutputWaveName(suggestedName, outputDataFolder, promptString)
	string suggestedName
	DFREF outputDataFolder
	string promptString
	// gets a new output wave name, even if there is another wave with the same name in 
	// outputDataFolder or anywhere else
	
	string outputWaveName = suggestedName
	variable outputWaveFound = 0
	
	do
		Prompt outputWaveName, "Wave name:"
		DoPrompt promptString, outputWaveName
		if (V_flag == 1)	// the user canceled
			return ""
		endif
		
		// check if the name is not a blank string
		if (strlen(outputWaveName) == 0)
			continue
		endif
		
		// check if the name is not reserved
		if ((StringMatch(outputWaveName, "* No Positions *") == 1) || (StringMatch(outputWaveName, "* No Tracks *") == 1) || (StringMatch(outputWaveName, "* No Calibrations *") == 1) || (StringMatch(outputWaveName, "* No Weights *") == 1))
			continue
		endif
		
		// check if the name is acceptable to Igor
		// do not allow liberal names
		if (stringmatch(outputWaveName, CleanupName(outputWaveName, 0)) != 1)
			outputWaveName = CleanupName(outputWaveName, 0)
			continue
		endif
		
		// check if the wave already exists
		wave /Z outputWave = outputDataFolder:$outputWaveName
		if (WaveExists(outputWave) == 1)
			// the requested wave already exists
			// if the pre-existing wave is a wave create by this package, then offer the user to overwrite it
			// otherwise require that the user specifies a new name
				if (IsPositionsWave(outputWave) || IsTrackingWave(outputWave) || IsRegistrationMap(outputWave) || IsAstigmatism3DCalibrationMap(outputWave) || IsSOFICombinationWeightsWave(outputWave))
				DoAlert 2, "A wave called \"" + outputWaveName + "\" already exists! Do you want to replace it with the new wave?"
				switch (V_flag)
					case 1:	// OK to replace
							// we're all set
						outputWaveFound = 1
						break
					case 2:	// Don't replace, use a different name instead
						outputWaveName = UniqueName(outputWaveName, 1, 0)
						continue
						break
					case 3:	// user cancelled, escape
						return ""
						break
					default:	// shouldn't happen
						Abort "Unknown value for V_flag in DoAlert"
						break
				endswitch
			else
				// don't overwrite a wave that is not a positions wave
				outputWaveName = UniqueName(outputWaveName, 1, 0)
				continue
			endif
		endif
		
		// if we arrive here then the waveName looks okay
		outputWaveFound = 1
	while (outputWaveFound == 0)
	
	return outputWaveName

End

Function /S GetNewPositionsWaveName(suggestedName, outputDataFolder, promptString)
	string suggestedName
	DFREF outputDataFolder
	string promptString
	
	do
		string outputWaveName = GetNewOutputWaveName(suggestedName, outputDataFolder, promptString)
		if (strlen(outputWaveName) == 0)
			return ""
		endif
		
		// the output wave does not exist within the requested data folder yet
		// however, it could still be name of a positions wave stored somewhere else
		if (IsUniquePositionsWaveName(outputWaveName) == 0)
			suggestedName = UniqueName(outputWaveName, 1, 0)
			continue
		endif
		
		break
	while (1)
	
	return outputWaveName
End

Function /S GetNewTrackingWaveName(suggestedName, outputDataFolder, promptString)
	string suggestedName
	DFREF outputDataFolder
	string promptString
	// prompt the user for a wave name in dataFolder
	// check that the name is a valid wave name and return it
	// if the user canceled then return an empty string
	// the promptstring is used in the dialogs to let the user know what happens
	
	do
		string outputWaveName = GetNewOutputWaveName(suggestedName, outputDataFolder, promptString)
		if (strlen(outputWaveName) == 0)
			return ""
		endif
		
		// the output wave does not exist within the requested data folder yet
		// however, it could still be name of a tracking wave stored somewhere else
		if (IsUniqueTrackingWaveName(outputWaveName) == 0)
			suggestedName = UniqueName(outputWaveName, 1, 0)
			continue
		endif
		
		// if we arrive here then the waveName looks okay
		break
	while (1)
	
	return outputWaveName
End

Function IsUniquePositionsWaveName(suggestedName)
	string suggestedName
	
	string listOfPositionsWaves = GetPossiblePositionsWaves()
	variable nWaves = ItemsInList(listOfPositionsWaves)
	
	string thisWaveName
	variable i
	for (i = 0; i < nWaves; i+=1)
		thisWaveName = StringFromList(i, listOfPositionsWaves)
		if (StringMatch(thisWaveName, suggestedName) == 1)
			return 0
		endif
	endfor
	
	return 1
End

Function IsUniqueTrackingWaveName(suggestedName)
	string suggestedName
	
	string listOfTrackingWaves = GetPossibleTrackingWaves()
	variable nWaves = ItemsInList(listOfTrackingWaves)
	
	string thisWaveName
	variable i
	for (i = 0; i < nWaves; i+=1)
		thisWaveName = StringFromList(i, listOfTrackingWaves)
		if (StringMatch(thisWaveName, suggestedName) == 1)
			return 0
		endif
	endfor
	
	return 1
End

Function /S GetNew3DCalibrationWaveName(suggestedName, outputDataFolder, promptString)
	string suggestedName
	DFREF outputDataFolder
	string promptString
	
	string outputWaveName = GetNewOutputWaveName(suggestedName, outputDataFolder, promptString)
	if (strlen(outputWaveName) == 0)
		return ""
	endif
	
	return outputWaveName
End

Function /S GetNewPixelCombinationsWaveName(suggestedName, outputDataFolder, promptString)
	string suggestedName
	DFREF outputDataFolder
	string promptString
	
	string outputWaveName = GetNewOutputWaveName(suggestedName, outputDataFolder, promptString)
	if (strlen(outputWaveName) == 0)
		return ""
	endif
	
	return outputWaveName
End

Function /WAVE MakeImageHistogram(image, windowName)
	wave image
	string windowName
	
	DFREF savDF = GetDataFolderDFR()
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	SetDataFolder windowDataFolder
	ImageHistogram image
	SetDataFolder savDF
	
	wave hist = windowDataFolder:W_ImageHist
	
	return hist
End

Function GetGUIFitSettingsAndOptions(fitParams, windowName)
	STRUCT FitData &fitParams
	string windowName
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	string baseWindowName = GetBaseWindowName(windowName)
	
	NVAR currentThreshold = windowDataFolder:V_currentThreshold
	NVAR currentImage = windowDataFolder:V_currentImage
	NVAR PFA = windowDataFolder:V_PFA
	NVAR smoothSigmaFactor = windowDataFolder:V_smoothSigmaFactor
	NVAR gaussianWidth = windowDataFolder:V_gaussianWidth
	NVAR nImages = windowDataFolder:V_nImages
	NVAR firstFrameInLocalization = windowDataFolder:V_firstFrameInLocalization
	NVAR lastFrameInLocalization = windowDataFolder:V_lastFrameInLocalization
	SVAR CCDFilePath = windowDataFolder:S_filePath
	
	variable fileType
	NVAR /Z gFileType = windowDataFolder:V_fileType
	if (!NVAR_Exists(gFileType))
		fileType = -1
	else
		fileType = gFileType
	endif
	
	NVAR particleVerifierOverlap = windowDataFolder:V_particleVerifierOverlap
	NVAR particleVerifierSymm = windowDataFolder:V_particleVerifierSymm
	NVAR particleVerifierEllipse = windowDataFolder:V_particleVerifierEllipse
	NVAR particleVerifierEllipseAstig = windowDataFolder:V_particleVerifierEllipseAstig
	
	variable thresholdMethod, preprocessing, postprocessing, particlefinder, localizationMethod
	
	ControlInfo /W=$baseWindowName#AnalysisControls PMThresholdMethod
	StrSwitch (S_Value)
		case "GLRT":
			thresholdMethod = THRESHOLD_METHOD_GLRT
			break
		case "Isodata":
			thresholdMethod = THRESHOLD_METHOD_ISODATA
			break
		case "Triangle":
			thresholdMethod = THRESHOLD_METHOD_TRIANGLE
			break
		case "Absolute":
			thresholdMethod = THRESHOLD_METHOD_DIRECT
			break
		case "SmoothSigma":
			thresholdMethod = THRESHOLD_METHOD_SMOOTHSIGMA
			break
		default:
			Abort "Unknown segmentation method"
			break
	EndSwitch
	
	ControlInfo /W=$baseWindowName#AnalysisControls PMThresholdPreprocessing
	StrSwitch (S_Value)
		case "None":
			preprocessing = PREPROCESSOR_NONE
			break
		case "3x3 Median Filter":
			preprocessing = PREPROCESSOR_3X3MEDIAN
			break
		case "5x5 Median Filter":
			preprocessing = PREPROCESSOR_5X5MEDIAN
			break
		case "1x1 Gaussian Smoothing":
			preprocessing = PREPROCESSOR_1X1GAUSSIAN
			break
		case "2x2 Gaussian Smoothing":
			preprocessing = PREPROCESSOR_2X2GAUSSIAN
			break
		case "3x3 Mean Filter":
			preprocessing = PREPROCESSOR_3X3MEAN
			break
		case "5x5 Mean Filter":
			preprocessing = PREPROCESSOR_5X5MEAN
			break
		default:
			Abort "Unknown preprocessing method"
			break
	EndSwitch
	
	ControlInfo  /W=$baseWindowName#AnalysisControls PMThresholdPostProcessing
	StrSwitch (S_Value)
		case "None":
			postprocessing =  POSTPROCESSOR_NONE
			break
		case "Remove Isolated Pixels":
			postprocessing =  POSTPROCESSOR_REMOVE_ISOLATED
			break
		default:
			Abort "Unknown postprocessing method"
			break
	EndSwitch
	
	ControlInfo /W=$baseWindowName#AnalysisControls PMParticleFinder
	StrSwitch (S_Value)
		case "4-way adjacency":
			particlefinder = PARTICLEFINDER_ADJACENT4
			break
		case "8-way adjacency":
			particlefinder = PARTICLEFINDER_ADJACENT8
			break
		default:
			Abort "Unknown particle verifier method"
			break
	EndSwitch
	
	ControlInfo /W=$baseWindowName#AnalysisControls PMLocalizationMethod
	StrSwitch (S_Value)
		case "Gaussian Fitting":
			localizationMethod = LOCALIZATION_GAUSS_FITTING
			break
		case "Gaussian Fitting (Fixed Width)":
			localizationMethod = LOCALIZATION_GAUSS_FITTING_FIX
			break
		case "Iterative Multiplication":
			localizationMethod = LOCALIZATION_MULTIPLICATION
			break
		case "Center-of-Mass":
			localizationMethod = LOCALIZATION_CENTROID
			break
		case "Ellipsoidal Gaussian Fitting":
			localizationMethod = LOCALIZATION_ELLIPSOIDAL2DGAUSS
			break
		case "MLEwG":
			localizationMethod = LOCALIZATION_MLEWG
			break
		case "Ellipsoidal Gaussian Fitting (astigmatism)":
			localizationMethod = LOCALIZATION_ELLIPSGAUSS_ASTIG
			break
		default:
			Abort "Unknown localization method in GetGUIFitSettingsAndOptions()"
	EndSwitch
	
	// do a bunch of ugly crap with the particle verification since it's pretty hard to procedurally generate a multi-keyword flag
	if (particleVerifierSymm != 0)
		fitParams.PVerSymm = PARTICLEVERIFIER_2DGAUSS
	else
		fitParams.PVerSymm = PARTICLEVERIFIER_NONE
	endif
	
	if (particleVerifierEllipse != 0)
		fitParams.PVerEllipseSymm = PARTICLEVERIFIER_ELLIPS_SYMM
	else
		fitParams.PVerEllipseSymm = PARTICLEVERIFIER_NONE
	endif
	
	if (particleVerifierOverlap != 0)
		fitParams.PVerOverlap = PARTICLEVERIFIER_OVERLAP
	else
		fitParams.PVerOverlap = PARTICLEVERIFIER_NONE
	endif
	
	if (particleVerifierEllipseAstig != 0)
		fitParams.PVerEllipseAstig = PARTICLEVERIFIER_ELLIPS_ASTIG
	else
		fitParams.PVerEllipseAstig = PARTICLEVERIFIER_NONE
	endif
	
	fitParams.directThresholdLevel = currentThreshold
	fitParams.currentImage = currentImage
	fitParams.PFA = PFA
	fitParams.smoothSigmaFactor = smoothSigmaFactor
	fitParams.PSFWidth = gaussianWidth
	fitParams.firstFrameToAnalyze = firstFrameInLocalization
	fitParams.lastFrameToAnalyze = lastFrameInLocalization
	fitParams.numberOfImages = nImages
	fitParams.thresholdMethod = thresholdMethod
	fitParams.preprocessing = preprocessing
	fitParams.postprocessing = postprocessing
	fitParams.particlefinder = particlefinder
	fitParams.localizationMethod = localizationMethod
	fitParams.CCDFilePath = CCDFilePath
	fitParams.cameraType = fileType
End

Function GenerateGUIThreshold(windowName)
	string windowName
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	DFREF savDF = GetDataFolderDFR()
	string baseWindowName = GetBaseWindowName(windowName)
	
	NVAR currentImage = windowDataFolder:V_currentImage
	variable showThreshold, showPositions, showFittedPositions, showParticleTracks, showTrackIdentifiers
	wave image = windowDataFolder:M_CCDImage
	string tracesInGraph, imagesInGraph
	
	// check if we need to display the segmentation
	ControlInfo /W=$baseWindowName#AnalysisControls CBThresholdCheck
	showThreshold = V_value
	
	// check if we need to display the particles
	ControlInfo /W=$baseWindowName#AnalysisControls CBShowParticles
	showPositions = V_value
	
	// check if we need to display fitted positions
	ControlInfo /W=$baseWindowName#AnalysisControls CBDisplayFittedPositions
	if (V_value == 1)
		// does the wave with the fitted positions exist?
		// we can only show the positions if they exist
		ControlInfo /W=$baseWindowName#AnalysisControls PMSelectPositionsWave
		
		wave /Z fittedPositions = GetPositionsWaveReference(S_Value)
		
		if (WaveExists(fittedPositions) == 0)
			showFittedPositions = 0
		else
			showFittedPositions = 1
		endif
	else
		showFittedPositions = 0
	endif
	
	// check if we need to display particle tracks
	ControlInfo /W=$baseWindowName#TrackingControls CBDisplayTracks
	showParticleTracks = V_Value
	variable showFullTrack, UseDisplacementColors
	if (showParticleTracks == 1)
		// does the wave with the particles exists?
		ControlInfo /W=$baseWindowName#TrackingControls PMSelectTrackWave
		wave /Z trackingWave = GetTrackingWaveReference(S_Value)
		if (WaveExists(trackingWave) == 0)
			showParticleTracks = 0
		endif
		
		// do we need the full track or not?
		ControlInfo /W=$baseWindowName#TrackingControls PMDisplayTracksMode
		StrSwitch (S_Value)
			case "Display tracks in entirety":
				showFullTrack= 1
				break
			case "Tracks are constructed with playback":
				showFullTrack = 0
				break
			case "Show all calculated tracks":
				showFullTrack = 2
				break
			default:
				Abort "Missing track mode"
				break
		EndSwitch
		
		// and the track identifiers?
		ControlInfo /W=$baseWindowName#TrackingControls CBDisplayTrackIdentifier
		showTrackIdentifiers = V_Value
		
		// color based on displacement or randomly?
		ControlInfo /W=$baseWindowName#TrackingControls PMTracksColorMode
		StrSwitch (S_Value)
			case "Arbitrary Track Colors":
				UseDisplacementColors = 0
				break
			case "Color tracks by average displacement":
				UseDisplacementColors = 1
				break
			case "Color each step individually by distance":
				UseDisplacementColors = 2
				break
			default:
				Abort "Missing track color mode"
				break
		EndSwitch
	endif
	
	// do we need to calculate the threshold or the positions?
	// if not then there is no reason to run the calculation
	if ((showThreshold == 1) || (showPositions == 1) || (showFittedPositions == 1))
		SetDataFolder windowDataFolder
		Struct FitData fitParams 
		GetGUIFitSettingsAndOptions(fitParams, baseWindowName)
		Generatethreshold(image, fittedPositions, fitParams, showThreshold, showPositions, showFittedPositions, currentImage)
		SetDataFolder savDF
	endif
	
	// check if we need to show tracks
	if (showParticleTracks == 1)
		SetDataFolder windowDataFolder
		wave /WAVE W_TrackingWaves = GenerateTrackingDisplayWaves(trackingWave, showFullTrack, showTrackIdentifiers, currentImage, UseDisplacementColors)
		SetDataFolder savDF
		wave M_Tracks = W_TrackingWaves[0]
		wave M_TrackingColors = W_TrackingWaves[1]
		wave W_TrackIndices = W_TrackingWaves[2]
		wave M_TrackIndicesLocation = W_TrackingWaves[3]
		wave /Z M_TrackIndicesColor = W_TrackingWaves[4]
	endif
	
	// what traces are currently displayed?
	tracesInGraph = TraceNameList(baseWindowName + "#CCDViewer",";",1)
	imagesInGraph = ImageNameList(baseWindowName + "#CCDViewer",";")
	
	// check if we need to show threshold traces
	if (showThreshold == 1)	// show the threshold
		// is the threshold already on the graph? If they are not on the graph, then show them
		if  (stringmatch(imagesInGraph,"*M_SegmentedImage*") != 1)
			AppendImage /W=$baseWindowName#CCDViewer windowDataFolder:M_SegmentedImage
			ModifyImage  /W=$baseWindowName#CCDViewer M_SegmentedImage  ctab= {0.5,1,RedWhiteGreen,1}, minRGB=NaN,maxRGB=0
		endif
		
	else		// don't show the threshold
		if  (stringmatch(imagesInGraph,"*M_SegmentedImage*") == 1)
			RemoveImage /W=$baseWindowName#CCDViewer M_SegmentedImage
		endif
	endif
	
	// check if we need to show positions
	if (showPositions == 1)	// show the particles
		
		// are the particles already on the graph? If not, show them
		if  (stringmatch(tracesInGraph,"*W_particles*") != 1)
			AppendToGraph /W=$baseWindowName#CCDViewer windowDataFolder:W_particles_Y vs windowDataFolder:W_particles_X
			ModifyGraph /W=$baseWindowName#CCDViewer mode(W_particles_Y)=3,msize(W_particles_Y)=4,rgb(W_particles_Y)=(0,0,65535)
			ModifyGraph /W=$baseWindowName#CCDViewer textMarker(W_particles_Y)=0
		endif
		
	else		// don't show the particles. If they are on the graph then remove them
		if  (stringmatch(tracesInGraph,"*W_particles*") == 1)
			RemoveFromGraph /W=$baseWindowName#CCDViewer W_particles_Y
		endif
	endif
	
	// check if we need to show the fitted positions
	if (showFittedPositions == 1)
		// are the fitted positions already on the graph? If not, show them
		if  (stringmatch(tracesInGraph,"*ps_temp*") != 1)
			AppendToGraph /W=$baseWindowName#CCDViewer windowDataFolder:ps_temp_Y vs windowDataFolder:ps_temp_X
			ModifyGraph /W=$baseWindowName#CCDViewer mode(ps_temp_Y)=3,marker(ps_temp_Y)=8, rgb(ps_temp_Y)=(65535,65535,0)
		endif
			
	else		// we need to remove the positions plot if it's on the graph
		if  (stringmatch(tracesInGraph,"*ps_temp*") == 1)
			RemoveFromGraph /W=$baseWindowName#CCDViewer $"ps_temp_Y"
		endif
	endif
	
	// check if we need to show the tracks
	if (showParticleTracks == 1)
		// are the fitted positions already on the graph? If not, show them
		if  (stringmatch(tracesInGraph,"*M_ActiveTracks*") != 1)
			AppendToGraph /W=$baseWindowName#CCDViewer M_Tracks[][1] vs M_Tracks[][0]
			ModifyGraph /W=$baseWindowName#CCDViewer mode(M_ActiveTracks)=4,marker(M_ActiveTracks)=19,msize(M_ActiveTracks)=1
			ModifyGraph /W=$baseWindowName#CCDViewer zColor(M_ActiveTracks)={M_TrackingColors,*,*,directRGB}
		endif
		if (showTrackIdentifiers == 1)
			if (stringMatch(tracesInGraph, "*M_TrackIndicesLocation*") == 0)
				AppendToGraph /W=$baseWindowName#CCDViewer M_TrackIndicesLocation[][1] vs M_TrackIndicesLocation[][0]
				ModifyGraph /W=$baseWindowName#CCDViewer zColor(M_TrackIndicesLocation)={M_TrackIndicesColor,*,*,directRGB}, mode(M_TrackIndicesLocation)=3
				ModifyGraph /W=$baseWindowName#CCDViewer msize(M_TrackIndicesLocation)=3,textMarker(M_TrackIndicesLocation)={W_TrackIndices,"default",0,0,5,5,5}
			endif
		else
			if (stringMatch(tracesInGraph, "*M_TrackIndicesLocation*") == 1)
				RemoveFromGraph /W=$baseWindowName#CCDViewer /Z M_TrackIndicesLocation
			endif
		endif
	else		// remove the tracking wave if its on the graph
		if  (stringmatch(tracesInGraph,"*M_ActiveTracks*") == 1)
			RemoveFromGraph /W=$baseWindowName#CCDViewer M_ActiveTracks
			RemoveFromGraph /W=$baseWindowName#CCDViewer /Z M_TrackIndicesLocation
		endif
	endif
		
	// explicitly update the display if we're running on Windows
	#if (StringMatch(IgorInfo(2), "Windows") == 1)
		DoUpdate
	#endif
	
End

Function Generatethreshold(image, fittedPositions, fp, showThreshold, showPositions, showFittedPositions, frameNumber)
	wave /Z image, fittedPositions
	struct FitData &fp
	variable showThreshold, showPositions, showFittedPositions, frameNumber
	
	// generate all waves corresponding to the emitter segmentation
	// right now this function requires that the active data folder is set to
	// one of the package subfolders
	
	Assert(frameNumber >= 0)
	Assert(WaveExists(image))
	
	EmitterSegmentation /M=(fp.thresholdMethod) /FM=(fp.localizationMethod) /G={fp.preprocessing, fp.postprocessing} /PFA=(fp.PFA) /SSG=(fp.smoothSigmaFactor) /ABS=(fp.directThresholdLevel) /WDTH=(fp.PSFWidth) /S=(showPositions) /F=(fp.particlefinder) /PVER={fp.PVerSymm, fp.PVerEllipseSymm, fp.PVerOverlap, fp.PVerEllipseAstig} /R=(fp.minDistanceBetweenParticles) image
	
	if (showFittedPositions != 0)
		GenerateFittedPositionsWaves(fittedPositions, frameNumber)
	endif
	
	if (showPositions != 0)
		MakeParticleWavesAndLabel()
	endif
End

Function /WAVE GenerateTrackingDisplayWaves(trackWave, showFullTrack, makeTrackIndices, frameNumber, colorByDisplacement)
	wave /WAVE trackWave
	variable showFullTrack, makeTrackIndices, frameNumber, colorByDisplacement
	
	NewDataFolder /O :TrackingDisplayWaves
	DFREF TrackingDisplayWavesFolder = :TrackingDisplayWaves
	
	// make a wave that will display the particle tracks active in the current frame
	// if showFullTrack == 0 then show the tracks that include frameNumber in its run
	// up to the frameNumber.
	// if showFullTrack == 1 then show the complete track for any track that includes
	// this frameNumber in its run.
	// if showFullTrack == 2 then show all calculated tracks.
	
	// if colorByDisplacement == 0 then each track gets a random color
	//	if colorByDisplacement == 1 then each track is colored according to its average step distance over the entire track
	// if colorByDisplacement == 2 then each individual step is colored according to its step distance
	
	// if makeTrackIndices != 0 then also make another set of waves, containing the
	// an additional 2D wave with x,y coordinates, a  numeric wave providing markers,
	// and a wave providing colors for each trace
	
	// the return from this function is a wave containing two wave references:
	// the first is a 2D wave containing the x and y coordinates of the tracks
	// each track is separated from the next by a row of NaN's
	// the second is a 2D wave suitable for the color as f(z) display, to color the tracks
	
	// this function is assumed to be called within the the relevant root:Packages:Localizer:<windowFolder>
	// for the requested window
	
	// we need to know how to get the x and y data out of the positions
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(trackWave, xCol, yCol, zCol)
	
	variable nTracks = DimSize(trackWave, 0)
	
	// waves that will display the track positions
	Make /O/N=(0, 2) /D trackingDisplayWavesFolder:M_ActiveTracks
	wave M_ActiveTracks = trackingDisplayWavesFolder:M_ActiveTracks
	
	// wave containing colors for the tracks
	Make /O/N=(0, 3)/W/U TrackingDisplayWavesFolder:M_TrackColors
	wave M_TrackColors = TrackingDisplayWavesFolder:M_TrackColors
	Make /FREE/N=(1,3) /D M_CurrentColor
	
	// extract tracks that contain the present frame if needed
	if (showFullTrack != 2)
		wave /wave W_ExtractionResult = ExtractTracksInFrame(trackWave, frameNumber)
		wave /wave W_ExtractedTracks = W_ExtractionResult[0]
		wave W_ExtractedTracksIndices = W_ExtractionResult[1]
	else
		wave /wave W_ExtractedTracks = trackWave
		Make /FREE /N=(nTracks) W_ExtractedTracksIndices = p
	endif
	variable nTracksInFrame = DimSize(W_ExtractedTracks, 0)
	
	// waves for the track indices
	Make /O/N=(nTracksInFrame) /I/U trackingDisplayWavesFolder:W_TrackIndices
	Make /O/N=(nTracksInFrame, 2) /D trackingDisplayWavesFolder:M_TrackIndicesLocation
	Make /O/N=(nTracksInFrame, 3)/W/U trackingDisplayWavesFolder:M_TrackIndicesColor
	wave W_TrackIndices = trackingDisplayWavesFolder:W_TrackIndices
	wave M_TrackIndicesLocation = trackingDisplayWavesFolder:M_TrackIndicesLocation
	wave M_TrackIndicesColor = trackingDisplayWavesFolder:M_TrackIndicesColor
	
	// colors to use (accessed using a hash function)
	ColorTab2Wave Rainbow
	wave M_Colors = M_Colors
	variable nColors = DimSize(M_Colors, 0)
	
	variable minDisplacement, maxDisplacement
	if (colorByDisplacement == 1)
		// color by average displacement
		wave W_TrackDisplacements = EstimatePerTrackDisplacements(trackWave)
		WaveStats /Q W_TrackDisplacements
		minDisplacement = V_min
		maxDisplacement = V_max
	elseif (colorByDisplacement == 2)
		// color each step individually
		GetMinAndMaxDisplacements(trackWave, minDisplacement, maxDisplacement)
	endif
	
	variable nPositionsInTrack, thisTrackIndex
	variable colorIndex, offset
	variable i, j, stepIndex
	for (i = 0; i < nTracksInFrame; i+=1)
		wave currentTrack = W_ExtractedTracks[i]
		thisTrackIndex = W_ExtractedTracksIndices[i]
		nPositionsInTrack = DimSize(currentTrack, 0)
		if (nPositionsInTrack > DimSize(M_CurrentColor, 0))
			Redimension /N=(nPositionsInTrack, -1) M_CurrentColor
		endif
		getColumnsForEmitterPositions(currentTrack, xCol, yCol, zCol)
		
		// color that this track will have
		variable thisStepDistance
		if (colorByDisplacement == 0)
			colorIndex = mod(StringCRC(0, num2str(currentTrack[0][0] + currentTrack[0][xCol] + currentTrack[0][yCol] * 1e6)), nColors)	// arbitrary semi-unique identifier
			M_CurrentColor = M_Colors[colorIndex][q]
		elseif (colorByDisplacement == 1)
			colorIndex = limit((W_TrackDisplacements[thisTrackIndex] - minDisplacement) / (maxDisplacement - minDisplacement) * nColors, 0, nColors - 1)
			M_CurrentColor = M_Colors[colorIndex][q]
		elseif (colorByDisplacement == 2)
			for (stepIndex = 0; stepIndex < (nPositionsInTrack - 1); stepIndex += 1)
				thisStepDistance = sqrt((currentTrack[stepIndex + 1][xCol] - currentTrack[stepIndex][xCol])^2 + (currentTrack[stepIndex + 1][yCol] - currentTrack[stepIndex][yCol])^2)
				colorIndex = limit((thisStepDistance - minDisplacement) / (maxDisplacement - minDisplacement) * nColors, 0, nColors - 1)
				M_CurrentColor[stepIndex][] = M_Colors[colorIndex][q]
			endfor
		endif
		
		// set up an entry in the track indices waves
		W_TrackIndices[i] = thisTrackIndex
		M_TrackIndicesLocation[i][0] = currentTrack[0][xCol]
		M_TrackIndicesLocation[i][1] = currentTrack[0][yCol]
		M_TrackIndicesColor[i][] = M_CurrentColor[0][q]
		
		// the number of points that we have to include in the plot depends on whether the full track
		// should be shown, or whether it should be constructed as we step through the movie
		variable firstPositionIndex = 0, lastPositionIndex, nPositionsToInclude
		if (showFullTrack)
			lastPositionIndex = nPositionsInTrack - 1
		else
			for (j = 0; j < nPositionsInTrack; j+=1)
				if (currentTrack[j][0] > frameNumber)
					break
				endif
			endfor
			lastPositionIndex = j - 1
		endif
		nPositionsToInclude = lastPositionIndex - firstPositionIndex + 1
		offset = DimSize(M_ActiveTracks, 0)
		Redimension /N=(DimSize(M_ActiveTracks, 0) + nPositionsToInclude + 1, -1) M_ActiveTracks
		M_ActiveTracks[offset, offset + nPositionsToInclude - 1][0] = currentTrack[p - offset][xCol]
		M_ActiveTracks[offset, offset + nPositionsToInclude - 1][1] = currentTrack[p - offset][yCol]
		Redimension /N=(DimSize(M_TrackColors, 0) + nPositionsToInclude + 1, -1) M_TrackColors
		M_TrackColors[offset, offset + nPositionsToInclude - 1][] = M_CurrentColor[p - offset][q]
		
		// finish up each entry with a row of NaN's
		M_ActiveTracks[offset + nPositionsToInclude][] = NaN
		M_TrackColors[offset + nPositionsToInclude][] = NaN
		
	endfor
	
	KillWaves /Z M_Colors
	
	Make /N=(5) /WAVE /FREE W_Result
	W_Result[0] = M_ActiveTracks
	W_Result[1] = M_TrackColors
	W_Result[2] = W_TrackIndices
	W_Result[3] = M_TrackIndicesLocation
	W_Result[4] = M_TrackIndicesColor
	return W_Result
End

Function DeleteThisPos_TracePopup()
	GetLastUserMenuInfo
	string baseWindowName = GetBaseWindowName(S_graphName)
	string traceName = S_traceName
	
	if (StringMatch(baseWindowName, "LocalizerViewer*") != 1)
		Abort "This menu should only be called on one of the viewers created by the localizer package"
	endif
	if (StringMatch(traceName, "ps_temp_Y") != 1)
		Abort "This menu should only be called on the trace showing particle tracks"
	endif
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	NVAR currentImageIndex = windowDataFolder:V_currentImage
	NVAR rightMouseH = windowDataFolder:V_RightMouseH
	NVAR rightMouseV = windowDataFolder:V_RightMouseV
	
	ControlInfo /W=$baseWindowName#AnalysisControls PMSelectPositionsWave
	string posName = S_value
	wave /z positions = GetPositionsWaveReference(posName)
	if (!WaveExists(positions))
		Abort "Unable to find the selected positions"
	endif
	
	variable xx = AxisValFromPixel(baseWindowName + "#CCDViewer", "bottom", rightMouseH)
	variable yy = AxisValFromPixel(baseWindowName + "#CCDViewer", "left", rightMouseV)
	
	wave /wave W_ExtractedPositions = ExtractPositionsInFrame(positions, currentImageIndex)
	wave M_emittersInCurrentFrame = W_ExtractedPositions[0]
	wave M_PositionsIndices = W_ExtractedPositions[1]
	
	variable lastX, lastY, nearestX, nearestY, index
	variable err = returnPointNearestFitPositions(M_emittersInCurrentFrame, xx, yy, nearestX, nearestY, index)
	if (err != 0)
		return 0
	endif
	
	variable distance = sqrt((xx - nearestX)^2 + (yy - nearestY)^2)
	if (distance > 5)
		return 0
	endif
	
	DeletePoints /M=0 M_PositionsIndices[index], 1, positions
	GenerateGUIThreshold(baseWindowName)
End

Function CombineTracks_TracePopup()
	GetLastUserMenuInfo
	string baseWindowName = GetBaseWindowName(S_graphName)
	string traceName = S_traceName
	
	if (StringMatch(baseWindowName, "LocalizerViewer*") != 1)
		Abort "This menu should only be called on one of the viewers created by the localizer package"
	endif
	
	if ((StringMatch(traceName, "M_TrackIndicesLocation") != 1) && (StringMatch(traceName, "M_ActiveTracks") != 1))
		Abort "This menu should only be called on the trace showing particle tracks"
	endif
	
	// get the name of the track wave
	ControlInfo /W=$baseWindowName#TrackingControls PMSelectTrackWave
	string trackWaveName = S_Value
	wave /WAVE /Z W_TrackingWave = GetTrackingWaveReference(trackWaveName)
	if (WaveExists(W_TrackingWave) == 0)
		Abort "Unable the find the selected tracking wave"
	endif
	
	variable nTracks = DimSize(W_TrackingWave, 0)
	
	// ask the user for which tracks to combine
	variable track1Index, track2Index
	Prompt track1Index, "First track identifier:"
	Prompt track2Index, "Second track identifier:"
	DoPrompt "Select the tracks to combine", track1Index, track2Index
	if (V_flag != 0)
		return 0
	endif
	
	// some constraints need to be enforced
	if ((track1Index < 0) || (track2Index < 0) || (track1Index >= track2Index))
		Abort "Invalid track indices"
	endif
	if ((track1Index >= nTracks) || (track2Index >= nTracks))
		Abort "Invalid track indices, only " + num2str(nTracks) + " tracks are present"
	endif
	
	wave track1 = W_TrackingWave[track1Index]
	wave track2 = W_TrackingWave[track2Index]
	
	// the first position of the second track needs to come after this one
	if (track1[DimSize(track1, 0) - 1][0] >= track2[0][0])
		Abort "The second track needs to start after the first one is finished"
	endif
	
	// looks like we're all set
	// replace track1 with the combination of the emitters
	// and delete track2
	variable nInitialEmittersInTrack1 = DimSize(track1, 0)
	Redimension /N=(nInitialEmittersInTrack1 + DimSize(track2, 0), -1) track1
	track1[nInitialEmittersInTrack1, ][] = track2[p - nInitialEmittersInTrack1][q]
	
	DeletePoints track2Index, 1, W_TrackingWave
	
	// update the display
	GenerateGuiThreshold(baseWindowName)
End

Function SplitTracks_TracePopup()
	GetLastUserMenuInfo
	string baseWindowName = GetBaseWindowName(S_graphName)
	string traceName = S_traceName
	
	if (StringMatch(baseWindowName, "LocalizerViewer*") != 1)
		Abort "This menu should only be called on one of the viewers created by the localizer package"
	endif
	
	if ((StringMatch(traceName, "M_TrackIndicesLocation") != 1) && (StringMatch(traceName, "M_ActiveTracks") != 1))
		Abort "This menu should only be called on the trace showing particle tracks"
	endif
	
	// get the name of the track wave
	ControlInfo /W=$baseWindowName#TrackingControls PMSelectTrackWave
	string trackWaveName = S_Value
	wave /WAVE /Z W_TrackingWave = GetTrackingWaveReference(trackWaveName)
	if (WaveExists(W_TrackingWave) == 0)
		Abort "Unable the find the selected tracking wave"
	endif
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	NVAR currentImageIndex = windowDataFolder:V_currentImage
	NVAR rightMouseH = windowDataFolder:V_RightMouseH
	NVAR rightMouseV = windowDataFolder:V_RightMouseV
	
	variable xx = AxisValFromPixel(baseWindowName + "#CCDViewer", "bottom", rightMouseH)
	variable yy = AxisValFromPixel(baseWindowName + "#CCDViewer", "left", rightMouseV)
	
	// get the track that is closest to the requested position
	ControlInfo /W=$baseWindowName#TrackingControls RBDisplayAllKnownTracks
	variable limitToCurrentImage = !V_value
	variable trackIndex, posIndexInTrack, nearestX, nearestY
	GetTrackNearestPosition(W_TrackingWave, (limitToCurrentImage ? currentImageIndex : -1), xx, yy, nearestX, nearestY, trackIndex, posIndexInTrack)
	if (trackIndex == -1)
		return 0
	endif
	
	variable distance = sqrt((xx - nearestX)^2 + (yy - nearestY)^2)
	if (distance > 5)
		return 0
	endif
	
	// we now create the two tracks such that the point the user clicked on is the start of the new track
	wave M_TrackToSplit = W_TrackingWave[trackIndex]
	variable nPosInOldTrack = DimSize(M_TrackToSplit, 0)
	variable nPosBefore = posIndexInTrack
	variable nPosAfter = nPosInOldTrack - posIndexInTrack
	variable deleteFirstPart = 0, deleteSecondPart = 0
	if ((nPosBefore < 2) && (nPosAfter < 2))
		DoAlert 1, "Splitting the track here does not leave any valid tracks. Delete the track instead?"
		if (V_flag == 1)
			// delete the track instead
			deleteFirstPart = 1
			deleteSecondPart = 1
		else
			return 0
		endif
	elseif (nPosBefore < 2)
		DoAlert 1, "Splitting the track here does not leave a valid first track. Delete the first part?"
		if (V_flag == 1)
			deleteFirstPart = 1
		else
			return 0
		endif
	elseif (nPosAfter < 2)
		DoAlert 1, "Splitting the track here does not leave a valid second track. Delete the second part?"
		if (V_flag == 1)
			deleteSecondPart = 1
		else
			return 0
		endif
	endif
	
	// divide the track
	Duplicate /FREE M_TrackToSplit, M_TrackBefore, M_TrackAfter
	DeletePoints /M=0 nPosBefore, nPosAfter, M_TrackBefore
	DeletePoints /M=0 0, nPosBefore, M_TrackAfter
	
	// insert the first track, or delete it
	if (!deleteFirstPart)
		W_TrackingWave[trackIndex] = M_TrackBefore
	else
		DeletePoints trackIndex, 1, W_TrackingWave
	endif
	
	if (!deleteSecondPart)
		// now look for the correct place to insert the second track
		// such that the sorting is maintained
		variable nTracks = DimSize(W_TrackingWave, 0)
		variable firstFrameInTrack2 = M_TrackAfter[0][0]
		variable i
		for (i = 0; i < nTracks; i+=1)
			wave currentEmitter = W_TrackingWave[i]
			if (currentEmitter[0][0] > firstFrameInTrack2)
				break
			endif
		endfor
		
		InsertPoints i, 1, W_TrackingWave
		W_TrackingWave[i] = M_TrackAfter
	endif
	
	// update the display
	GenerateGuiThreshold(baseWindowName)
End

Function DeleteTrackPosition_Popup()
	GetLastUserMenuInfo
	string baseWindowName = GetBaseWindowName(S_graphName)
	string traceName = S_traceName
	
	if (StringMatch(baseWindowName, "LocalizerViewer*") != 1)
		Abort "This menu should only be called on one of the viewers created by the localizer package"
	endif
	
	if ((StringMatch(traceName, "M_TrackIndicesLocation") != 1) && (StringMatch(traceName, "M_ActiveTracks") != 1))
		Abort "This menu should only be called on the trace showing particle tracks"
	endif
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	NVAR currentImageIndex = windowDataFolder:V_currentImage
	NVAR rightMouseH = windowDataFolder:V_RightMouseH
	NVAR rightMouseV = windowDataFolder:V_RightMouseV
	
	// get the name of the track wave
	ControlInfo /W=$baseWindowName#TrackingControls PMSelectTrackWave
	string trackWaveName = S_Value
	wave /WAVE /Z W_TrackingWave = GetTrackingWaveReference(trackWaveName)
	if (WaveExists(W_TrackingWave) == 0)
		Abort "Unable the find the selected tracking wave"
	endif
	
	variable xx = AxisValFromPixel(baseWindowName + "#CCDViewer", "bottom", rightMouseH)
	variable yy = AxisValFromPixel(baseWindowName + "#CCDViewer", "left", rightMouseV)
	
	// get the track that is closest to the requested position
	ControlInfo /W=$baseWindowName#TrackingControls RBDisplayAllKnownTracks
	variable limitToCurrentImage = !V_value
	variable trackIndex, posIndexInTrack, nearestX, nearestY
	GetTrackNearestPosition(W_TrackingWave, (limitToCurrentImage ? currentImageIndex : -1), xx, yy, nearestX, nearestY, trackIndex, posIndexInTrack)
	if (trackIndex == -1)
		return 0
	endif
	
	variable distance = sqrt((xx - nearestX)^2 + (yy - nearestY)^2)
	if (distance > 5)
		return 0
	endif
	// all set, just remove the requested position from the track
	wave M_MatchingTrack = W_TrackingWave[trackIndex]
	DeletePoints /M=0 posIndexInTrack, 1, M_MatchingTrack
	
	// update the display
	GenerateGuiThreshold(baseWindowName)
End

Function DeleteTrack_TracePopup()
	GetLastUserMenuInfo
	string baseWindowName = GetBaseWindowName(S_graphName)
	string traceName = S_traceName
	
	if (StringMatch(baseWindowName, "LocalizerViewer*") != 1)
		Abort "This menu should only be called on one of the viewers created by the localizer package"
	endif
	
	if ((StringMatch(traceName, "M_TrackIndicesLocation") != 1) && (StringMatch(traceName, "M_ActiveTracks") != 1))
		Abort "This menu should only be called on the trace showing particle tracks"
	endif
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	NVAR currentImageIndex = windowDataFolder:V_currentImage
	NVAR rightMouseH = windowDataFolder:V_RightMouseH
	NVAR rightMouseV = windowDataFolder:V_RightMouseV
	
	// get the name of the track wave
	ControlInfo /W=$baseWindowName#TrackingControls PMSelectTrackWave
	string trackWaveName = S_Value
	wave /WAVE /Z W_TrackingWave = GetTrackingWaveReference(trackWaveName)
	if (WaveExists(W_TrackingWave) == 0)
		Abort "Unable the find the selected tracking wave"
	endif
	
	variable xx = AxisValFromPixel(baseWindowName + "#CCDViewer", "bottom", rightMouseH)
	variable yy = AxisValFromPixel(baseWindowName + "#CCDViewer", "left", rightMouseV)
	
	// get the track that is closest to the requested position
	ControlInfo /W=$baseWindowName#TrackingControls RBDisplayAllKnownTracks
	variable limitToCurrentImage = !V_value
	variable trackIndex, posIndexInTrack, nearestX, nearestY
	GetTrackNearestPosition(W_TrackingWave, (limitToCurrentImage ? currentImageIndex : -1), xx, yy, nearestX, nearestY, trackIndex, posIndexInTrack)
	if (trackIndex == -1)
		return 0
	endif
	
	variable distance = sqrt((xx - nearestX)^2 + (yy - nearestY)^2)
	if (distance > 5)
		return 0
	endif
	// all set, just remove the track
	DeletePoints trackIndex, 1, W_TrackingWave
	
	// update the display
	GenerateGuiThreshold(baseWindowName)
End
	
Function GenerateFittedPositionsWaves(PositionsWave, n)
	wave PositionsWave
	variable n
	
	Assert(WaveExists(PositionsWave))
	Assert(n >= 0)
	
	variable numberOfPositions, xCol, yCol, zCol
	
	wave /WAVE M_ExtractedPositions = ExtractPositionsInFrame(PositionsWave, n)
	wave positions = M_ExtractedPositions[0]
	
	
	numberOfPositions = DimSize(positions, 0)
	
	if (numberOfPositions == 0)	// no particles were found in this image
									// due to the (perhaps rather poor) design of the GUI
									// we need to create the waves anyway
									// we give them a size of 1 point, set to NaN
		Make /D/O/N=(numberOfPositions) ps_temp_X = NaN
		Make /D/O/N=(numberOfPositions) ps_temp_Y = NaN
	endif
	
	// get the indices of the columns containing the x and y coordinates
	getColumnsForEmitterPositions(PositionsWave, xCol, yCol, zCol)
	if ((xCol == -1) || (yCol == -1))
		Abort "The positions passed to GenerateFittedPositionsWaves() do not appear to contain any (x,y) information"
	endif
	
	// are these positions subsets starting from a particular origin?
	variable xStart, yStart
	if (NumType(NumberByKey("X START", note(PositionsWave))) != 2)
		xStart = NumberByKey("X START", note(PositionsWave))
	else
		xStart = 0
	endif
	
	if (NumType(NumberByKey("Y START", note(PositionsWave))) != 2)
		yStart = NumberByKey("Y START", note(PositionsWave))
	else
		yStart = 0
	endif
	
	Make /D/O/N=(numberOfPositions) ps_temp_X = positions[p][xCol] + xStart
	Make /D/O/N=(numberOfPositions) ps_temp_Y = positions[p][yCol] + yStart
End

Function MakeParticleWavesAndLabel()
	wave M_locatedParticles
	
	Assert(WaveExists(M_LocatedParticles))
	
	variable nParticles = DimSize(M_locatedParticles, 0)
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(M_LocatedParticles, xCol, yCol, zCol)
	
	Make /D/O/N=(nParticles) W_particles_X
	Make /D/O/N=(nParticles) W_particles_Y
	
	W_particles_X = M_locatedParticles[p][xCol]
	W_particles_Y = M_locatedParticles[p][yCol]
End

Function DoPALMFitting(fp, windowName, posName, abortOnError)
	STRUCT FitData &fp
	string windowName
	string &posName
	variable abortOnError	// if nonzero then allow LocalizationAnalysis to abort with an error
							// if zero then do not cause aborts, but report a possible error as a nonzero return value
	
	string baseWindowName = GetBaseWindowName(windowName)
	DFREF wDF = GetWindowDataFolder(baseWindowName)
	
	NVAR pixelSize = wDF:V_CCDPixelSize
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct pStruct
	pStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype pStruct.func = ProgressReporterFunc
	
	// fit the data but also add the waveNote
	string waveNote
	variable err
	
	// make an output name for the position
	posName = GetBasePosName(wDF, fp.CCDFilePath)
	
	variable cType = fp.cameraType	// stupid crap begause Igor rejects command lines that are too long
	variable lMeth = fp.localizationMethod
	variable threshMeth = fp.thresholdMethod
	variable minDistance = fp.minDistanceBetweenParticles
	variable firstFrame = fp.firstFrameToAnalyze
	variable lastFrame = fp.lastFrameToAnalyze
	if (abortOnError != 0)
		LocalizationAnalysis /Y=(cType)/M=(lMeth)/D=(threshMeth)/G={fp.preprocessing, fp.postprocessing}/F=(fp.particlefinder)/PVER={fp.PVerSymm,fp.PVerEllipseSymm,fp.PVerEllipseAstig,fp.PVerOverlap}/RNG={firstFrame,lastFrame}/PFA=(fp.PFA)/SSG=(fp.smoothSigmaFactor)/T=(fp.directThresholdLevel)/R=(minDistance)/W=(fp.PSFWidth)/PROG=pStruct/DEST=wDF:$posName fp.CCDFilePath
	else	// don't abort for errors, but return them as an error code (use the /Z flag for LocalizationAnalysis)
		LocalizationAnalysis /Y=(cType)/M=(lMeth)/D=(threshMeth)/G={fp.preprocessing, fp.postprocessing}/F=(fp.particlefinder)/PVER={fp.PVerSymm,fp.PVerEllipseSymm,fp.PVerEllipseAstig,fp.PVerOverlap}/RNG={firstFrame,lastFrame}/PFA=(fp.PFA)/SSG=(fp.smoothSigmaFactor)/T=(fp.directThresholdLevel)/R=(minDistance)/W=(fp.PSFWidth)/Z/PROG=pStruct/DEST=wDF:$posName fp.CCDFilePath
		err = V_flag
	endif
	
	// check if the output wave exists
	// it might not exist if the user immediately aborted the procedure
	// since it is by request, do not treat a user abort as an error
	wave /Z POS_out = wDF:$posName
	
	if (WaveExists(POS_out) != 0)
		// append a note with information on the fit
		// most of the information is now appended by the XOP
		waveNote = ""
		
		if (pixelSize != 0)
			waveNote += "X PIXEL SIZE:" + num2str(pixelSize) + ";"
			waveNote += "Y PIXEL SIZE:" + num2str(pixelSize) + ";"
		endif
		waveNote += "CALCULATION DATE:" + date() + ";"
		
		Note /NOCR POS_out, waveNote
	endif
	
	if (abortOnError == 0)
		return err
	else
		return 0
	endif
End

Function /S GetBasePosName(inDataFolder, filePath)
	DFREF inDataFolder
	string filePath
	
	// get a temp name to use for the manual fit
	// base it on the file name, but add a suffix ("_pos")
	// and make sure it is a legal name
	string macFilePath = ParseFilePath(5, filePath, ":", 0, 1)
	string fileName = ParseFilePath(3, macFilePath, ":", 0, 1)
	
	string posName = GetValidWaveNameFromFileName(fileName, "_pos")
	
	// if there already is a positions wave with the same name, then we need to change the proposed name
	// unless that wave is in the same folder (i.e. we've updated the analysis settings and now
	// the earlier results need to be overwritten)
	if (IsUniquePositionsWaveName(posName) == 0)
		if (WaveExists(inDataFolder:$posName) && IsPositionsWave(inDataFolder:$posName))
			// posName is fine, we will overwrite a previous positions wave
		else
			posName = UniqueName(posName, 1, 0)
		endif
	endif
	
	return posName
End

Function /S GetValidWaveNameFromFileName(fileName, suffix)
	string fileName, suffix
	
	variable suffixLength = strlen(suffix)
	string wName
	if (strlen(fileName) + suffixLength < kMaxWaveName)
		wName =  CleanUpName(fileName + suffix, 0)
	else
		variable partLength = floor((kMaxWaveName - 2 - suffixLength) / 2)
		fileName = fileName[0, partLength] + "_" + fileName[strlen(fileName) - partLength, strlen(fileName) - 1] + suffix
		wName = Cleanupname(fileName, 0)
	endif
	
	return wName
End

Function /WAVE GetPositionsWaveReference(posName)
	string posName
	
	Assert(strlen(posName) > 0)
	
	// it's possible that the positions name is quoted
	// so take care of that by removing the quotes
	string noQuoteName = ReplaceString("'", posName, "")
	
	// look for a positions wave with the given name
	// only look in root:'Localized Positions' and in
	// the package folder, but do so recursively
	DFREF packageFolder = root:Packages:Localizer
	DFREF positionsFolder = root:'Localized Positions'
	
	string wavesInPackageFolder = RecursiveWaveList(packageFolder, noQuoteName, "DIMS:2")
	string wavesInPositionsFolder = RecursiveWaveList(positionsFolder, noQuoteName, "DIMS:2")
	
	// if the same positions wave appears twice then we'll keep that as an error for now
	Assert((ItemsInList(wavesInPackageFolder) == 0) || (ItemsInList(wavesInPositionsFolder) == 0))
	Assert(ItemsInList(wavesInPackageFolder) <= 1)
	Assert(ItemsInList(wavesInPositionsFolder) <= 1)
	
	if (ItemsInList(wavesInPackageFolder) == 1)
		return $StringFromList(0, wavesInPackageFolder)
	else
		return $StringFromList(0, wavesInPositionsFolder)
	endif
End

Function /WAVE GetTrackingWaveReference(trackName)
	string trackName
	
	Assert(strlen(trackName) > 0)
	
	// it's possible that the positions name is quoted
	// so take care of that by removing the quotes
	string noQuoteName = ReplaceString("'", trackName, "")
	
	// look for a tracking wave with the given name
	// only look in root:'Particle Tracks', but do so recursively
	DFREF trackFolder = root:'Particle Tracks'
	
	string wavesInTrackFolder = RecursiveWaveList(trackFolder, noQuoteName, "DIMS:1,BYTE:0,INTEGER:0,SP:0,DP:0,CMPLX:0,TEXT:0")
	
	// if the same positions wave appears twice then we'll keep that as an error for now
	Assert(ItemsInList(wavesInTrackFolder) <= 1)
	
	return $StringFromList(0, wavesInTrackFolder)
End

Function /WAVE GetPosOrTrackWaveFromGraph(posOrTracksString, windowName, userCancelled)
	string posOrTracksString, windowName
	variable &userCancelled
	// assume that the specified window contains an image or plot created with this package.
	// Try to get the positions wave that is referenced in this graph.
	
	userCancelled = 0
	variable wantPositions
	StrSwitch (posOrTracksString)
		case "positions":
			wantPositions = 1
			break
		case "tracks":
			wantPositions = 0
			break
		default:
			Abort "Unknown type passed to GetPosOrTrackWaveFromGraph()"
			break
	EndSwitch
	// we start by getting a special case out of the way: if we're looking at the CCD viewer
	// then we can simply look at the relevant control.
	string baseWindowName = GetBaseWindowName(windowName)
	if (GrepString(baseWindowName, "^LocalizerViewer[0-9]*$") == 1)
		// this case is not very clear-cut
		// for now we assume that the positions wave of interest
		// is the one selected in the popupmenu
		if (wantPositions)
			ControlInfo /W=$baseWindowName#AnalysisControls PMSelectPositionsWave
			return GetPositionsWaveReference(S_Value)
		else
			ControlInfo /W=$baseWindowName#TrackingControls PMSelectTrackWave
			return GetTrackingWaveReference(S_Value)
		endif
	endif
	
	// try to determine a positions wave automatically based on the paths of the waves that are shown
	do
		// get a list of traces in the graph
		string traceNames, imageNames
		traceNames = TraceNameList("", ";", 1)
		imageNames = ImageNameList("", ";")
		// how many traces and images do we have?
		variable nTraces = ItemsInList(traceNames), i
		variable nImages = ItemsInList(imageNames)
		if ((nTraces == 0) && (nImages == 0))
			break
		endif
		
		// get the the full wave paths (names including the data folder)
		string fullWavePaths = ""
		for (i = 0; i < nTraces; i+=1)
			fullWavePaths += GetWavesDataFolder(TraceNameToWaveRef("", StringFromList(i, traceNames)), 2) + ";"
		endfor
		for (i = 0; i < nImages; i+=1)
			fullWavePaths += GetWavesDataFolder(ImageNameToWaveRef("", StringFromList(i, imageNames)), 2) + ";"
		endfor
		
		// find the wave paths that can provide us with a position
		// we assume that the generated image or traces are in 
		// the magic folders.
		string possibleWaves
		variable nPossibleWaves
		if (wantPositions)
			possibleWaves = ListMatch(fullWavePaths, "root:'Localizer Images':*")
		else
			possibleWaves = ListMatch(fullWavePaths, "root:'Tracking Images':*")
		endif
		nPossibleWaves = ItemsInList(possibleWaves)
		// no candidates?
		if (nPossibleWaves == 0)
			break
		endif
		
		// we found at least one possibility. The name of the positions wave
		// is now contained in possibleWaves as :root:'Localizer Images':<name>
		string candidateName = StringFromList(2, StringFromList(0, possibleWaves, ";"), ":")
		// if there is more than one match then require that these are the same
		string alternatePositionsName
		for (i = 1; i < nPossibleWaves; i+=1)	
			alternatePositionsName = StringFromList(2, StringFromList(i, possibleWaves, ";"), ":")
			if (StringMatch(candidateName, alternatePositionsName) != 1)
				break
			endif
		endfor
		
		// the wave name may be quoted, in which case we have to remove them
		candidateName = ReplaceString("'", candidateName, "")
		
		if (wantPositions)
			return GetPositionsWaveReference(candidateName)
		else
			return GetTrackingWaveReference(candidateName)
		endif
	while (0)
	
	// if we're here then we failed to get a match automatically. Prompt to user for a manual selection.
	string selectedName
	string promptString
	if (wantPositions)
		Prompt selectedName, "Positions wave:", popup, GetPossiblePositionsWaves()
	else
		Prompt selectedName,  "Tracks wave:", popup, GetPossibleTrackingWaves()
	endif
	
	DoPrompt "Choose a wave", selectedName
	if (V_flag == 1)
		userCancelled = 1
		return $""
	endif
	
	if (wantPositions)
		return GetPositionsWaveReference(selectedName)
	else
		return GetTrackingWaveReference(selectedName)
	endif
End

Function /WAVE GetPositionsWaveFromGraph(windowName, userCancelled)
	string windowName
	variable &userCancelled
	
	return GetPosOrTrackWaveFromGraph("positions", windowName, userCancelled)
End

Function /WAVE GetTrackingWaveFromGraph(windowName, userCancelled)
	string windowName
	variable &userCancelled
	
	return GetPosOrTrackWaveFromGraph("tracks", windowName, userCancelled)
End

Function /S GetPossiblePositionsWaves()
	string savDF = GetDataFolder(1)
	string listOfWaves = ""
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	NewDataFolder /O root:'Localized Positions'
	
	DFREF packageFolder = root:Packages:Localizer
	DFREF positionsFolder = root:'Localized Positions'
	
	// make a list of all 2D waves in the usual locations
	listOfWaves += RecursiveWaveList(packageFolder, "*", "DIMS:2")
	listOfWaves += RecursiveWaveList(positionsFolder, "*", "DIMS:2")
	
	// now select only those waves that really are positions
	variable nItems = ItemsInList(listOfWaves), nItemsInPath
	variable i
	string currentWaveName, validWaves = ""
	for (i = 0; i < nItems; i+=1)
		currentWaveName = StringFromList(i, listOfWaves)
		if (IsPositionsWave($currentWaveName))
			// this wave looks valid.
			// but reject the wave if it is in the Localizer package folder but not in one of the viewer folders.
			// This way intermediate positions waves created as intermediates in calculations (e.g. drift correction) are discarded.
			if (StringMatch(currentWaveName, "root:Packages:Localizer*") && !StringMatch(currentWaveName, "root:Packages:Localizer:LocalizerViewer*"))
				continue
			endif
			nItemsInPath = ItemsInList(currentWaveName, ":")
			// return only the name, not the complete path
			currentWaveName = StringFromList(nItemsInPath - 1, currentWaveName, ":")
			validWaves += currentWaveName + ";"
		endif
	endfor
	
	if (strlen(validWaves) == 0)
		validWaves = "* No Positions *;"
		return validWaves
	endif
	
	// sort the list alphabetically
	validWaves = SortList(validWaves, ";", 16)
	
	return validWaves
End

Function /S GetPossiblePositionsWavesForPM()
	// get all possible positions waves, but assume that this function
	// will only be called in the context of populating the popupmenu
	// in the analysis controls. Only return those positions that are global
	// (i.e. in the positions DF) or within the DF specific to the calling interface panel
	string windowName = GetNameOfTopInterfacePanel()
	string baseWindowName = GetBaseWindowName(windowName)
	string listOfWaves = ""
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	NewDataFolder /O root:'Localized Positions'
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	DFREF positionsFolder = root:'Localized Positions'
	
	// make a list of all 2D waves in the usual locations
	listOfWaves += RecursiveWaveList(windowDataFolder, "*", "DIMS:2")
	listOfWaves += RecursiveWaveList(positionsFolder, "*", "DIMS:2")
	
	// now select only those waves that really are positions
	// determined by the fact that they have a "LOCALIZATION METHOD" entry in the wave note
	variable nItems = ItemsInList(listOfWaves), nItemsInPath
	variable i
	string currentWaveName, validWaves = ""
	for (i = 0; i < nItems; i+=1)
		currentWaveName = StringFromList(i, listOfWaves)
		if (IsPositionsWave($currentWaveName))
			// this wave looks valid
			// return only the name, not the complete path
			nItemsInPath = ItemsInList(currentWaveName, ":")
			currentWaveName = StringFromList(nItemsInPath - 1, currentWaveName, ":")
			validWaves += currentWaveName + ";"
		endif
	endfor
	
	if (strlen(validWaves) == 0)
		validWaves = "* No Positions *;"
		return validWaves
	endif
	
	// sort the list alphabetically
	validWaves = SortList(validWaves, ";", 16)
	
	return validWaves
End

Function /S GetGaussAstigPositionsWaves()
	string positionsWaves = GetPossiblePositionsWaves()
	
	// check for the "* no positions *" return
	if (StringMatch(positionsWaves, "*No Positions*"))
		positionsWaves = ""
	endif
	variable nPositionsWaves = ItemsInList(positionsWaves)
	
	// only positions fitted using astigmatism are eligible
	string eligiblePositionsNames = ""
	variable i
	for (i = 0; i < nPositionsWaves; i+=1)
		wave thesePositions = GetPositionsWaveReference(StringFromList(i, positionsWaves))
		if (NumberByKey("LOCALIZATION METHOD", note(thesePositions)) == LOCALIZATION_ELLIPSGAUSS_ASTIG)
			eligiblePositionsNames += StringFromList(i, positionsWaves) + ";"
		endif
	endfor
	
	return positionsWaves
End

Function /S GetPossible3DPositionsWaves()
	string positionsWaves = GetPossiblePositionsWaves()
	string possible3DPositionsWaves = ""
	
	if (StringMatch(positionsWaves, "*No Positions*"))
		return ""
	endif
	
	variable nPositionsWaves = ItemsInList(positionsWaves)
	variable i, xCol, yCol, zCol
	for (i = 0; i < nPositionsWaves; i+=1)
		string thesePosName = StringFromList(i, positionsWaves)
		wave thesePos = GetPositionsWaveReference(thesePosName)
		GetColumnsForEmitterPositions(thesePos, xCol, yCol, zCol)
		if (zCol >= 0)
			possible3DPositionsWaves += thesePosName + ";"
		endif
	endfor
	
	return possible3DPositionsWaves
End

Function /S GetSavedPositionsWaves()
	string listOfWaves = ""
	
	NewDataFolder /O root:'Localized Positions'
	
	DFREF positionsFolder = root:'Localized Positions'
	
	listOfWaves += RecursiveWaveList(positionsFolder, "*", "DIMS:2")
	
	// now select only those waves that really are positions
	variable nItems = ItemsInList(listOfWaves), nItemsInPath
	variable i
	string currentWaveName, validWaves = ""
	for (i = 0; i < nItems; i+=1)
		currentWaveName = StringFromList(i, listOfWaves)
		if (IsPositionsWave($currentWaveName))
			// this wave looks valid
			// return only the name, not the complete path
			nItemsInPath = ItemsInList(currentWaveName, ":")
			currentWaveName = StringFromList(nItemsInPath - 1, currentWaveName, ":")
			validWaves += currentWaveName + ";"
		endif
	endfor
	
	if (strlen(validWaves) == 0)
		validWaves = "* No Positions *;"
		return validWaves
	endif
	
	// sort the list alphabetically
	validWaves = SortList(validWaves, ";", 16)
	
	return validWaves
End

Function /S GetPossibleTrackingWaves()
	string listOfWaves = ""
	
	NewDataFolder /O root:'Particle Tracks'
	
	DFREF trackingFolder = root:'Particle Tracks'
	
	// make a list of all 2D waves in the usual locations
	listOfWaves += RecursiveWaveList(trackingFolder, "*", "DIMS:1,BYTE:0,CMPLX:0,DP:0,INTEGER:0,SP:0,TEXT:0")
	
	// now make use of the fact that the tracking waves
	// have a wavenote based on that of the positions
	variable nItems = ItemsInList(listOfWaves), nItemsInPath
	variable i
	string currentWaveName, validWaves = ""
	for (i = 0; i < nItems; i+=1)
		currentWaveName = StringFromList(i, listOfWaves)
		if (IsTrackingWave($currentWaveName))
			// this wave looks valid
			// return only the name, not the complete path
			nItemsInPath = ItemsInList(currentWaveName, ":")
			currentWaveName = StringFromList(nItemsInPath - 1, currentWaveName, ":")
			validWaves += currentWaveName + ";"
		endif
	endfor
	
	if (strlen(validWaves) == 0)
		validWaves = "* No Tracks *;"
		return validWaves
	endif
	
	// sort the list alphabetically
	validWaves = SortList(validWaves, ";", 16)
	
	return validWaves
End

Function /S GetPossibleRegistrationMaps()
	string listOfWaves = ""
	
	NewDataFolder /O root:'Registration Maps'
	
	DFREF registrationFolder = root:'Registration Maps'
	
	listOfWaves += RecursiveWaveList(registrationFolder, "*", "DP:1,CMPLX:0")
	
	// now select only those waves that really are registration maps
	variable nItems = ItemsInList(listOfWaves), nItemsInPath
	variable i
	string currentWaveName, validWaves = ""
	for (i = 0; i < nItems; i+=1)
		currentWaveName = StringFromList(i, listOfWaves)
		if (IsRegistrationMap($currentWaveName))
			// this wave looks valid
			// return only the name, not the complete path
			nItemsInPath = ItemsInList(currentWaveName, ":")
			currentWaveName = StringFromList(nItemsInPath - 1, currentWaveName, ":")
			validWaves += currentWaveName + ";"
		endif
	endfor
	
	if (strlen(validWaves) == 0)
		validWaves = "* No Positions *;"
		return validWaves
	endif
	
	// sort the list alphabetically
	validWaves = SortList(validWaves, ";", 16)
	
	return validWaves
End

Function /S GetPossibleAstig3DCalibrMaps()
	string listOfWaves = ""
	
	NewDataFolder /O root:'Calibration 3D'
	DFREF registrationFolder = root:'Calibration 3D'
	
	listOfWaves += RecursiveWaveList(registrationFolder, "*", "DP:1,CMPLX:0")
	
	// now select only those waves that really are registration maps
	variable nItems = ItemsInList(listOfWaves), nItemsInPath
	variable i
	string currentWaveName, validWaves = ""
	for (i = 0; i < nItems; i+=1)
		currentWaveName = StringFromList(i, listOfWaves)
		if (IsAstigmatism3DCalibrationMap($currentWaveName))
			// this wave looks valid
			// return only the name, not the complete path
			nItemsInPath = ItemsInList(currentWaveName, ":")
			currentWaveName = StringFromList(nItemsInPath - 1, currentWaveName, ":")
			validWaves += currentWaveName + ";"
		endif
	endfor
	
	if (strlen(validWaves) == 0)
		validWaves = "* No Calibrations *;"
		return validWaves
	endif
	
	// sort the list alphabetically
	validWaves = SortList(validWaves, ";", 16)
	
	return validWaves
End

Function /S GetPossibleSOFICombWeights()
	string listOfWaves = ""
	
	NewDataFolder /O root:'SOFI Combination Weights'
	DFREF combinationWeightsFolder = root:'SOFI Combination Weights'
	
	listOfWaves += RecursiveWaveList(combinationWeightsFolder, "*", "")
	
	// now select only those waves that really are registration maps
	variable nItems = ItemsInList(listOfWaves), nItemsInPath
	variable i
	string currentWaveName, validWaves = ""
	for (i = 0; i < nItems; i+=1)
		currentWaveName = StringFromList(i, listOfWaves)
		if (IsSOFICombinationWeightsWave($currentWaveName))
			// this wave looks valid
			// return only the name, not the complete path
			nItemsInPath = ItemsInList(currentWaveName, ":")
			currentWaveName = StringFromList(nItemsInPath - 1, currentWaveName, ":")
			validWaves += currentWaveName + ";"
		endif
	endfor
	
	if (strlen(validWaves) == 0)
		validWaves = "* No Weights *;"
		return validWaves
	endif
	
	// sort the list alphabetically
	validWaves = SortList(validWaves, ";", 16)
	
	return validWaves
End

Function GetOutputStorageTypeAndPath(dataFilePath, outputFormat, outputPath, [offerPresentationOutput])
	string dataFilePath	// used to provide some default output name
	variable &outputFormat
	string &outputPath
	variable offerPresentationOutput
	
	// remember the selected format for future calls
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageDF = root:Packages:Localizer
	variable /G packageDF:V_selectedOutputFormat
	NVAR gSelectedOutputFormat = packageDF:V_selectedOutputFormat
	
	variable havePresentation = 0
	if (!ParamIsDefault(offerPresentationOutput))
		havePresentation = offerPresentationOutput
	endif
	outputPath = ""
	
	string macFilePath = ParseFilePath(5, dataFilePath, ":", 0, 1)
	string baseFileName = ParseFilePath(3, macFilePath, ":", 0 , 1)
	
	// ask how the user wants to save the data
	string allowedFormats = "TIFF file;Deflate-Compressed TIFF File;Multi-file TIFF;Igor wave;"
	if (havePresentation)
		allowedFormats += "Presentation Movie;"
	endif
	variable selectedFormat = gSelectedOutputFormat
	Prompt selectedFormat, "Format: ", popup, allowedFormats
	DoPrompt "Select the desired output format", selectedFormat
	if (V_flag == 1)
		return 0
	endif
	gSelectedOutputFormat = selectedFormat
	
	string saveFileFilter = ""
	string extension	 = ""	// in case we need to manually add it
	string openDialogMessage = ""
	variable refNum
	
	// determine the output format
	string strOutputFormat = StringFromList(selectedFormat - 1, allowedFormats)
	StrSwitch (strOutputFormat)
		case "Igor wave":
			outputFormat = IMAGE_OUTPUT_TYPE_IGOR
			break
		case "TIFF file":
			outputFormat = IMAGE_OUTPUT_TYPE_TIFF
			openDialogMessage = "Specify an output file"
			saveFileFilter = "TIFF File (*.tif):.tif;"
			extension = ".tif"
			break
		case "Multi-file TIFF":
			outputFormat = IMAGE_OUTPUT_TYPE_MULTITIFF
			openDialogMessage = "Specify the output folder and the base file name (e.g. \"base\" to get files like \"base000000.tif\")"
			saveFileFilter = ""
			break
		case "Deflate-Compressed TIFF File":
			outputFormat = IMAGE_OUTPUT_TYPE_COMPR_TIFF
			openDialogMessage = "Specify an output file"
			saveFileFilter = "TIFF File (*.tif):.tif;"
			extension = ".tif"
			break
		case "PDE file":
			outputFormat = IMAGE_OUTPUT_TYPE_PDE
			saveFileFilter = "PDE File (*.pde):.pde;"
			extension = ".pde"
			break
		case "Presentation Movie":
			outputFormat = IMAGE_OUTPUT_TYPE_MOV
			// this needs to be handled by the caller
			outputPath = "DUMMY"	// so that we can distinguish this from cancel
			return 0
			break
		default:
			Abort "Bug: unknown argument to strswitch in GetOutputStorageTypeAndPath"
			break
	EndSwitch
	
	if (outputFormat == IMAGE_OUTPUT_TYPE_IGOR)
		NewDataFolder /O root:'Localizer Movies'
		DFREF movieFolder = root:'Localizer Movies'
		outputPath = "root:'Localizer Movies':" + GetNewOutputWaveName(CleanupName(baseFileName, 0), movieFolder, "")
		if (strlen(outputPath) == 0)
			return 0
		endif
	else
		Open /D /F=saveFileFilter /M=openDialogMessage refNum
		if (strlen(S_fileName) == 0)
			outputPath = ""
			return 0
		endif
		outputPath = S_fileName
		
		// did we get an extension and/or do we need one?
		string leafName = ParseFilePath(0, outputPath, ":", 1, 0)
		variable haveExtension = (StringMatch(leafName, "*.*") != 0)
		variable needExtension = (outputFormat != IMAGE_OUTPUT_TYPE_MULTITIFF)
		if (!haveExtension && needExtension)
			outputPath += extension
		elseif (haveExtension && !needExtension)
			outputPath = outputPath[0, StrSearch(outputPath, ".", strlen(outputPath) - 1, 1) - 1]	// drops the extension
		endif
	endif
End

Function DoClusteringAnalysis(doLFunction, isBivariate)
	variable doLFunction
	variable isBivariate
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	
	Variable /G root:Packages:Localizer:V_lClusteringRange
	Variable /G root:Packages:Localizer:V_lClusteringNBins
	Variable /G root:Packages:Localizer:V_lCluseringNPointsToSample
	Variable /G root:Packages:Localizer:V_lClusteringNMonteCarloSamples
	NVAR calculationRange = root:Packages:Localizer:V_lClusteringRange
	NVAR nBins = root:Packages:Localizer:V_lClusteringNBins
	NVAR nPointsToSample = root:Packages:Localizer:V_lCluseringNPointsToSample
	NVAR nMonteCarloSamples= root:Packages:Localizer:V_lClusteringNMonteCarloSamples
	
	if (calculationRange <= 0)
		calculationRange = 30
	endif
	if (nBins <= 0)
		nBins = 100
	endif
	if (nPointsToSample <= 0)
		nPointsToSample = 50000
	endif
	if (nMonteCarloSamples == 0)
		nMonteCarloSamples = 9
	endif
	
	string posName
	string posName2
	variable localCalculationRange = calculationRange
	variable localNBins = nBins
	variable localNPointsToSample = nPointsToSample
	variable localNMonteCarloSamples = nMonteCarloSamples
	
	Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
	if (isBivariate)
		Prompt posName2, "Positions wave:", popup, GetPossiblePositionsWaves()
	endif
	Prompt localCalculationRange, "Calculation range (pixels):"
	Prompt localNBins, "Number of bins:"
	Prompt localNPointsToSample, "Number of points to randomly sample (-1 for all points):"
	Prompt localNMonteCarloSamples, "Monte Carlo trials for confidence analysis:"
	
	if (!isBivariate)
		DoPrompt "Calculation parameters", posName, localCalculationRange, localNBins, localNPointsToSample, localNMonteCarloSamples
	else
		DoPrompt "Calculation parameters", posName, posName2, localCalculationRange, localNBins, localNPointsToSample, localNMonteCarloSamples
	endif
	if (V_flag == 1)	// cancel
		return 0
	endif
	
	calculationRange = localCalculationRange
	nBins = localNBins
	nPointsToSample = localNPointsToSample
	nMonteCarloSamples = localNMonteCarloSamples
	
	// get positions waves
	wave /Z pos = GetPositionsWaveReference(posName)
	if (WaveExists(pos) == 0)
		Abort "The selected wave doesn't seem to exist!"
	endif
	if (isBivariate)
		wave /Z pos2 = GetPositionsWaveReference(posName2)
		if (WaveExists(pos2) == 0)
			Abort "The selected wave doesn't seem to exist!"
		endif
	else
		wave /Z pos2 = $""
	endif
	
	// get active region
	variable xStart, xEnd, yStart, yEnd
	GetMinMaxCoordsFromPositions(pos, xStart, xEnd, yStart, yEnd)
	if (isBivariate)
		variable xStart2, xEnd2, yStart2, yEnd2
		GetMinMaxCoordsFromPositions(pos2, xStart2, xEnd2, yStart2, yEnd2)
		xStart = min(xStart, xStart2)
		xEnd = max(xEnd, xEnd2)
		yStart = min(yStart, yStart2)
		yEnd = max(yEnd, yEnd2)
	endif
	
	// do the actual calculation
	CalculateAndDisplayClustering(doLFunction, pos, calculationRange, nBins, nPointsToSample, nMonteCarloSamples, xStart, xEnd, yStart, yEnd, pos2=pos2)
End

Function CalculateAndDisplayClustering(doLFunction, pos, calculationRange, nBins, nPointsToSample, nMonteCarloSamples, xStart, xEnd, yStart, yEnd, [pos2])
	variable doLFunction
	wave pos
	variable calculationRange, nBins, nPointsToSample, nMonteCarloSamples
	variable xStart, xEnd, yStart, yEnd
	wave /Z pos2
	
	DFREF savDF = GetDataFolderDFR()
	variable isBivariate = (!ParamIsDefault(pos2) && WaveExists(pos2))
	
	// try to get the CCD pixel size
	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(pos, xSize, ySize, pixelSize)
	xSize = xEnd - xStart
	ySize = yEnd - yStart
	
	variable nPositions1 = DimSize(pos, 0)
	if (isBivariate)
		variable nPositions2 = DimSize(pos2, 0)
	endif
	
	if (nPointsToSample <= 0)
		nPointsToSample = nPositions1
	endif
	SetDataFolder root:Packages:Localizer
	if (nPointsToSample < nPositions1)
		// calculate only some points
		StatsSample /N=(nPointsToSample) /MC pos
		wave M_Sampled
		Duplicate /O M_Sampled, positionsToAnalyze
		Note /K positionsToAnalyze, note(pos)
		if (isBivariate)
			StatsSample /N=(min(nPointsToSample * nPositions2 / nPositions1, nPositions2)) /MC pos2
			wave M_Sampled
			Duplicate /O M_Sampled, positionsToAnalyze2
			Note /K positionsToAnalyze2, note(pos2)
		endif
	else
		// calculate over all points
		wave positionsToAnalyze = pos
		if (isBivariate)
			wave positionsToAnalyze2 = pos2
		endif
	endif
	
	if (!isBivariate)
		if (doLFunction)
			RipleyLFunctionClustering /RNGE={calculationRange, nBins} /REGN={xStart, xEnd, yStart, yEnd} positionsToAnalyze
		else
			PairwiseCorrelationClustering /RNGE={calculationRange, nBins} /REGN={xStart, xEnd, yStart, yEnd} positionsToAnalyze
		endif
	else
		if (doLFunction)
			RipleyLFunctionClustering /RNGE={calculationRange, nBins} /REGN={xStart, xEnd, yStart, yEnd} positionsToAnalyze, positionsToAnalyze2
		else
			PairwiseCorrelationClustering /RNGE={calculationRange, nBins} /REGN={xStart, xEnd, yStart, yEnd} positionsToAnalyze, positionsToAnalyze2
		endif
	endif
	SetDataFolder savDF
	
	string windowName
	if (doLFunction)
		wave W_LFunction = root:Packages:Localizer:W_lFunction
		Duplicate /O W_LFunction, root:Packages:Localizer:W_LFunctionSubX
		wave W_LFunctionSubX = root:Packages:Localizer:W_LFunctionSubX
		
		W_LFunctionSubX = W_LFunction - x
		wave W_Result = root:Packages:Localizer:W_LFunctionSubX
		windowName = "LFunction"
	else
		wave W_Result = root:Packages:Localizer:W_PairwiseCorrelation
		windowName = "PairwiseCorrelation"
	endif
	
	// adjust the size to nm if possible
	if (pixelSize != 0)
		SetScale /I x, 0, calculationRange * pixelSize, W_Result
	endif
	
	DoWindow /F $windowName
	if (V_flag != 1)
		Display /K=1/N=$windowName W_Result
		if (doLFunction)
			Label /W=$windowName Left, "L Function - x"
		else
			Label /W=$windowName Left, "Pairwise correlation"
		endif
	endif
	if (pixelSize == 0)
		Label /W=$windowName Bottom, "Distance (pixel)"
	else
		Label /W=$windowName Bottom, "Distance (nm)"
	endif
	
	// also simulate confidence bands
	string lowerBoundWaveName, upperBoundWaveName
	if (doLFunction)
		lowerBoundWaveName = "W_LFunctionLowerBound"
		upperBoundWaveName = "W_LFunctionUpperBound"
	else
		lowerBoundWaveName = "W_PairwiseCorrelationLowerBound"
		upperBoundWaveName = "W_PairwiseCorrelationUpperBound"
	endif
	if (nMonteCarloSamples > 0)
		variable generateRandomPositions = !isBivariate
		wave /WAVE W_SimulationResult = SimulateRandomClusteringBounds(doLFunction, generateRandomPositions, pos, nMonteCarloSamples, nPointsToSample, calculationRange, nBins, xSize, ySize, positions2=pos2)
		wave lowerBound = W_SimulationResult[0]
		wave upperBound = W_SimulationResult[1]
		if (pixelSize != 0)
			SetScale /I x, 0, calculationRange * pixelSize, lowerBound, upperBound
		endif
		Duplicate /O lowerBound, root:Packages:Localizer:$lowerBoundWaveName
		Duplicate /O upperBound, root:Packages:Localizer:$upperBoundWaveName
		wave savedLowerBound = root:Packages:Localizer:$lowerBoundWaveName
		wave savedUpperBound = root:Packages:Localizer:$upperBoundWaveName
		string traceList = TraceNameList(windowName, ";", 1)
		if (FindListItem(lowerBoundWaveName, traceList) == -1)
			AppendToGraph /W=$windowName savedLowerBound
			ModifyGraph /W=$windowName lstyle($lowerBoundWaveName)=3
		endif
		if (FindListItem(upperBoundWaveName, traceList) == -1)
			AppendToGraph /W=$windowName savedUpperBound
			ModifyGraph /W=$windowName lstyle($upperBoundWaveName)=3
		endif
	else
		// remove the confidence bounds if they are being shown
		RemoveFromGraph /W=$windowName /Z $lowerBoundWaveName
		RemoveFromGraph /W=$windowName /Z $upperBoundWaveName
	endif
End

Function CalculateMSDs_menu([tracksWave, providedTracksWaveName])
	wave /WAVE tracksWave			// if specified then do not prompt the user for a tracks wave but use this instead
	string providedTracksWaveName	// 'real' name of the provided tracks wave. Needed to avoid "_free_" being selected as the name of the tracks wave.
	
	if (!ParamIsDefault(tracksWave) && ParamIsDefault(providedTracksWaveName))
		Abort "If tracksWave is specified then providedTracksWaveName must also be specified in CalculateMSDs_menu()"
	endif
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	NVAR /Z gMaxTimeLag = packageFolder:V_MSDMaxTimeLag
	if (!NVAR_Exists(gMaxTimeLag))
		variable /G packageFolder:V_MSDMaxTimeLag
		NVAR gMaxTimeLag = packageFolder:V_MSDMaxTimeLag
		gMaxTimeLag = 10
	endif
	NVAR /Z gMinTrackLength = packageFolder:V_MSDMinTrackLength
	if (!NVAR_Exists(gMinTrackLength))
		variable /G packageFolder:V_MSDMinTrackLength
		NVAR gMinTrackLength = packageFolder:V_MSDMinTrackLength
		gMinTrackLength = -1
	endif
	NVAR /Z gMaxTrackLength = packageFolder:V_MSDMaxTrackLength
	if (!NVAR_Exists(gMaxTrackLength))
		variable /G packageFolder:V_MSDMaxTrackLength
		NVAR gMaxTrackLength = packageFolder:V_MSDMaxTrackLength
		gMaxTrackLength = -1
	endif
	
	string tracksWaveName
	variable maxTimeLag = gMaxTimeLag, minTrackLength = gMinTrackLength, maxTrackLength = gMaxTrackLength
	
	Prompt tracksWaveName, "Tracking wave:", popup, GetPossibleTrackingWaves()
	Prompt maxTimeLag, "Longest time lag (frames):"
	Prompt minTrackLength, "Min track length:"
	Prompt maxTrackLength, "Max track length:"
	
	if (ParamIsDefault(tracksWave))
		DoPrompt "Calculation parameters", tracksWaveName, maxTimeLag, minTrackLength, maxTrackLength
	else
		DoPrompt "Calculation parameters", maxTimeLag, minTrackLength, maxTrackLength
	endif
	if (V_flag == 1)	// cancel
		return 0
	endif
	if (maxTimeLag <= 0)
		Abort "Invalid max time lag"
	endif
	
	gMaxTimeLag = maxTimeLag
	gMinTrackLength = minTrackLength
	gMaxTrackLength = maxTrackLength
	
	if (ParamIsDefault(tracksWave))
		wave /wave tracksWave = GetTrackingWaveReference(tracksWaveName)
	else
		tracksWaveName = providedTracksWaveName
	endif
	if (!WaveExists(tracksWave))
		Abort "Unable to find the requested track wave"
	endif
	
	variable nTracks = DimSize(tracksWave, 0)
	NewDataFolder /O root:'Tracking Results'
	NewDataFolder /O root:'Tracking Results':$tracksWaveName
	DFREF outputFolder = root:'Tracking Results':$tracksWaveName
	
	wave M_MSD = CalculateMSDs(tracksWave, maxTimeLag, minTrackLength, maxTrackLength)
	variable nLags = DimSize(M_MSD, 0)
	Make /O /N=(nLags) /D outputFolder:W_MSD, outputFolder:W_MSDError, outputFolder:W_MSDNTrackContributions
	wave W_MSD = outputFolder:W_MSD
	wave W_MSDError = outputFolder:W_MSDError
	wave W_MSDNTrackContributions = outputFolder:W_MSDNTrackContributions
	SetScale /P x, DimDelta(M_MSD, 0), DimOffset(M_MSD, 0), W_MSD, W_MSDError, W_MSDNTrackContributions
	W_MSD = M_MSD[p][0]
	W_MSDError = M_MSD[p][1]
	W_MSDNTrackContributions = M_MSD[p][2]
	
	string windowName = CleanupName("MSD_" + tracksWaveName, 0)
	DoWindow /F $windowName
	if (!V_flag)
		Display /K=1/N=$windowName W_MSD
		ErrorBars /W=$windowName W_MSD Y,wave=(W_MSDError,W_MSDError)
		AppendToGraph /R /W=$windowName W_MSDNTrackContributions
		ModifyGraph /W=$windowName rgb(W_MSDNTrackContributions)=(0,0,65535)
	endif
	DoWindow /T $windowName, "MSDs from " + tracksWaveName + " (using " + num2str(nTracks) + " tracks)"
	
	Label right, "Number of contributing tracks"
	// if we have information on the pixel size then use it
	variable pixelSize = NumberByKey("X PIXEL SIZE", note(tracksWave))
	if (NumType(pixelSize) == 0)
		W_MSD *= pixelSize^2 / 1e6		// converts nm^2 to µm^2
		W_MSDError *= pixelSize^2 / 1e6
		Label left, "MSD (µm\\S2\\M)"
	else
		Label left, "MSD (pixel\\S2\\M)"
	endif
	
	// if we have information on the time delta then use it
	variable timeDelta = NumberByKey("TRACKING TIME DELTA", note(tracksWave))
	if ((NumType(timeDelta) == 0) && (timeDelta > 0))
		SetScale /P x, DimOffset(W_MSD, 0) * timeDelta, DimDelta(W_MSD, 0) * timeDelta, W_MSD, W_MSDError
		Label bottom, "Time lag (s)"
	else
		Label bottom, "Time lag (frames)"
	endif
	
End

Function CalculateCDFs_menu([tracksWave, providedTracksWaveName])
	wave /wave tracksWave				// if specified then use these tracks instead of prompting the user
	string providedTracksWaveName	// 'real' name of the provided tracks wave. Needed to avoid "_free_" being selected as the name of the tracks wave.
	
	if (!ParamIsDefault(tracksWave) && ParamIsDefault(providedTracksWaveName))
		Abort "If tracksWave is specified then providedTracksWaveName must also be specified in CalculateMSDs_menu()"
	endif
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	NVAR /Z gMaxTimeLag = packageFolder:V_CDFMaxTimeLag
	if (!NVAR_Exists(gMaxTimeLag))
		variable /G packageFolder:V_CDFMaxTimeLag
		NVAR gMaxTimeLag = packageFolder:V_CDFMaxTimeLag
		gMaxTimeLag = 20
	endif
	NVAR /Z gMinTrackLength = packageFolder:V_CDFMinTrackLength
	if (!NVAR_Exists(gMinTrackLength))
		variable /G packageFolder:V_CDFMinTrackLength
		NVAR gMinTrackLength = packageFolder:V_CDFMinTrackLength
		gMinTrackLength = -1
	endif
	NVAR /Z gMaxTrackLength = packageFolder:V_CDFMaxTrackLength
	if (!NVAR_Exists(gMaxTrackLength))
		variable /G packageFolder:V_CDFMaxTrackLength
		NVAR gMaxTrackLength = packageFolder:V_CDFMaxTrackLength
		gMaxTrackLength = -1
	endif
	NVAR /Z gMaxSqDisplacement = packageFolder:V_CDFMaxSqDisplacement
	if (!NVAR_Exists(gMaxSqDisplacement))
		variable /G packageFolder:V_CDFMaxSqDisplacement
		NVAR gMaxSqDisplacement = packageFolder:V_CDFMaxSqDisplacement
		gMaxSqDisplacement = 50
	endif
	NVAR /Z gNPointsInCDF = packageFolder:V_CDFNPoints
	if (!NVAR_Exists(gNPointsInCDF))
		variable /G packageFolder:V_CDFNPoints
		NVAR gNPointsInCDF = packageFolder:V_CDFNPoints
		gNPointsInCDF = 20
	endif
	
	string tracksWaveName
	variable maxTimeLag = gMaxTimeLag, minTrackLength = gMinTrackLength, maxTrackLength = gMaxTrackLength
	variable maxSqDisplacement = gMaxSqDisplacement, nPointsInCDF = gNPointsInCDF
	
	Prompt tracksWaveName, "Tracking wave:", popup, GetPossibleTrackingWaves()
	Prompt maxTimeLag, "Longest time lag (frames):"
	Prompt nPointsInCDF, "N points in CDF:"
	Prompt maxSqDisplacement, "Largest displacement (pixels^2):"
	Prompt minTrackLength, "Min track length:"
	Prompt maxTrackLength, "Max track length:"
	
	if (ParamIsDefault(tracksWave))
		DoPrompt "Calculation parameters", tracksWaveName, maxTimeLag, maxSqDisplacement, nPointsInCDF, minTrackLength, maxTrackLength
	else
		DoPrompt "Calculation parameters", maxTimeLag, maxSqDisplacement, nPointsInCDF, minTrackLength, maxTrackLength
	endif
	if (V_flag == 1)	// cancel
		return 0
	endif
	if (maxTimeLag <= 0)
		Abort "Invalid max time lag"
	endif
	if ((nPointsInCDF <= 0) || (maxSqDisplacement <= 0))
		Abort "Largest displacement and number of points in the CDF must be positive"
	endif
	
	gMaxTimeLag = maxTimeLag
	gMinTrackLength = minTrackLength
	gMaxTrackLength = maxTrackLength
	gMaxSqDisplacement = maxSqDisplacement
	gNPointsInCDF = nPointsInCDF
	
	if (ParamIsDefault(tracksWave))
		wave /wave tracksWave = GetTrackingWaveReference(tracksWaveName)
	else
		tracksWaveName = providedTracksWaveName
	endif
	if (!WaveExists(tracksWave))
		Abort "Unable to find the requested track wave"
	endif
	
	variable nTracks = DimSize(tracksWave, 0)
	NewDataFolder /O root:'Tracking Results'
	NewDataFolder /O root:'Tracking Results':$tracksWaveName
	DFREF outputFolder = root:'Tracking Results':$tracksWaveName
	
	wave M_CDF = CalculateCDFs(tracksWave, maxTimeLag, maxSqDisplacement, nPointsInCDF, minTrackLength, maxTrackLength)
	Duplicate /O M_CDF, outputFolder:M_CDF
	wave M_CDF = outputFolder:M_CDF
	
	string windowName = CleanupName("CDF_" + tracksWaveName, 0)
	DoWindow /F $windowName
	if (!V_flag)
		Display /N=$windowName /K=1
	endif
	DoWindow /T $windowName, "CDFs from " + tracksWaveName + " (using " + num2str(nTracks) + " tracks)"

	string tracesInGraph = TraceNameList(windowName, ";", 1), instanceName
	variable nTracesInGraph = ItemsInList(tracesInGraph)
	variable maxInstanceOnGraph = 0
	variable i
	for (i = 0; i < nTracesInGraph; i+=1)
		instanceName = StringFromList(i, tracesInGraph)
		if (StringMatch(instanceName, "M_CDF#*"))
			maxInstanceOnGraph = max(maxInstanceOnGraph, ExtractInstanceNumber(instanceName))
		endif
	endfor
	for (i = 0; i < DimSize(M_CDF, 0); i+=1)
		if (i == 0)
			instanceName = "M_CDF"
		else
			instanceName = "M_CDF#" + num2istr(i)
		endif
		if (FindListItem(instanceName, tracesInGraph) == -1)
			AppendToGraph /W=$windowName M_CDF[i][]
		endif
	endfor
	// remove all unneeded instance names
	for (i = maxInstanceOnGraph; i >= DimSize(M_CDF, 0); i-=1)
		instanceName = "M_CDF#" + num2istr(i)
		RemoveFromGraph /W=$windowName $instanceName
	endfor
	
	// if we have information on the pixel size then use it
	variable pixelSize = NumberByKey("X PIXEL SIZE", note(tracksWave))
	if (NumType(pixelSize) == 0)
		SetScale /P y, DimOffset(M_CDF, 0) * pixelSize^2 / 1e6, DimDelta(M_CDF, 0) * pixelSize^2 / 1e6, M_CDF
		Label /W=$windowName bottom, "r\\S2\\M (µm\\S2\\M)"
	else
		Label /W=$windowName bottom, "r\\S2\\M (pixels\\S2\\M)"
	endif
	Label /W=$windowName left, "CDF"
	
	// if we have information on the time delta then use it
	variable timeDelta = NumberByKey("TRACKING TIME DELTA", note(tracksWave))
	if ((NumType(timeDelta) == 0) && (timeDelta > 0))
		SetScale /P x, DimOffset(M_CDF, 0) * timeDelta, DimDelta(M_CDF, 0) * timeDelta, M_CDF
	endif
	
	string legendString = ""
	for (i = 0; i < DimSize(M_CDF, 0); i+=1)
		if (i == 0)
			if ((NumType(timeDelta) == 0) && (timeDelta > 0))
				legendString += "(lag in seconds)\r"
			else
				legendString += "(lag in frames)\r"
			endif
			legendString += "\\s(M_CDF) lag = " + num2str(DimOffset(M_CDF, 0) + i * DimDelta(M_CDF, 0)) + "\r"
		else
			legendString += "\\s(M_CDF#" + num2istr(i) + ") lag = " + num2str(DimOffset(M_CDF, 0) + i * DimDelta(M_CDF, 0)) + "\r"
		endif
	endfor
	Legend /W=$windowName/C/N=text0/J  legendString
End

Function DisplacementHistograms_menu([tracksWave, providedTracksWaveName])
	wave /wave tracksWave				// if specified then use these tracks instead of prompting the user
	string providedTracksWaveName	// 'real' name of the provided tracks wave. Needed to avoid "_free_" being selected as the name of the tracks wave.
	
	if (!ParamIsDefault(tracksWave) && ParamIsDefault(providedTracksWaveName))
		Abort "If tracksWave is specified then providedTracksWaveName must also be specified in DisplacementHistograms_menu()"
	endif
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	NVAR /Z gMaxTimeLag = packageFolder:V_DispHistsMaxTimeLag
	if (!NVAR_Exists(gMaxTimeLag))
		variable /G packageFolder:V_DispHistsMaxTimeLag
		NVAR gMaxTimeLag = packageFolder:V_DispHistsMaxTimeLag
		gMaxTimeLag = 5
	endif
	NVAR /Z gMinTrackLength = packageFolder:V_DispHistsMinTrackLength
	if (!NVAR_Exists(gMinTrackLength))
		variable /G packageFolder:V_DispHistsMinTrackLength
		NVAR gMinTrackLength = packageFolder:V_DispHistsMinTrackLength
		gMinTrackLength = -1
	endif
	NVAR /Z gMaxTrackLength = packageFolder:V_DispHistsMaxTrackLength
	if (!NVAR_Exists(gMaxTrackLength))
		variable /G packageFolder:V_DispHistsMaxTrackLength
		NVAR gMaxTrackLength = packageFolder:V_DispHistsMaxTrackLength
		gMaxTrackLength = -1
	endif
	NVAR /Z gMaxDisplacement = packageFolder:V_DispHistsMaxDisplacement
	if (!NVAR_Exists(gMaxDisplacement))
		variable /G packageFolder:V_DispHistsMaxDisplacement
		NVAR gMaxDisplacement = packageFolder:V_DispHistsMaxDisplacement
		gMaxDisplacement = 25
	endif
	NVAR /Z gNHistogramBins = packageFolder:V_DispHistsNPoints
	if (!NVAR_Exists(gNHistogramBins))
		variable /G packageFolder:V_DispHistsNPoints
		NVAR gNHistogramBins = packageFolder:V_DispHistsNPoints
		gNHistogramBins = 50
	endif
	
	string tracksWaveName
	variable maxTimeLag = gMaxTimeLag, minTrackLength = gMinTrackLength, maxTrackLength = gMaxTrackLength
	variable maxDisplacement = gMaxDisplacement, nHistogramBins = gNHistogramBins
	
	Prompt tracksWaveName, "Tracking wave:", popup, GetPossibleTrackingWaves()
	Prompt maxTimeLag, "Longest time lag (frames):"
	Prompt nHistogramBins, "N bins in histogram:"
	Prompt maxDisplacement, "Largest displacement (pixels):"
	Prompt minTrackLength, "Min track length:"
	Prompt maxTrackLength, "Max track length:"
	
	if (ParamIsDefault(tracksWave))
		DoPrompt "Calculation parameters", tracksWaveName, maxTimeLag, maxDisplacement, nHistogramBins, minTrackLength, maxTrackLength
	else
		DoPrompt "Calculation parameters", maxTimeLag, maxDisplacement, nHistogramBins, minTrackLength, maxTrackLength
	endif
	if (V_flag == 1)	// cancel
		return 0
	endif
	if (maxTimeLag <= 0)
		Abort "Invalid max time lag"
	endif
	if ((nHistogramBins <= 1) || (maxDisplacement <= 0))
		Abort "Largest displacement and number of points in the CDF must be positive"
	endif
	
	gMaxTimeLag = maxTimeLag
	gMinTrackLength = minTrackLength
	gMaxTrackLength = maxTrackLength
	gMaxDisplacement = maxDisplacement
	gNHistogramBins = nHistogramBins
	
	if (ParamIsDefault(tracksWave))
		wave /wave tracksWave = GetTrackingWaveReference(tracksWaveName)
	else
		tracksWaveName = providedTracksWaveName
	endif
	if (!WaveExists(tracksWave))
		Abort "Unable to find the requested track wave"
	endif
	
	variable nTracks = DimSize(tracksWave, 0)
	NewDataFolder /O root:'Tracking Results'
	NewDataFolder /O root:'Tracking Results':$tracksWaveName
	DFREF outputFolder = root:'Tracking Results':$tracksWaveName
	
	wave M_DisplacementHistograms = DisplacementHistograms(tracksWave, maxDisplacement, nHistogramBins, maxTimeLag, minTrackLength, maxTrackLength)
	Duplicate /O M_DisplacementHistograms, outputFolder:M_DisplacementHistograms
	wave M_DisplacementHistograms = outputFolder:M_DisplacementHistograms
	
	string windowName = CleanupName("Hists_" + tracksWaveName, 0)
	DoWindow /F $windowName
	if (!V_flag)
		Display /N=$windowName /K=1
	endif
	DoWindow /T $windowName, "Displacement Histograms from " + tracksWaveName + " (using " + num2str(nTracks) + " tracks)"

	string tracesInGraph = TraceNameList(windowName, ";", 1), instanceName
	variable nTracesInGraph = ItemsInList(tracesInGraph)
	variable maxInstanceOnGraph = 0
	variable i
	for (i = 0; i < nTracesInGraph; i+=1)
		instanceName = StringFromList(i, tracesInGraph)
		if (StringMatch(instanceName, "M_DisplacementHistograms#*"))
			maxInstanceOnGraph = max(maxInstanceOnGraph, ExtractInstanceNumber(instanceName))
		endif
	endfor
	for (i = 0; i < DimSize(M_DisplacementHistograms, 0); i+=1)
		if (i == 0)
			instanceName = "M_DisplacementHistograms"
		else
			instanceName = "M_DisplacementHistograms#" + num2istr(i)
		endif
		if (FindListItem(instanceName, tracesInGraph) == -1)
			AppendToGraph /W=$windowName M_DisplacementHistograms[i][]
		endif
	endfor
	// remove all unneeded instance names
	for (i = maxInstanceOnGraph; i >= DimSize(M_DisplacementHistograms, 0); i-=1)
		instanceName = "M_DisplacementHistograms#" + num2istr(i)
		RemoveFromGraph /W=$windowName $instanceName
	endfor
	
	// if we have information on the pixel size then use it
	variable pixelSize = NumberByKey("X PIXEL SIZE", note(tracksWave))
	if (NumType(pixelSize) == 0)
		SetScale /P y, DimOffset(M_DisplacementHistograms, 0) * pixelSize / 1e3, DimDelta(M_DisplacementHistograms, 0) * pixelSize / 1e3, M_DisplacementHistograms
		Label /W=$windowName bottom, "r (µm)"
	else
		Label /W=$windowName bottom, "r (pixels)"
	endif
	Label /W=$windowName left, "Occurrence"
	
	// if we have information on the time delta then use it
	variable timeDelta = NumberByKey("TRACKING TIME DELTA", note(tracksWave))
	if ((NumType(timeDelta) == 0) && (timeDelta > 0))
		SetScale /P x, DimOffset(M_DisplacementHistograms, 0) * timeDelta, DimDelta(M_DisplacementHistograms, 0) * timeDelta, M_DisplacementHistograms
	endif
	
	string legendString = ""
	for (i = 0; i < DimSize(M_DisplacementHistograms, 0); i+=1)
		if (i == 0)
			if ((NumType(timeDelta) == 0) && (timeDelta > 0))
				legendString += "(lag in seconds)\r"
			else
				legendString += "(lag in frames)\r"
			endif
			legendString += "\\s(M_DisplacementHistograms) lag = " + num2str(DimOffset(M_DisplacementHistograms, 0) + i * DimDelta(M_DisplacementHistograms, 0)) + "\r"
		else
			legendString += "\\s(M_DisplacementHistograms#" + num2istr(i) + ") lag = " + num2str(DimOffset(M_DisplacementHistograms, 0) + i * DimDelta(M_DisplacementHistograms, 0)) + "\r"
		endif
	endfor
	Legend /W=$windowName/C/N=text0/J  legendString
End


























// D. Dibble created functions----------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
Function ExportScatterDataMenu()

	

	variable /G FilteringType
	// create a window to check which kind of initial filtering to do first

	NewPanel /W=(536,238,686,346)/N=FilterPanel
	ShowTools/A
	CheckBox NoFilter,pos={25.00,10.00},size={75.00,16.00}, value= 1, proc=NoFilter, title="No Filtering"
	CheckBox DuplicateSearch,pos={25.00,35.00},size={98.00,16.00}, value=0,  proc=DupFilter, title="Filter Duplicates"	
	CheckBox FilterTracks,pos={25.00,60.00},size={103.00,16.00}, value=0,  proc=TrackFilter, title="Track Points Only"
	Button Execute,pos={50.00,80.00},size={50.00,20.00},proc=ExecuteButton, title="Execute"




	


End


Function ExecuteButton(CtrlName) : ButtonControl
	String CtrlName

	svar imageFileLocation
	string posName
	nvar FilteringType = root:FilteringType


	KillWindow/Z FilterPanel

	if (FilteringType != 2)
		// get the localized points (not tracks)
		Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
		DoPrompt "Select the positions wave", posName
		if (V_flag == 1)	// cancel
			return 0
		endif
		// get the position wave
		wave /Z pos = GetPositionsWaveReference(posName)
		if (WaveExists(pos) == 0)
			Abort "The selected wave doesn't seem to exist!"
		endif
	else
		// get the tracks we want to save
		String tracksName
		Prompt tracksName, "Tracks wave: ", popup, GetPossibleTrackingWaves()
		DoPrompt "Select the tracks to be saved", tracksName
		if (V_flag == 1)
			return 0
		endif
		// in case tracking was done, get the info to use as a filter for the data:
		wave /Z W_Tracks = GetTrackingWaveReference(tracksName)
		if (WaveExists(W_Tracks) == 0)
			Abort "The selected wave doesn't seem to exist!"
		endif
		// also need the point data to get the dimensions
		// get the localized points (not tracks)
		Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
		DoPrompt "Select the positions wave", posName
		if (V_flag == 1)	// cancel
			return 0
		endif
		// get the position wave
		wave /Z pos = GetPositionsWaveReference(posName)
		if (WaveExists(pos) == 0)
			Abort "The selected wave doesn't seem to exist!"
		endif
	endif
				



	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(pos, xSize, ySize, pixelSize)
	

	// the wave we feed the next function depends on the filter type
	// no tracks
	if (FilteringType != 2)
		ExportScatterData(pos, xSize, ySize, pixelSize, FilteringType, imageFileLocation)
	else
	//tracks
		ExportScatterData(W_Tracks, xSize, ySize, pixelSize, FilteringType, imageFileLocation)
	endif	
End


Function NoFilter(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	nvar FilteringType = root:FilteringType

	switch( cba.eventCode )
		case 2: // mouse up
						
			FilteringType = 0
			CheckBox NoFilter, value = 1
			CheckBox DuplicateSearch, value = 0
			CheckBox FilterTracks, value = 0

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function DupFilter(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	nvar FilteringType = root:FilteringType

	switch( cba.eventCode )
		case 2: // mouse up
						
			FilteringType = 1
			CheckBox NoFilter, value = 0
			CheckBox DuplicateSearch, value = 1
			CheckBox FilterTracks, value = 0

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function TrackFilter(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	nvar FilteringType = root:FilteringType

	switch( cba.eventCode )
		case 2: // mouse up
						
			FilteringType = 2
			CheckBox NoFilter, value = 0
			CheckBox DuplicateSearch, value = 0
			CheckBox FilterTracks, value = 1

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End



// Function to toggle menu item and act on it if exported points are available the above function 
Function SubtractBackground()
	String PointsExported = StrVarOrDefault("root:ExportedStr","(Subtract Background")
	if(CmpStr(PointsExported, "(Subtract Background")==0)
		String/G root:ExportedStr = "(Subtract Background"
		//do nothing, there is no data to work with
	else
		// execute subtraction routine 
		SubtractBackgroundPanel()

	endif	
		
End


Function Decimalize()
	
	string posName
	Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
	DoPrompt "Select the positions wave", posName
	if (V_flag == 1)	// cancel
		return 0
	endif
	wave /Z pos = GetPositionsWaveReference(posName)
	if (WaveExists(pos) == 0)
		Abort "The selected wave doesn't seem to exist!"
	endif
	DecimalizePanel(pos)
		
End


Function CompareBackgroundChannelImport()
	CompareBackgroundChannel()
End


Function ExtractTrackHistogram()

	DFREF packagekageFolder = root:Packages:Localizer

	string trackingWaveList = GetPossibleTrackingWaves()
	if (ItemsInList(trackingWaveList) == 0)
		Abort "No particle tracks seem to have been calculated"
	endif
	
	string trackingWaveName
	Prompt trackingWaveName, "Track wave:", popup, trackingWaveList
	DoPrompt "Select the tracks to filter", trackingWaveName
	if (V_flag != 0)
		return 0
	endif
	
	wave /WAVE /Z trackWave = GetTrackingWaveReference(trackingWaveName)
	if (WaveExists(trackWave) == 0)
		Abort "Unable to find the requested track"
	endif



	Wave HistogramPlot = GetTrackLengthHistogram(trackWave)

	Duplicate /O HistogramPlot LocalPlot

	Make /O/T/N=(numpnts(LocalPlot)) HistogramXAxis

	variable itor
	for(itor = 0; itor < numpnts(LocalPlot); itor++)
		HistogramXAxis[itor] = num2str(itor)		
	endfor

	Display /W=(20,20,520,430)/N=Histogram LocalPlot vs HistogramXAxis

End


Function ExtractVelocityHistograms()
	DFREF packagekageFolder = root:Packages:Localizer

	string trackingWaveList = GetPossibleTrackingWaves()
	if (ItemsInList(trackingWaveList) == 0)
		Abort "No particle tracks seem to have been calculated"
	endif
	
	string trackingWaveName
	Prompt trackingWaveName, "Track wave:", popup, trackingWaveList
	DoPrompt "Select the tracks to filter", trackingWaveName
	if (V_flag != 0)
		return 0
	endif
	
	wave /WAVE /Z trackWave = GetTrackingWaveReference(trackingWaveName)
	if (WaveExists(trackWave) == 0)
		Abort "Unable to find the requested track"
	endif


	VelocityHistogramDialog(trackWave)
End

Function ExtractMSDHistograms()
	DFREF packagekageFolder = root:Packages:Localizer

	string trackingWaveList = GetPossibleTrackingWaves()
	if (ItemsInList(trackingWaveList) == 0)
		Abort "No particle tracks seem to have been calculated"
	endif
	
	string trackingWaveName
	Prompt trackingWaveName, "Track wave:", popup, trackingWaveList
	DoPrompt "Select the tracks to filter", trackingWaveName
	if (V_flag != 0)
		return 0
	endif
	
	wave /WAVE /Z trackWave = GetTrackingWaveReference(trackingWaveName)
	if (WaveExists(trackWave) == 0)
		Abort "Unable to find the requested track"
	endif


	MSDHistogramDialog(trackWave)
End


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------


























Function GenerateScatterPlot_Menu()
	string posName
	string windowName
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	SVAR /Z gScatterColorTable = packageFolder:V_scatterColorTable
	if (SVAR_Exists(gScatterColorTable) != 1)
		string /G packageFolder:V_scatterColorTable
	endif
	SVAR gScatterColorTable = packageFolder:V_scatterColorTable
	
	string colorScaleString = gScatterColorTable
	
	Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
	Prompt colorScaleString, "Color scale (shows when positions were localized):", popup, "None;" + CTabList()
	
	DoPrompt "Select the positions wave", posName, colorScaleString
	if (V_flag == 1)	// cancel
		return 0
	endif
	
	if (StringMatch(colorScaleString, "None") == 1)
		colorScaleString = ""
	endif
	gScatterColorTable = colorScaleString
	
	wave /Z pos = GetPositionsWaveReference(posName)
	if (WaveExists(pos) == 0)
		Abort "The selected wave doesn't seem to exist!"
	endif
	
	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(pos, xSize, ySize, pixelSize)
	
	DFREF savDF = GetDataFolderDFR()
	NewDataFolder /O root:'Localizer Images'
	NewDataFolder /S/O root:'Localizer Images':$posName
	
	// do the actual work
	GenerateScatterPlot(pos, xSize, ySize, pixelSize)
	
	wave Scatter_positions_X =  Scatter_positions_X 
	wave Scatter_positions_Y =  Scatter_positions_Y
	wave Scatter_positions_color = Scatter_positions_color
	SetDataFolder savDF
	
	windowName = "Scatter plot from " + posName
	
	Display /K=1 /N=Scatter_Viewer Scatter_positions_Y vs Scatter_positions_X as windowName
	if (strlen(colorScaleString) != 0)
		ModifyGraph zColor(Scatter_positions_Y)={Scatter_positions_color,*,*,$colorScaleString,0}
	endif
	ModifyGraph mode=2,tickUnit=1
	ModifyGraph /W=$S_name width={Plan,1,bottom,left}
	
	if (pixelSize == 0)
		SetAxis /W=$S_name bottom, 0, xSize - 1
		SetAxis /W=$S_name left, 0, ySize - 1
		Label /W=$S_name left "Distance (pixel)"
		Label /W=$S_name bottom "Distance (pixel)"
	else
		SetAxis /W=$S_name bottom, 0, (xSize - 1) * pixelSize / 1000
		SetAxis /W=$S_name left, 0, (ySize - 1) * pixelSize / 1000
		Label /W=$S_name left "Distance (µm)"
		Label /W=$S_name bottom "Distance (µm)"
	endif
	
End

Function Generate3DScatterPlot_Menu()

	string posName
	string windowName
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	SVAR /Z g3DScatterColorTable = packageFolder:S_3DSscatterColorTable
	if (SVAR_Exists(g3DScatterColorTable) != 1)
		string /G packageFolder:S_3DSscatterColorTable
	endif
	SVAR g3DScatterColorTable = packageFolder:S_3DSscatterColorTable
	string colorScaleString = g3DScatterColorTable
	
	Prompt posName, "Positions wave:", popup, GetPossible3DPositionsWaves()
	Prompt colorScaleString, "Color scale (shows z position):", popup, CTabList()
	
	DoPrompt "Select the positions wave", posName, colorScaleString
	if (V_flag == 1)	// cancel
		return 0
	endif
	
	g3DScatterColorTable = colorScaleString
	
	wave /Z pos = GetPositionsWaveReference(posName)
	if (WaveExists(pos) == 0)
		Abort "The selected wave doesn't seem to exist!"
	endif
	
	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(pos, xSize, ySize, pixelSize)
	
	DFREF savDF = GetDataFolderDFR()
	NewDataFolder /O root:'Localizer Images'
	NewDataFolder /S/O root:'Localizer Images':$posName
	
	// do the actual work
	Generate3DScatterPlot(pos, xSize, ySize, pixelSize, colorScaleString)
	
	wave Scatter_positions_X =  Scatter_positions_X 
	wave Scatter_positions_Y =  Scatter_positions_Y
	wave Scatter_positions_Z =  Scatter_positions_Z
	SetDataFolder savDF
	
	windowName = "3D scatter plot from " + posName
	
	Display /K=1 /N=Scatter_Viewer Scatter_positions_Y vs Scatter_positions_X as windowName
	ModifyGraph mode=2,tickUnit=1
	ModifyGraph /W=$S_name zColor(Scatter_positions_Y)={Scatter_positions_Z,*,*,Rainbow,0}
	ModifyGraph /W=$S_name width={Plan,1,bottom,left}, height=0
	
	if (pixelSize == 0)
		SetAxis /W=$S_name bottom, 0, xSize - 1
		SetAxis /W=$S_name left, 0, ySize - 1
		Label /W=$S_name left "Distance (pixel)"
		Label /W=$S_name bottom "Distance (pixel)"
		ColorScale /W=$S_name/C/N=text0 "Z position (pixel)"
	else
		SetAxis /W=$S_name bottom, 0, (xSize - 1) * pixelSize / 1000
		SetAxis /W=$S_name left, 0, (ySize - 1) * pixelSize / 1000
		Label /W=$S_name left "Distance (µm)"
		Label /W=$S_name bottom "Distance (µm)"
		ColorScale /W=$S_name/C/N=text0 "Z position (µm)"
	endif
	
	//ColorScale /W=$S_name 
End

Function MakeAccumulatedImage_Menu()
	Variable /G root:Packages:Localizer:V_AccumulatedBinSize
	
	NVAR binSize = root:Packages:Localizer:V_AccumulatedBinSize
	if (binSize == 0)
		binSize = 0.2
	endif
	
	string posName
	string windowTitle
	variable colorScaleIndex, localBinSize
	string colorScaleString
	
	localBinSize = binSize
	Prompt localBinSize, "Bin size (pixels):"
	Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
	Prompt colorScaleString, "Color scale:", popup, CTabList()
	DoPrompt "Select the pixel size and positions wave", posName, localBinSize, colorScaleString
	if (V_flag == 1)	// cancel
		return 0
	endif
	binSize = localBinSize
	
	wave /Z pos = GetPositionsWaveReference(posName)
	
	if (WaveExists(pos) == 0)
		Abort "The selected wave doesn't seem to exist!"
	endif
	
	NewDataFolder /O root:'Localizer Images'
	NewDataFolder /O root:'Localizer Images':$posName
	DFREF outputDataFolder = root:'Localizer Images':$posName
	
	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(pos, xSize, ySize, pixelSize)
	
	wave M_FreeAccumulatedImage = MakeAccumulatedImage(pos, xSize, ySize, binSize, pixelSize=pixelSize)
	Duplicate /O M_FreeAccumulatedImage, outputDataFolder:M_AccumulatedImage
	wave M_AccumulatedImage = outputDataFolder:M_AccumulatedImage
	
	windowTitle = "Accumulated image from " + posName
	
	// Typically the histogram for the resulting image is very inhomogeneous
	// a reasonable guess for an initial colorscale appears to be to set min to auto
	// and then take 10% of the maximum for the upper limit
	// this is similar to the bitmap image
	WaveStats /Q /M=1 M_AccumulatedImage
	
	Display /K=1 /N=AccumulatedImageViewer as windowTitle
	AppendImage /W=$S_name M_AccumulatedImage
	AdjustGraphForImageDisplay(S_name)
	ModifyImage /W=$S_name M_AccumulatedImage ctab= {*,V_max / 10,$colorScaleString,0}
	ModifyGraph /W=$S_name width={Plan,1,bottom,left},tickUnit=1
	
	if (pixelSize == 0)
		Label /W=$S_name left "Distance (pixels)"
		Label /W=$S_name bottom "Distance (pixels)"
	else
		Label /W=$S_name left "Distance (µm)"
		Label /W=$S_name bottom "Distance (µm)"
	endif
	
End

Function MakePALMImage_Menu()
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	
	Variable /G root:Packages:Localizer:V_PALMBitmapConstantWidth
	Variable /G root:Packages:Localizer:V_PALMBitmapFitUncertainty
	Variable /G root:Packages:Localizer:V_PALMBitmapCamMagnification
	Variable /G root:Packages:Localizer:V_PALMBitmapCamOffset
	Variable /G root:Packages:Localizer:V_PALMCCDPixelSize
	Variable /G root:Packages:Localizer:V_PALMBitmapImageWidth
	Variable /G root:Packages:Localizer:V_gaussianWidth
	Variable /G root:Packages:Localizer:V_PALMBitmapSelectedMethod
	
	NVAR constantWidth = root:Packages:Localizer:V_PALMBitmapConstantWidth
	NVAR fitUncertainty = root:Packages:Localizer:V_PALMBitmapFitUncertainty
	NVAR camMagnification = root:Packages:Localizer:V_PALMBitmapCamMagnification
	NVAR camOffset = root:Packages:Localizer:V_PALMBitmapCamOffset
	NVAR imageWidth = root:Packages:Localizer:V_PALMBitmapImageWidth
	NVAR PSFWidth = root:Packages:Localizer:V_gaussianWidth
	NVAR SelectedMethod = root:Packages:Localizer:V_PALMBitmapSelectedMethod
	
	// provide initial values if this is the first time that we've used this function
	if (constantWidth == 0)
		constantWidth = 0.2
	endif
	if (fitUncertainty == 0)
		fitUncertainty = 3
	endif
	if (imageWidth == 0)
		imageWidth = 2096
	endif
	if (PSFWidth == 0)
		PSFWidth = 2
	endif
	if (camMagnification == 0)
		camMagnification = 1
	endif
	if (camOffset == 0)
		camOffset = 1000
	endif
	
	string listPositionsFunc = GetIndependentModuleName() + "#GetPossiblePositionsWaves()"
	
	DoWindow /F PALMBitmapWindow
	if (V_flag == 0)
		NewPanel /N=PALMBitmapWindow /K=1 /W=(137,135,397,490) as "Enter calculation parameters"
		GroupBox GBPALMBitmapColorBox,pos={4,65},size={250,25}
		PopupMenu PMPALMBitmapSelectPositionsWave,pos={14,8},size={220,20},bodyWidth=150,title="Positions wave:"
		PopupMenu PMPALMBitmapSelectPositionsWave,mode=1,value= #listPositionsFunc
		PopupMenu PMPALMBitmapSelectColorTable,pos={35,33},size={198,20},bodyWidth=145,title="Color table:"
		PopupMenu PMPALMBitmapSelectColorTable,mode=7,value= #"\"*COLORTABLEPOPNONAMES*\""
		
		GroupBox GBPALMBitmapEmitterWeighing,pos={4,102},size={250,57},title="Relative position weight:"
		CheckBox CBPALMBitmapSameEmitterWeight,pos={11,119},size={209,14},proc=CBPALMBitmapImage,title="Each localized position has the same weight"
		CheckBox CBPALMBitmapSameEmitterWeight,value= 1,mode=1
		CheckBox CBPALMBitmapPreserveIntegral,pos={11,139},size={188,14},proc=CBPALMBitmapImage,title="Preserve the fitted integrated intensity"
		CheckBox CBPALMBitmapPreserveIntegral,value= 0,mode=1
		
		GroupBox GBPALMBitmapDeviationBox,pos={3,171},size={250,150},title="Spot deviation:"
		CheckBox CBPALMBitmapConstantWidth,pos={10,191},size={168,14},proc=CBPALMBitmapImage,title="Constant width:\t\t\t\tCCD pixels"
		CheckBox CBPALMBitmapConstantWidth,value= 0,mode=1
		CheckBox CBPALMBitmapFitUncertainty,pos={10,214},size={166,14},proc=CBPALMBitmapImage,title="Scale fit uncertainty\t\t\t\ttimes"
		CheckBox CBPALMBitmapFitUncertainty,value= 0,mode=1
		CheckBox CBPALMBitmapGME,pos={10,239},size={130,14},proc=CBPALMBitmapImage,title="Gaussian mask estimation"
		CheckBox CBPALMBitmapGME,value= 0,mode=1
		SetVariable SVPALMBitmapImageWidth,pos={7,69},size={140,15},title="Image width (pixels)"
		SetVariable SVPALMBitmapImageWidth,limits={32,inf,0},value=imageWidth
		SetVariable SVPALMBitmapConstantWidth,pos={103,191},size={30,15},title=" "
		SetVariable SVPALMBitmapConstantWidth,limits={0,inf,0},value=constantWidth
		SetVariable SVPALMBitmapFitUncertainty,pos={121,214},size={30,15},title=" "
		SetVariable SVPALMBitmapFitUncertainty,limits={0,inf,0},value=fitUncertainty
		SetVariable SVPALMCameraMagnification,pos={29,258},size={190,15},title="Camera amplification factor:"
		SetVariable SVPALMCameraMagnification,limits={1e-4,inf,0},value=camMagnification
		SetVariable SVPALMCameraOffset,pos={28,278},size={190,15},title="Camera offset:"
		SetVariable SVPALMCameraOffset,limits={0,inf,0},value=camOffset
		SetVariable SVPALMBitmapPSFWidth,pos={28,298},size={190,15},title="PSF width (in CCD pixels):"
		SetVariable SVPALMBitmapPSFWidth,limits={0,inf,0},value=PSFWidth
		Button BTPALMBitmapImageCalculate,pos={179,329},size={75,20},proc=BTPALMBitmapImageCalculateProc,title="Calculate!"
		
		// if the window was opened previously then restore its size and position
		WC_WindowCoordinatesRestore("PALMBitmapWindow")
		SetWindow PALMBitmapWindow hook(WindowCoordsHook)=WC_WindowCoordinatesNamedHook
	endif
	
	// restore the selected method if we used it previously
	Switch (SelectedMethod)
		case 0:
			CheckBox CBPALMBitmapConstantWidth value= 1
			break
		case 1:
			CheckBox CBPALMBitmapFitUncertainty value= 1
			break
		case 2:
			CheckBox CBPALMBitmapGME value = 1
			break
	EndSwitch
	
End

Function BTPALMBitmapImageCalculateProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			NVAR constantWidth = root:Packages:Localizer:V_PALMBitmapConstantWidth
			NVAR uncertaintyScale = root:Packages:Localizer:V_PALMBitmapFitUncertainty
			NVAR camOffset = root:Packages:Localizer:V_PALMBitmapCamOffset
			NVAR camMagnification = root:Packages:Localizer:V_PALMBitmapCamMagnification
			NVAR imageWidth = root:Packages:Localizer:V_PALMBitmapImageWidth
			NVAR PSFWidth = root:Packages:Localizer:V_gaussianWidth
			string posName
			string colorTable
			variable method = -1
			variable weighingMethod = 0
			ControlInfo /W=PALMBitmapWindow PMPALMBitmapSelectPositionsWave
			posName = S_Value
			
			ControlInfo /W=PALMBitmapWindow PMPALMBitmapSelectColorTable
			colorTable = S_value
			
			ControlInfo /W=PALMBitmapWindow CBPALMBitmapSameEmitterWeight
			if (V_value == 1)
				weighingMethod = 0
			endif
			
			ControlInfo /W=PALMBitmapWindow CBPALMBitmapPreserveIntegral
			if (V_value == 1)
				weighingMethod = 1
			endif
			
			ControlInfo /W=PALMBitmapWindow CBPALMBitmapConstantWidth
			if (V_value == 1)
				method = PALMBITMAP_DEV_SAME
			endif
			ControlInfo /W=PALMBitmapWindow CBPALMBitmapFitUncertainty
			if (V_value == 1)
				method = PALMBITMAP_DEV_FITUNCERTAINTY
			endif
			ControlInfo /W=PALMBitmapWindow CBPALMBitmapGME
			if (V_value == 1)
				method = PALMBITMAP_DEV_GAUSSIANMASK
			endif
				
			
			
			wave /Z pos = GetPositionsWaveReference(posName)
			
			if (WaveExists(pos) == 0)
				Abort "The selected wave doesn't seem to exist!"
			endif
			
			DFREF savDF = GetDataFolderDFR()
			NewDataFolder /O root:'Localizer Images'
			NewDataFolder /O/S root:'Localizer Images':$posName
			
			// if there is already a PALM image in this data folder then the calculation will overwrite it
			// try to kill it before the actual calculation to maybe free up some memory
			KillWaves /Z $"M_LocalizationBitmap"
			
			variable xSize, ySize, pixelSize
			GetCCDDimensionsFromPositions(pos, xSize, ySize, pixelSize)
			
			variable imageScaleFactor = round(imageWidth / xSize)
			variable upperLimit = xSize
			
			// add a structure with a reference to the progress reporting function
			STRUCT LocalizerProgStruct progressStruct
			progressStruct.version = kLocalizerProgStructVersion
			FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
			
			switch (method)
				case PALMBITMAP_DEV_SAME:
					LocalizationBitmap /M=(method) /WGHT=(weighingMethod) /S=(constantWidth) /W={xSize, ySize, imageScaleFactor} /PROG=progressStruct pos
					break
				case PALMBITMAP_DEV_FITUNCERTAINTY:
					LocalizationBitmap /M=(method) /WGHT=(weighingMethod) /S=(uncertaintyScale) /L=(upperLimit) /W={xSize, ySize, imageScaleFactor} /PROG=progressStruct pos
					break
				case PALMBITMAP_DEV_GAUSSIANMASK:
					LocalizationBitmap /M=(method) /WGHT=(weighingMethod) /WDTH=(PSFWidth) /CAL={camOffset, camMagnification} /W={xSize, ySize, imageScaleFactor} /PROG=progressStruct pos
					break
				default:
					SetDataFolder savDF
					Abort "Unknown method"
					break
			endswitch
			
			wave M_LocalizationBitmap
			string windowTitle = "PALM bitmap image from " + posName
			
			// Typically the histogram for the resulting image is very inhomogeneous
			// a reasonable guess for an initial colorscale appears to be to set min to auto
			// and then take 10% of the maximum for the upper limit
			WaveStats /Q /M=1 M_LocalizationBitmap
			
			Display /K=1 /N=PALMImageViewer as windowTitle
			AppendImage /W=$S_name M_LocalizationBitmap
			AdjustGraphForImageDisplay(S_name)
			ModifyImage /W=$S_name M_LocalizationBitmap ctab= {*, V_max / 10,$colorTable,0}
			ModifyGraph /W=$S_name width={Plan,1,bottom,left}
			
			SetDataFolder savDF
			
			if (pixelSize == 0)
				SetScale /I x, 0, xSize -1, "pixel", M_LocalizationBitmap
				SetScale /I y, 0, ySize -1, "pixel", M_LocalizationBitmap
				Label /W=$S_name left "Distance (pixels)"
				Label /W=$S_name bottom "Distance (pixels)"
			else
				SetScale /I x, 0, (xSize -1) * pixelSize / 1000, "um", M_LocalizationBitmap
				SetScale /I y, 0, (ySize -1) * pixelSize / 1000, "um", M_LocalizationBitmap
				Label /W=$S_name left "Distance (µm)"
				Label /W=$S_name bottom "Distance (µm)"
			endif
			
			AppendColorRangeSliders()
			
			DoWindow /K PALMBitmapWindow
					
			break
	endswitch

	return 0
End

Function CBPALMBitmapImage(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			String controlName = cba.ctrlName
			Variable checked = cba.checked
			NVAR selectedMethod = root:Packages:Localizer:V_PALMBitmapSelectedMethod
			
			StrSwitch(controlName)
				// checkboxes that affect how different positions are weighed
				case "CBPALMBitmapSameEmitterWeight":
					CheckBox CBPALMBitmapPreserveIntegral value = 0
					break
				case "CBPALMBitmapPreserveIntegral":
					CheckBox CBPALMBitmapSameEmitterWeight value = 0
					break
				
				
				// checkboxes that affect the calculation of the spot standard deviation
				case "CBPALMBitmapConstantWidth":
					CheckBox CBPALMBitmapFitUncertainty value=0
					CheckBox CBPALMBitmapGME value=0
					selectedMethod = 0
					break
				case "CBPALMBitmapFitUncertainty":
					CheckBox CBPALMBitmapConstantWidth value=0
					CheckBox CBPALMBitmapGME value=0
					selectedMethod = 1
					break
				case "CBPALMBitmapGME":
					CheckBox CBPALMBitmapConstantWidth value=0
					CheckBox CBPALMBitmapFitUncertainty value=0
					selectedMethod = 2
					break
			EndSwitch
			
			break
	endswitch

	return 0
End

Function SaveCCDImagesAs()
	
	string windowName = GetNameOfTopInterfacePanel()
	if (strlen(windowName) == 0)
		Abort "No data file seems to have been opened, or the top panel does not belong tot this package"
	endif
	
	DFREF windowDataFolder = GetWindowDataFolder(windowName)
	SVAR CCDFilePath = windowDataFolder:S_filePath
	variable fileType
	NVAR /Z gFileType = windowDataFolder:V_fileType
	if (!NVAR_Exists(gFileType))
		fileType = -1
	else
		fileType = gFileType
	endif
	
	string outputPath
	variable outputFormat
	GetOutputStorageTypeAndPath(CCDFilePath, outputFormat, outputPath, offerPresentationOutput = 1)
	if (strlen(outputPath) == 0)
		return 0
	endif
	if (outputFormat == IMAGE_OUTPUT_TYPE_MOV)
		// presentation movie output is handled by a different function, not by ProcessCCDImages
		MakeMovieFromOpenFile(windowName)
		return 0
	endif
	
	// add a structure with a reference to the progress reporting function
	STRUCT LocalizerProgStruct progressStruct
	progressStruct.version = kLocalizerProgStructVersion
	FUNCREF ProgressFunctionPrototype progressStruct.func = ProgressReporterFunc
	
	ProcessCCDImages /O /Y=(fileType) /M=(PROCESSCCDIMAGE_CONVERTFORMAT) /OUT=(outputFormat) /PROG=progressStruct CCDFilePath, outputPath
	
End

Function MakeParticleTrackingScatterPlot()
	string trackingWaveList = GetPossibleTrackingWaves()
	if (ItemsInList(trackingWaveList) == 0)
		Abort "No particle tracks seem to have been calculated"
	endif
	
	string trackingWaveName, colorTableName
	Prompt trackingWaveName, "Track wave:", popup, trackingWaveList
	Prompt colorTableName, "Color table:", popup, "None;" + CTabList()
	DoPrompt "Select the tracks to display", trackingWaveName, colorTableName
	if (V_flag != 0)
		return 0
	endif
	
	wave /WAVE /Z trackWave = GetTrackingWaveReference(trackingWaveName)
	if (WaveExists(trackWave) == 0)
		Abort "Unable to find the requested track"
	endif
	
	variable nTracks = DimSize(trackWave, 0)
	
	variable xSize, ySize, pixelSize
	GetCCDDimensionsFromPositions(trackWave, xSize, ySize, pixelSize)
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(trackWave, xCol, yCol, zCol)
	
	NewDataFolder /O root:'Tracking Images'
	NewDataFolder /O root:'Tracking Images':$trackingWaveName
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF imageFolder = root:'Tracking Images':$trackingWaveName
	DFREF savDF = GetDataFolderDFR()
	
	Make /O/N=0 /D imageFolder:W_trackScatter_X, imageFolder:W_trackScatter_Y
	wave W_trackScatter_X = imageFolder:W_trackScatter_X
	wave W_trackScatter_Y = imageFolder:W_trackScatter_Y
	
	variable useColorTable = 0
	variable nColors
	if (StringMatch(colorTableName, "None") == 0)	
		useColorTable = 1
	endif
	
	wave M_ScatterCoordinates = ParticleTrackToScatterPairs(trackWave)
	variable nCoordinates = DimSize(M_ScatterCoordinates, 0)
	Make /O/N=(nCoordinates) /D imageFolder:W_trackScatter_X, imageFolder:W_trackScatter_Y, imageFolder:W_trackScatter_Colors
	wave W_trackScatter_X = imageFolder:W_trackScatter_X
	wave W_trackScatter_Y = imageFolder:W_trackScatter_Y
	wave W_trackScatter_Colors = imageFolder:W_trackScatter_Colors
	W_trackScatter_X = M_ScatterCoordinates[p][0]
	W_trackScatter_Y = M_ScatterCoordinates[p][1]
	W_trackScatter_Colors = M_ScatterCoordinates[p][2]

	Display /K=1 /N=TrackingScatter W_trackScatter_Y vs W_trackScatter_X
	DoWindow /T $S_name, "All tracks from " + NameOfWave(trackWave)
	ModifyGraph /W=$S_name width={Plan,1,bottom,left}
	if (useColorTable != 0)
		ModifyGraph /W=$S_name zColor(W_trackScatter_Y)={W_trackScatter_Colors,*,*,$colorTableName,0}
	endif
	
	// take the pixel size into account
	if (pixelSize == 0)
		SetAxis /W=$S_name bottom, 0, xSize - 1
		SetAxis /W=$S_name left, 0, ySize - 1
		Label /W=$S_name left "Distance (pixel)"
		Label /W=$S_name bottom "Distance (pixel)"
	else
		W_trackScatter_X *= pixelSize / 1000
		W_trackScatter_Y *= pixelSize / 1000
		SetScale /P x, 0, 1, "um", W_trackScatter_X, W_trackScatter_Y
		SetAxis /W=$S_name bottom, 0, (xSize - 1) * pixelSize / 1000
		SetAxis /W=$S_name left, 0, (ySize - 1) * pixelSize / 1000
		Label /W=$S_name left "Distance (µm)"
		Label /W=$S_name bottom "Distance (µm)"
	endif
End

Function MakeMovieFromOpenFile(windowName)
	String windowName
	
	Assert(strlen(windowName) > 0)
	
	string baseWindowName = GetBaseWindowName(windowName)
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	
	SVAR dataFilePath = windowDataFolder:S_filePath
	NVAR nCCDFrames = windowDataFolder:V_nImages
	NVAR currentImage = windowDataFolder:V_currentImage
	
	variable firstFrameInMovie = 0, lastFrameInMovie = nCCDFrames - 1
	variable frameDelta = 1
	variable originalFrame = currentImage
	string movieFormat
	variable fps = 25, refnum, i
	string outputFilePath
	string error
	
	Prompt firstFrameInMovie, "First frame in movie:"
	Prompt lastFrameInMovie, "Last frame in movie:"
	Prompt frameDelta, "Display every nth frame:"
	#ifdef MACINTOSH
		Prompt movieFormat, "Movie format:", popup, "Quicktime;"
	#else
		Prompt movieFormat, "Movie format:", popup, "Quicktime;AVI;"
	#endif
	
	DoPrompt "Enter movie parameters", firstFrameInMovie, lastFrameInMovie, frameDelta, movieFormat
	if (V_Flag == 1)
		return 1
	endif
	
	string movieExtension
	StrSwitch(movieFormat)
		case "Quicktime":
			movieExtension = ".mov"
			break
		case "AVI":
			movieExtension = ".avi"
			break
		default:
			Abort "Unknown extension for movie format"
			break
	EndSwitch
	
	// remove the extension of the data file and replace it with .mov
	outputFilePath  = ParseFilePath(3, dataFilePath, ":", 0, 0)
	outputFilePath += movieExtension
	
	if (firstFrameInMovie < 0)
		Abort "Invalid first frame"
	elseif (lastFrameInMovie > nCCDFrames)
		Abort "The data contains less frames than requested"
	elseif (firstFrameInMovie >= lastFrameInMovie)
		Abort "Invalid frame range"
	endif
	
	Open /D /T=movieExtension refnum as outputFilePath
	if (strlen(S_fileName) == 0)
		return 1
	endif
	
	outputFilePath = S_fileName
	
	// remove a pre-exisiting movie viewer if it exists
	DoWindow /K MovieViewer
	
	// get a recreation macro for the graph with the movie
	// and parse it to get the macro for only the graph
	String recreationMacro = WinRecreation(baseWindowName, 0)
	recreationMacro = ReplaceString("\r\t", recreationMacro, ";")
	
	// don't get all of the macro, leave out the display command
	variable startOfMacro = strsearch(recreationMacro, ";Display", 0)
	startOfMacro = strsearch(recreationMacro, ";", startOfMacro + 1) + 1
	variable endOfMacro = strsearch(recreationMacro, "RenameWindow", startOfMacro) - 1
	
	string relevantMacro = recreationMacro[startOfMacro, endOfMacro]
	
	// parsing the recreation macro is surprisingly annoying. It turns out that when e.g.
	// particle tracks are displayed, a SetDataFolder call is inserted before the AppendImage
	// and all other appends are relative to this one. So try to kludge around that
	variable firstSetDataFolder = strsearch(recreationMacro, "SetDataFolder", 0)
	if ((firstSetDataFolder != -1) && (firstSetDataFolder < startOfMacro))
		string dataFolderSetTo = recreationMacro[strsearch(recreationMacro, "root:", firstSetDataFolder), strsearch(recreationMacro, ";Display", firstSetDataFolder) - 1]
		// even more annoying is the fact that the Display command will list a wave to append directly
		// so get that one out, too
		string displayCommand = recreationMacro[strsearch(recreationMacro, ";Display", firstSetDataFolder), strsearch(recreationMacro, ";Append", firstSetDataFolder) - 1]
		string waveToAppend
		SplitString /E=".[^ ]* (.*)$" displayCommand, waveToAppend
		
		// and generate the append command and prepend it to the macro
		string savDF = GetDataFolder(1)
		relevantMacro = "SetDataFolder " + dataFolderSetTo + "; AppendToGraph " + waveToAppend + "; " + relevantMacro
		relevantMacro = ReplaceString("SetDataFolder fldrSav0;", relevantMacro, "SetDataFolder " + savDF + ";")
	endif
		
	
	// we'll record the movie with the correct aspect size, so get hold of the image dimension
	wave image = windowDataFolder:M_CCDImage
	variable xSize = DimSize(image, 0)
	variable ySize = DimSize(image, 1)
	
	// display the graph
	DoWindow /K MovieViewer
	Display /N=MovieViewer /W=(67,128, 67 + xSize * ScreenResolution / 72, 128 + ySize * ScreenResolution / 72)
	
	// execute the recreation macro in a piecewise fashion
	// otherwise the string to execute can become too long
	for (i = 0; i < ItemsInList(relevantMacro); i += 1)
		string Dummy = StringFromList(i, relevantMacro)
		Execute StringFromList(i, relevantMacro)
	endfor
	
	// get rid of margins and axes
	ModifyGraph /W=MovieViewer margin=-1
	ModifyGraph /W=MovieViewer noLabel=2,axThick=0
	
	if (StringMatch(movieFormat, "AVI") == 1)
		NewMovie /A /Z /O /L /F=(fps) /I as outputFilePath
	else
		NewMovie /Z /O /L /F=(fps) /I as outputFilePath
	endif
	if (V_flag == -1)
		// user canceled
		DoWindow /K MovieViewer
		return 0
	endif
	if (V_flag != 0)
		error = GetErrMessage(V_flag, 3)
		DoWindow /K MovieViewer
		Abort error
	endif
	
	for (i = firstFrameInMovie; i < lastFrameInMovie; i += frameDelta)
		ChangeCurrentImage(i, baseWindowName)
		currentImage = i
		DoUpdate
		AddMovieFrame
	endfor
	
	CloseMovie
	DoWindow /K MovieViewer
	
	// restore the image we were originally looking at
	ChangeCurrentImage(originalFrame, baseWindowName)
End

Function SetupBatchLocalizationPanel()
	
	// check if there is a viewer panel open
	// this is required because we take the analysis settings
	// from the topmost panel
	string topWinName = GetNameOfTopInterfacePanel()
	if (strlen(topWinName) == 0)
		Abort "Please open a data file and set the desired processing options in the panel"
	endif
	string baseWindowName = GetBaseWindowName(topWinName)
	
	Make /O/N=(0, 2) /T root:Packages:Localizer:M_BatchProcessListWave
	Make /O/B/U/N=(0, 2, 2) root:Packages:Localizer:M_BatchProcessSelectionWave
	Make /O/W/U/N=(2, 3) root:Packages:Localizer:M_BatchProcessColorWave
	Make /O/N=2 /T root:Packages:Localizer:M_BatchProcessColumnHeadings
	wave /T M_BatchProcessColumnHeadings = root:Packages:Localizer:M_BatchProcessColumnHeadings
	wave M_BatchProcessColorWave = root:Packages:Localizer:M_BatchProcessColorWave
	wave M_BatchProcessSelectionWave = root:Packages:Localizer:M_BatchProcessSelectionWave
	
	M_BatchProcessColumnHeadings[0] = "CCD Files"
	M_BatchProcessColumnHeadings[1] = "Output position names"
	
	M_BatchProcessColorWave = 0
	M_BatchProcessColorWave[1][0] = 65535
	SetDimLabel 2,1, backColors, M_BatchProcessSelectionWave
	M_BatchProcessSelectionWave[][][1] = 0	// provide the default Igor colors
	
	DoWindow /F BatchProcess
	if (V_flag == 0)
		NewPanel /K=1 /N=BatchProcess /W=(367,181,920,534) as "Select the files to localize (double click to remove)"
		SetDrawLayer UserBack
		SetDrawEnv fsize= 11,fstyle= 2
		DrawText 9,20,"Note: all files will be processed with the settings in the top viewer panel"
		ListBox LBBatchListLocalization,pos={7,27},size={541,246}, proc=LBDeleteDoubleClickedRowProc, listWave=root:Packages:Localizer:M_BatchProcessListWave
		ListBox LBBatchListLocalization,frame=2,userColumnResize=1,titleWave=root:Packages:Localizer:M_BatchProcessColumnHeadings
		ListBox LBBatchListLocalization,mode=1,selWave=root:Packages:Localizer:M_BatchProcessSelectionWave
		ListBox LBBatchListLocalization,colorWave = root:Packages:Localizer:M_BatchProcessColorWave
		ListBox LBBatchListLocalization,userdata(ResizeControlsInfo)= A"!!,@C!!#=;!!#Cl5QF/pz!!#](Aon\"Qzzzzzzzzzzzzzz!!#o2B4uAezz"
		ListBox LBBatchListLocalization,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#u:Du]k<zzzzzzzzzzz"
		ListBox LBBatchListLocalization,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		GroupBox GBSpecialOptions,pos={7,282},size={541,43},title="Special options"
		GroupBox GBSpecialOptions,userdata(ResizeControlsInfo)= A"!!,@C!!#BG!!#Cl5QF,%z!!#](Aon\"Qzzzzzzzzzzzzzz!!#o2B4uAezz"
		GroupBox GBSpecialOptions,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#?(FEDG<zzzzzzzzzzz"
		GroupBox GBSpecialOptions,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		CheckBox CBBatchAnalyzeAverage,pos={15,303},size={303,14},title="Analyze average image for each dataset (for e.g. multicolor or 3D calibration)"
		CheckBox CBBatchAnalyzeAverage,value= 0
		CheckBox CBBatchAnalyzeAverage,userdata(ResizeControlsInfo)= A"!!,B)!!#BQJ,hs'J,hlCz!!#](Aon\"Qzzzzzzzzzzzzzz!!#o2B4uAezz"
		CheckBox CBBatchAnalyzeAverage,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#?(FEDG<zzzzzzzzzzz"
		CheckBox CBBatchAnalyzeAverage,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		CheckBox CBBatchSaveExperiment,pos={6,332},size={175,14},title="Autosave experiment after each run"
		CheckBox CBBatchSaveExperiment,value= 1
		CheckBox CBBatchSaveExperiment,userdata(ResizeControlsInfo)= A"!!,@#!!#B`!!#A>!!#;mz!!#](Aon\"Qzzzzzzzzzzzzzz!!#o2B4uAezz"
		CheckBox CBBatchSaveExperiment,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#?(FEDG<zzzzzzzzzzz"
		CheckBox CBBatchSaveExperiment,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		Button BTBatchStartAnalysis,pos={434,329},size={100,20},title="Start!",proc=BTStartBatchLocalization
		Button BTBatchStartAnalysis,userdata(ResizeControlsInfo)= A"!!,I?!!#B^J,hpW!!#<Xz!!#o2B4uAezzzzzzzzzzzzzz!!#o2B4uAezz"
		Button BTBatchStartAnalysis,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#?(FEDG<zzzzzzzzzzz"
		Button BTBatchStartAnalysis,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		Button BTBatchAddLocalizationFiles,pos={320,329},size={100,20},title="Add files",proc=BTBatchAddDataFiles
		Button BTBatchAddLocalizationFiles,userdata(ResizeControlsInfo)= A"!!,H[!!#B^J,hpW!!#<Xz!!#o2B4uAezzzzzzzzzzzzzz!!#o2B4uAezz"
		Button BTBatchAddLocalizationFiles,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#?(FEDG<zzzzzzzzzzz"
		Button BTBatchAddLocalizationFiles,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		
		// the controls in this panel will need to know the name of the panel
		// with the settings. Encode this in the userdata
		SetWindow kwTopWin, userdata = baseWindowName
		SetWindow kwTopWin,userdata(ResizeControlsInfo)= A"!!*'\"z!!#Co5QF0UJ,fQLzzzzzzzzzzzzzzzzzzzz"
		SetWindow kwTopWin,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzzzzzzzzzzzzzzz"
		SetWindow kwTopWin,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzzzzzzzzz!!!"
		
		// make the controls resize dynamically when the window is resized thanks to the <Resize Controls> package
		SetWindow BatchProcess,hook(ResizeControls)=ResizeControls#ResizeControlsHook
		
		// if the window was opened previously then restore its size and position
		WC_WindowCoordinatesRestore("BatchProcess")
		SetWindow BatchProcess hook(WindowCoordsHook)=WC_WindowCoordinatesNamedHook
	endif
	
End

Function BTBatchAddDataFiles(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	string baseFileName, extension, outputFilePath, macintoshFilePath
	variable nFiles, nNewFiles, i, offset
	
	switch( ba.eventCode )
		case 2: // mouse up
				// click code here
			wave /T M_BatchProcessListWave = root:Packages:Localizer:M_BatchProcessListWave
			wave M_BatchProcessSelectionWave = root:Packages:Localizer:M_BatchProcessSelectionWave
			wave M_BatchProcessColorWave = root:Packages:Localizer:M_BatchProcessColorWave
			
			Open /D/R /MULT=1/T="????" refNum
			if (strlen(S_fileName) == 0)
				// user canceled
				return 0
			endif
			
			offset = DimSize(M_BatchProcessListWave, 0)
			nNewFiles = ItemsInList(S_fileName, "\r") 
			nFiles = DimSize(M_BatchProcessListWave, 0) + nNewFiles
			Redimension /N=(nFiles, 2) M_BatchProcessListWave
			Redimension /N=(nFiles, 2, 2) M_BatchProcessSelectionWave
			for (i = 0; i < nNewFiles; i+=1)
				M_BatchProcessListWave[offset][0] = StringFromList(i, S_fileName, "\r")
				// the output data will be saved as igor waves
				baseFileName = ParseFilePath(3, StringFromList(i, S_fileName, "\r"), ":", 0, 0)
				//provide a valid Igor output name
				baseFileName = CleanupName(baseFileName, 0)
				M_BatchProcessListWave[offset][1] = baseFileName
				offset += 1
			endfor
			
			M_BatchProcessSelectionWave[][1][0] = 2	// be sure that the output names are editable
			
			break
	endswitch	
	
	return 0
End

Function BTStartBatchLocalization(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	variable nFiles, nNewFiles, i, j, startTime
	variable invalidOutputNames = 0
	variable identicalNames = 0
	Struct FitData fitParams
	variable err
	
	switch( ba.eventCode )
		case 2: // mouse up
				// click code here
			// get the name of the panel with the settings
			string interfacePanelName = GetUserData("BatchProcess", "", "")
			string baseWindowName = GetBaseWindowName(interfacePanelName)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			
			wave /T M_BatchProcessListWave = root:Packages:Localizer:M_BatchProcessListWave
			wave M_BatchProcessSelectionWave = root:Packages:Localizer:M_BatchProcessSelectionWave
			wave M_BatchProcessColorWave = root:Packages:Localizer:M_BatchProcessColorWave
			
			Make /T/O/N=(0, 2) root:Packages:Localizer:M_BatchProcessErrors	// keep track of data files that had errors
			wave /T M_BatchProcessErrors = root:Packages:Localizer:M_BatchProcessErrors
			
			
			nFiles = DimSize(M_BatchProcessListWave, 0)
			
			// make sure that none of the cells are highlighted
			M_BatchProcessSelectionWave[][][1] = 0	// provide the default Igor colors
				
			// are any of the output file names invalid?
			for (i = 0; i < nFiles; i+=1)
				if ((StringMatch(M_BatchProcessListWave[i][1], CleanupName(M_BatchProcessListWave[i][1], 0)) != 1) || (strlen(M_BatchProcessListWave[i][1]) == 0))
					M_BatchProcessSelectionWave[i][1][1] = 1	// give the cell a red color
					invalidOutputNames = 1
				endif
			endfor
			
			if (invalidOutputNames == 1)
				Abort "Some of the output names are invalid, please correct these"
			endif
			
			// do any of the positions have identical names?
			for (i = 0; i < nFiles - 1; i+=1)
				for (j = i + 1; j < nFiles; j+=1)
					if (StringMatch(M_BatchProcessListWave[i][1], M_BatchProcessListWave[j][1]) == 1)
						M_BatchProcessSelectionWave[i][1][1] = 1// give the cell a red color
						M_BatchProcessSelectionWave[j][1][1] = 1
						identicalNames = 1
					endif
				endfor
			endfor
			
			if (identicalNames == 1)
				Abort "Some of the output names are identical, please correct these"
			endif
			
			// if we're here then we're all set to go
			ControlInfo /W=BatchProcess CBBatchSaveExperiment
			variable autoSave = V_Value
			
			ControlInfo /W=BatchProcess CBBatchAnalyzeAverage
			variable analyzeAverageImage = V_Value
			
			// if the user selected to autosave the experiment, save it now
			// this avoids the situation where the experiment was not saved before and the 'SaveExperiment'
			// operation opens up a save dialog, but only after the first run has finished
			// and the user might have gone home
			if (autoSave != 0)
				SaveExperiment
			endif
			
			startTime = StopMSTimer(-2)
			string outputWaveName
			
			for (i = 0; i < nFiles; i+=1)
				LocalizerProgressWindow("file", i + 1, nFiles)
				// get all the fit options and parameters the user defined in the GUI
				GetGUIFitSettingsAndOptions(fitParams, baseWindowName)
				
				string thisFilePath = M_BatchProcessListWave[i][0]
				// does the user want us to fit the average?
				if (analyzeAverageImage)
					// create an average image from the dataset, and analyze that instead
					AnalyzeCCDImages /Q/M=(ANALYZECCD_AVERAGE_IMAGE) /DEST=M_BatchAvgImage thisFilePath
					fitParams.CCDFilePath = "M_BatchAvgImage"
					fitParams.firstFrameToAnalyze=0
					fitParams.lastFrameToAnalyze=0
				else
					// set the data filepath to the requested file
					fitParams.CCDFilePath = M_BatchProcessListWave[i][0]
				endif
				fitParams.cameraType = -1
				
				// run the actual fit
				err = DoPALMFitting(fitParams, baseWindowName, outputWaveName, 0)
				if (err != 0)
					if (err == kUserAbort)
						// the user wants to abort the analysis
						Print "Batch fitting stopped due to user abort"
						KillWaves /Z M_BatchProcessListWave, M_BatchProcessSelectionWave, M_BatchProcessColorWave, M_BatchProcessErrors
						return 0
					endif
					
					// if we got here then it means that there was an error fitting the current data file
					// don't abort now, but remember the filePath that caused the error and report it to the user
					Redimension /N=(DimSize(M_BatchProcessErrors, 0) + 1, DimSize(M_BatchProcessErrors, 1)) M_BatchProcessErrors
					M_BatchProcessErrors[DimSize(M_BatchProcessErrors, 0) - 1][0] = M_BatchProcessListWave[i][0]
					M_BatchProcessErrors[DimSize(M_BatchProcessErrors, 0) - 1][1] = M_BatchProcessListWave[i][1]
					
					// move on to the next file
					continue
				endif
				
				// the output positions wave will have been created with a name that depends on the filepath 
				wave POS_Out = windowDataFolder:$outputWaveName
				// copy the wave to the requested output positions name
				Duplicate /O POS_Out, root:'Localized Positions':$M_BatchProcessListWave[i][1]
				// and delete the original
				KillWaves /Z POS_Out
				
				// if the users selected to autosave the experiment, do so now
				if (autoSave != 0)
					SaveExperiment
				endif
			endfor
			
			// if we're here then the fitting is finished
			Printf "Total duration is %s\r", PrettyPrintTimeElapsed((StopMSTimer(-2) - startTime) / 1e6)
			
			// if there were any errors then report them now
			if (DimSize(M_BatchProcessErrors, 0) != 0)
				Print "WARNING: the following data file(s) reported errors and should not be trusted:"
				for (i = 0; i < DimSize(M_BatchProcessErrors, 0); i+=1)
					Printf "%d. %s,  %s %s\r", i + 1, M_BatchProcessErrors[i][0], "positions were to be written in", M_BatchProcessErrors[i][1]
				endfor
			endif
			
			
			// clean up
			KillWaves /Z M_BatchProcessListWave, M_BatchProcessSelectionWave, M_BatchProcessColorWave, M_BatchProcessErrors, M_BatchAvgImage
			
			DoWindow /K BatchProcess
			
			break
	endswitch
	
	return 0
End

Function MergePositions()
	
	string savDF = GetDataFolder(1)
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	
	DFREF packageFolder = root:Packages:Localizer
	
	string listPositionsFunc = GetIndependentModuleName() + "#GetPossiblePositionsWaves()"
	
	
	Make /T/O/N=(0) packageFolder:candidatesListBoxWave
	Make /O/N=(0) packageFolder:candidatesListBoxSelWave
	
	DoWindow /F positionsWavesViewer
	if (V_flag == 0)
		NewPanel /K=1/W=(101,286,335,497) /N=positionsWavesViewer as "Select positions waves (double click to remove)"
		ListBox ListBoxCandidateWaves,pos={1,2},size={231,183},proc=LBDeleteDoubleClickedRowProc
		ListBox ListBoxCandidateWaves,listWave=packageFolder:candidatesListBoxWave,mode=4
		ListBox ListBoxCandidateWaves,selRow= 0,selWave=packageFolder:candidatesListBoxSelWave
		ListBox ListBoxCandidateWaves,userdata(ResizeControlsInfo)= A"!!,<7!!#7a!!#B!!!#AFz!!#](Aon\"Qzzzzzzzzzzzzzz!!#o2B4uAezz"
		ListBox ListBoxCandidateWaves,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#u:Du]k<zzzzzzzzzzz"
		ListBox ListBoxCandidateWaves,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		PopupMenu PMAddPositionsWaves,pos={2,190},size={109,20},title="Add positions",proc=PMAddPositionsWavesProc
		PopupMenu PMAddPositionsWaves,mode=0,value= #listPositionsFunc
		PopupMenu PMAddPositionsWaves,userdata(ResizeControlsInfo)= A"!!,=b!!#AM!!#@<!!#<Xz!!#](Aon\"Qzzzzzzzzzzzzzz!!#](Aon\"Qzz"
		PopupMenu PMAddPositionsWaves,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#?(FEDG<zzzzzzzzzzz"
		PopupMenu PMAddPositionsWaves,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		Button BTCreateImageFromMultiplePos,pos={164,189},size={50,20},proc=BTMergePositionsProc,title="OK"
		Button BTCreateImageFromMultiplePos,userdata(ResizeControlsInfo)= A"!!,G4!!#AL!!#>V!!#<Xz!!#o2B4uAezzzzzzzzzzzzzz!!#o2B4uAezz"
		Button BTCreateImageFromMultiplePos,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzz!!#?(FEDG<zzzzzzzzzzz"
		Button BTCreateImageFromMultiplePos,userdata(ResizeControlsInfo) += A"zzz!!#?(FEDG<zzzzzzzzzzzzzz!!!"
		SetWindow kwTopWin,hook(closeControl)=SetupImageViewer_IgorWaves_hook
		SetWindow kwTopWin,hook(ResizeControls)=ResizeControls#ResizeControlsHook
		SetWindow kwTopWin,userdata(ResizeControlsInfo)= A"!!*'\"z!!#B$!!#Abzzzzzzzzzzzzzzzzzzzzz"
		SetWindow kwTopWin,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzzzzzzzzzzzzzzz"
		SetWindow kwTopWin,userdata(ResizeControlsInfo) += A"zzzzzzzzzzzzzzzzzzz!!!"
		
		// if the window was opened previously then restore its size and position
		WC_WindowCoordinatesRestore("positionsWavesViewer")
		SetWindow positionsWavesViewer hook(WindowCoordsHook)=WC_WindowCoordinatesNamedHook
	endif
End

Function PMAddPositionsWavesProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			DFREF packageFolder = root:Packages:Localizer
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			wave/T candidatesListBoxWave = packageFolder:candidatesListBoxWave
			wave candidatesListBoxSelWave = packageFolder:candidatesListBoxSelWave
			
			if (StringMatch(popStr, "* No Positions *") == 1)
				return 0
			endif
			
			wave /Z newPos = GetPositionsWaveReference(popStr)
			
			if (WaveExists(newPos) == 0)
				Abort "Internal error: the positions wave does not exist"
			endif
			
			Redimension /N=(DimSize(candidatesListBoxWave, 0) + 1) candidatesListBoxWave
			Redimension /N=(DimSize(candidatesListBoxSelWave, 0) + 1) candidatesListBoxSelWave
			
			candidatesListBoxWave[DimSize(candidatesListBoxWave, 0) - 1] = popStr
			candidatesListBoxSelWave[DimSize(candidatesListBoxSelWave, 0) - 1] = 0
			
			ControlUpdate /W=positionsWavesViewer ListBoxCandidateWaves
			break
	endswitch

	return 0
End


Function BTMergePositionsProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			DFREF positionsFolder = root:'Localized Positions'
			DFREF packageFolder = root:Packages:Localizer
			
			wave/T waveNames = root:Packages:Localizer:candidatesListBoxWave
			string outputName = ""
			variable nWaves = DimSize(waveNames, 0)
			variable i, j
			variable offset = 0
			variable nPoints, totalNumberOfPoints = 0
			variable lastFrameIndex = 0
			string waveNote
			string originalFilePaths
			string originalCalculationDates
			variable totalCalculationTime
			
			if (nWaves < 2)
				DoWindow /K positionsWavesViewer
				return 0
			endif
			
			outputName = GetNewPositionsWaveName("", positionsFolder, "Enter a name for the merged positions wave")
			if (strlen(outputName) == 0)	// the user canceled the dialog
				return -1
			endif
			
			for (i = 0; i < nWaves; i+=1)
				wave pos = GetPositionsWaveReference(waveNames[i])
				totalNumberOfPoints += DimSize(pos, 0)
			endfor
			
			// did the positions originate from the same localization algorithm?
			variable localizationMethod
			localizationMethod = NumberByKey("LOCALIZATION METHOD", note(pos))
			
			for (i = 1; i < nWaves; i+=1)
				wave pos = GetPositionsWaveReference(waveNames[i])
				if (localizationMethod != NumberByKey("LOCALIZATION METHOD", note(pos)))
					Abort "These positions were acquired using different localization algorithms"
				endif
			endfor 
			
			Make /D/O/N=(totalNumberOfPoints, DimSize(pos, 1)) positionsFolder:$outputName
			
			wave concatenatedPositions = positionsFolder:$outputName
			
			for (i = 0; i < nWaves; i+=1)
				if (i > 0)
					lastFrameIndex += pos[DimSize(pos, 0) - 1][0]
				endif
				wave pos = GetPositionsWaveReference(waveNames[i])
				nPoints = DimSize(pos, 0)
				for (j = 0; j <  nPoints; j+=1)
					concatenatedPositions[offset][] = pos[j][q]
					concatenatedPositions[offset][0] += lastFrameIndex
					offset += 1
				endfor
			endfor
			
			// set up an appropriate wave note that contains the calculation information
			waveNote = note(pos)
			originalFilePaths = StringByKey("ORIGINAL FILE PATH", waveNote)
			originalCalculationDates = StringByKey("CALCULATION DATE", waveNote)
			totalCalculationTime = NumberByKey("CALCULATION DURATION", waveNote)
			for (i = 1; i < nWaves; i+=1)
				originalFilePaths += " and "
				originalFilePaths += StringByKey("ORIGINAL FILE PATH", note(GetPositionsWaveReference(waveNames[i])))
				
				originalCalculationDates += " and "
				originalCalculationDates += StringByKey("CALCULATION DATE", note(GetPositionsWaveReference(waveNames[i])))
				
				totalCalculationTime += NumberByKey("CALCULATION DURATION", note(GetPositionsWaveReference(waveNames[i])))
			endfor
			
			waveNote = ReplaceStringByKey("ORIGINAL FILE PATH", waveNote, originalFilePaths)
			waveNote = ReplaceStringByKey("CALCULATION DATE", waveNote, originalCalculationDates)
			waveNote = ReplaceNumberByKey("CALCULATION DURATION", waveNote, totalCalculationTime)
			
			Note /K concatenatedPositions, waveNote
			
			DoWindow /K positionsWavesViewer
			break
	endswitch

	return 0
End

Function ConsolidateSameEmitters_menu()
	// present a graphical interface to ConsolidateSameEmitters
	string posName
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	NewDataFolder /O root:'Localized Positions'
	
	Variable /G root:Packages:Localizer:V_consolidateMaxDistance
	Variable /G root:Packages:Localizer:V_consolidateMaxBlinking
	Variable /G root:Packages:Localizer:V_consolidateMinTrajLength
	Variable /G root:Packages:Localizer:V_consolidateMaxTrajLength
	
	NVAR maxPositionDifference = root:Packages:Localizer:V_consolidateMaxDistance
	NVAR maxBlinking = root:Packages:Localizer:V_consolidateMaxBlinking
	NVAR consolidateMinTrajLength = root:Packages:Localizer:V_consolidateMinTrajLength
	NVAR consolidateMaxTrajLength = root:Packages:Localizer:V_consolidateMaxTrajLength
	
	maxPositionDifference = (maxPositionDifference == 0) ? 1.0 : maxPositionDifference
	consolidateMinTrajLength = (consolidateMinTrajLength == 0)?  -1 : consolidateMinTrajLength
	consolidateMaxTrajLength = (consolidateMaxTrajLength == 0) ? -1 : consolidateMaxTrajLength
	
	variable maxPositionDifferenceLocal = maxPositionDifference
	variable maxBlinkingLocal = maxBlinking
	variable MinTrajLengthLocal = consolidateMinTrajLength
	variable MaxTrajLengthLocal = consolidateMaxTrajLength
	
	Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
	Prompt maxPositionDifferenceLocal, "Maximum difference in position (in pixels)"
	Prompt maxBlinkingLocal, "Maximum allowed gap due to blinking (in frames)"
	Prompt MinTrajLengthLocal, "Min number of observations for localization error analysis:"
	Prompt MaxTrajLengthLocal, "Max number of observations for localization error analysis:"
	
	DoPrompt "Enter calculation parameters", posName, maxPositionDifferenceLocal, maxBlinkingLocal, MinTrajLengthLocal, MaxTrajLengthLocal
	if (V_flag == 1)	// cancel
		return 0
	endif
	
	maxBlinking = maxBlinkingLocal
	maxPositionDifference = maxPositionDifferenceLocal
	consolidateMinTrajLength = MinTrajLengthLocal
	consolidateMaxTrajLength = MaxTrajLengthLocal
	
	if (MinTrajLengthLocal == -1)
		MinTrajLengthLocal = 2
	endif
	if (MaxTrajLengthLocal == -1)
		MaxTrajLengthLocal = inf
	endif
	
	if ((MinTrajLengthLocal <= 1) || (MaxTrajLengthLocal <= 1) || (MinTrajLengthLocal > MaxTrajLengthLocal))
		Abort "Invalid input for the min and/or max trajectory length for error analysis"
	endif
	
	wave /Z pos = GetPositionsWaveReference(posName)
	if (WaveExists(pos) == 0)
		Abort "The selected wave doesn't seem to exist!"
	endif
	
	// get a wave name for the output positions
	string suggestedOutputName = CleanupName(NameOfWave(pos) + "_cons", 0)
	string outputName = GetNewPositionsWaveName(suggestedOutputName, root:'Localized Positions', "Enter a name for the output wave")
	if (strlen(outputName) == 0)
		return 0
	endif
	
	// do the actual calculation
	ConsolidateSameEmitters(pos, "root:'Localized Positions':" + outputName, maxPositionDifference, maxBlinking, minLocalizationsForError = consolidateMinTrajLength, maxLocalizationsForError = MaxTrajLengthLocal)
	
End

Function RemoveOutlierPositions_Menu()
	DFREF packageFolder = root:Packages:Localizer
	
	NVAR /Z gOutlierRadius = packageFolder:V_outlierRadius
	if (!NVAR_Exists(gOutlierRadius))
		variable /G packageFolder:V_outlierRadius
		NVAR gOutlierRadius = packageFolder:V_outlierRadius
		gOutlierRadius = 10
	endif
	NVAR /Z gNPositionsInRadius = packageFolder:V_nPositionsInRadius
	if (!NVAR_Exists(gNPositionsInRadius))
		variable /G packageFolder:V_nPositionsInRadius
		NVAR gNPositionsInRadius = packageFolder:V_nPositionsInRadius
		gNPositionsInRadius = 25
	endif
	NVAR /Z gOutlierRemovalProcedure = packageFolder:V_OutlierRemovalProcedure
	if (!NVAR_Exists(gOutlierRemovalProcedure))
		variable /G packageFolder:V_OutlierRemovalProcedure
		NVAR gOutlierRemovalProcedure = packageFolder:V_OutlierRemovalProcedure
		gOutlierRemovalProcedure = 1
	endif
	
	string posName
	variable radius = gOutlierRadius
	variable nPosInRadius = gNPositionsInRadius
	variable operation = gOutlierRemovalProcedure
	
	Prompt posName, "Positions wave:", popup, GetPossiblePositionsWaves()
	Prompt radius, "Radius of neighborhood (pixels):"
	Prompt nPosInRadius, "# positions required in neighborhood:"
	Prompt operation "Operation", popup, "Keep only points that meet criteria (=\"nuclei\");Also keep points within radius of nuclei;"
	DoPrompt "Select the positions wave", posName, radius, nPosInRadius, operation
	if (V_flag == 1)	// cancel
		return 0
	endif
	
	if ((radius <= 0) || (nPosInRadius < 1))
		Abort "Expected positive non-zero numbers"
	endif
	
	gOutlierRadius = radius
	gNPositionsInRadius = nPosInRadius
	gOutlierRemovalProcedure = operation
	
	wave /Z positions = GetPositionsWaveReference(posName)
	if (WaveExists(positions) == 0)
		Abort "The selected wave doesn't seem to exist!"
	endif
	// get a wave name for the output positions
	string suggestedOutputName = CleanupName(NameOfWave(positions) + "_filt", 0)
	string outputName = GetNewPositionsWaveName(suggestedOutputName, root:'Localized Positions', "Enter a name for the output wave")
	if (strlen(outputName) == 0)
		return 0
	endif
	
	wave/Z filteredPositions = RemoveOutliers(positions, positions, radius, nPosInRadius)
	if (!WaveExists(filteredPositions))
		return 0	// user abort
	endif
	
	if (operation == 1)
		// remove all non-nuclei
		Duplicate /O filteredPositions, root:'Localized Positions':$outputName
	else
		// only remove non-nuclei and not within radius of nucleus
		wave/Z alsoInRadius = RemoveOutliers(positions, filteredPositions, radius, 1)
		if (!WaveExists(alsoInRadius))
		return 0	// user abort
	endif
		Duplicate /O alsoInRadius, root:'Localized Positions':$outputName
	endif
End

Function CreateRegistrationMap_Menu()
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	variable /G packageFolder:V_RegistrationMaxJump
	string /G packageFolder:S_RegistrationAlgorithm
	variable /G packageFolder:V_nReferencesInLocalMean
	variable /G packageFolder:V_SystematicRegShiftX
	variable /G packageFolder:V_SystematicRegShiftY
	NVAR registrationMaxJump = packageFolder:V_RegistrationMaxJump
	SVAR gRegistrationAlgorithm = packageFolder:S_RegistrationAlgorithm
	NVAR gNReferencesInLocalMean = packageFolder:V_nReferencesInLocalMean
	NVAR gSystematicRegShiftX = packageFolder:V_SystematicRegShiftX
	NVAR gSystematicRegShiftY = packageFolder:V_SystematicRegShiftY
	
	if (registrationMaxJump == 0)
		registrationMaxJump = 10
	endif
	if (gNReferencesInLocalMean == 0)
		gNReferencesInLocalMean = 10
	endif
	
	DoWindow RegistrationMapPanel
	if (V_flag != 0)
		DoWindow /K RegistrationMapPanel
	endif
	
	string possiblePositionsWaves = GetSavedPositionsWaves()
	variable nPositionsWaves = ItemsInList(possiblePositionsWaves)
	variable i
	
	Make /O/N=(nPositionsWaves)/T packageFolder:W_PossiblePosForRegistration1
	Make /O/N=(nPositionsWaves)/T packageFolder:W_PossiblePosForRegistration2
	wave /T W_PossiblePosForRegistration1 = packageFolder:W_PossiblePosForRegistration1
	wave /T W_PossiblePosForRegistration2 = packageFolder:W_PossiblePosForRegistration2
	for (i = 0; i < nPositionsWaves; i+=1)
		W_PossiblePosForRegistration1[i] = StringFromList(i, possiblePositionsWaves)
		W_PossiblePosForRegistration2[i] = StringFromList(i, possiblePositionsWaves)
	endfor
	Make /O/N=(0,2)/T packageFolder:M_RegistrationCombinations
	wave /T M_RegistrationCombinations = packageFolder:M_RegistrationCombinations
	Make /O/T=2 packageFolder:W_RegistrationTitleWave = {"Channel 1", "Channel 2"}
	wave /T W_RegistrationTitleWave = packageFolder:W_RegistrationTitleWave
	
	DoWindow /F RegistrationMapPanel
	if (V_flag != 1)
		NewPanel /K=1/N=RegistrationMapPanel /W=(225,102,1090,423) as "Create registration map"
		SetDrawLayer UserBack
		DrawText 8,24,"Channel 1 (e.g. blue, cyan, green)"
		DrawText 218,24,"Channel 2 (e.g. yellow, orange, red)"
		DrawText 459,24,"Calibration pairs"
		ListBox LBChannel1PositionsWaves,pos={6,26},size={199,176},mode= 2,selRow= 0,listWave=W_PossiblePosForRegistration1
		ListBox LBChannel2PositionsWaves,pos={216,28},size={193,175},mode= 2,selRow= 0,listWave=W_PossiblePosForRegistration2
		ListBox LBCalibrationPairs,pos={458,28},size={396,176},mode= 1,selRow= 0,listWave=M_RegistrationCombinations,titleWave=W_RegistrationTitleWave
		Button BTCreateRegistrationMap,pos={778,294},size={75,20},title="Do it",proc=BTCreateRegistrationMap
		Button BTAddPair,pos={414,93},size={40,20},title=">>",proc=BTAddRegistrationPair
		Button BTRemovePair,pos={414,126},size={40,20},title="<<",proc=BTRemoveRegistrationPair
		TitleBox TBErrorString,pos={464,126},size={50,20}
		TitleBox TBErrorString,pos={9,219},size={146,20},title="",fColor=(65535,0,0)
		SetVariable SVMaxRegistrationJump,pos={562,213},size={296,15},bodyWidth=125,title="Max difference in position between channels (pixels)",value=registrationMaxJump
		PopupMenu BTAlgorithmToUse,pos={576,239},size={279,20},bodyWidth=200,title="Algorithm to use:",proc=PMRegistrationAlgorithmChanged
		PopupMenu BTAlgorithmToUse,mode=1,value= #"\"Local Weighted Mean;Global Poly;\""
		PopupMenu BTAlgorithmToUse,popmatch=gRegistrationAlgorithm
		SetVariable SVReferencesInLocalRegion,pos={594,268},size={264,15},bodyWidth=100,title="Number of references in local region:"
		SetVariable SVReferencesInLocalRegion,value=gNReferencesInLocalMean,limits={6,inf,1}
		SetVariable SVSystematicShiftX,pos={196,270},size={300,15},title="Approx. systematic shift from 1 to 2 along X (pixels)",limits={-inf,inf,0.25},value=gSystematicRegShiftX
		SetVariable SVSystematicShiftY,pos={196,290},size={300,15},title="Approx. systematic shift from 1 to 2 along Y (pixels)",limits={-inf,inf,0.25},value=gSystematicRegShiftY
		SetWindow RegistrationMapPanel, hook(killHook)=RegistrationMapPanel_Hook
		
		ControlInfo BTAlgorithmToUse
		if (StringMatch(S_Value, "Local Weighted Mean") != 1)
			SetVariable SVReferencesInLocalRegion, disable=2
		endif
		
		string error = CheckRegistrationPanelForError()
		TitleBox TBErrorString win=RegistrationMapPanel,title=error
	endif
End

Function RegistrationMapPanel_Hook(s) 
	STRUCT WMWinHookStruct &s
	// this function attaches to a window of choice and appends the X- and Y-coordinates where the user clicked the mouse to a specific wave
	// both a numeric and a text version of the wave is made (for use with listbox controls)
	
	switch (s.eventcode)
		case 2:	// kill
			DFREF packageFolder = root:Packages:Localizer
			wave /T W_PossiblePosForRegistration1 = packageFolder:W_PossiblePosForRegistration1
			wave /T W_PossiblePosForRegistration2 = packageFolder:W_PossiblePosForRegistration2
			wave /T W_RegistrationTitleWave = packageFolder:W_RegistrationTitleWave
			KillWaves /Z W_PossiblePosForRegistration1, W_PossiblePosForRegistration2, W_RegistrationTitleWave
			return 1
			break
		default:
			return 0
			break
	endswitch
End

Function BTAddRegistrationPair(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			wave /T W_PossiblePosForRegistration1 = packageFolder:W_PossiblePosForRegistration1
			wave /T W_PossiblePosForRegistration2 = packageFolder:W_PossiblePosForRegistration2
			wave /T M_RegistrationCombinations = packageFolder:M_RegistrationCombinations
			
			variable selectedRowInChannel1
			variable selectedRowInChannel2
			
			ControlInfo /W=RegistrationMapPanel LBChannel1PositionsWaves
			selectedRowInChannel1 = V_Value
			ControlInfo /W=RegistrationMapPanel LBChannel2PositionsWaves
			selectedRowInChannel2 = V_Value
			
			if ((selectedRowInChannel1 == -1) || (selectedRowInChannel2 == -1) || (DimSize(W_PossiblePosForRegistration1, 0) == 0) || (DimSize(W_PossiblePosForRegistration2, 0) == 0))
				return 1
			endif
			if ((selectedRowInChannel1 >= DimSize(W_PossiblePosForRegistration1, 0)) || (selectedRowInChannel2 >= DimSize(W_PossiblePosForRegistration2, 0)))
				return 1
			endif
			
			string posName1 = W_PossiblePosForRegistration1[selectedRowInChannel1]
			string posName2 = W_PossiblePosForRegistration2[selectedRowInChannel2]
			Redimension /N=(DimSize(M_RegistrationCombinations, 0) + 1, -1) M_RegistrationCombinations
			M_RegistrationCombinations[DimSize(M_RegistrationCombinations, 0) - 1][0] = posName1
			M_RegistrationCombinations[DimSize(M_RegistrationCombinations, 0) - 1][1] = posName2
			
			DeletePoints selectedRowInChannel1, 1, W_PossiblePosForRegistration1
			DeletePoints selectedRowInChannel2, 1, W_PossiblePosForRegistration2
			
			// check for error
			string error = CheckRegistrationPanelForError()
			TitleBox TBErrorString win=RegistrationMapPanel,title=error
			if (strlen(error) > 0)
				Button BTCreateRegistrationMap win=RegistrationMapPanel, disable=2
			else
				Button BTCreateRegistrationMap win=RegistrationMapPanel, disable=0
			endif
			return 1
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTRemoveRegistrationPair(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			wave /T W_PossiblePosForRegistration1 = packageFolder:W_PossiblePosForRegistration1
			wave /T W_PossiblePosForRegistration2 = packageFolder:W_PossiblePosForRegistration2
			wave /T M_RegistrationCombinations = packageFolder:M_RegistrationCombinations
			
			variable selectedCombinationRow
			ControlInfo /W=RegistrationMapPanel LBCalibrationPairs
			selectedCombinationRow = V_Value
			
			if ((selectedCombinationRow == -1) || (DimSize(M_RegistrationCombinations, 0) == 0))
				return 1
			endif
			
			string combinationPosName1 = M_RegistrationCombinations[selectedCombinationRow][0]
			string combinationPosName2 = M_RegistrationCombinations[selectedCombinationRow][1]
			
			Redimension /N=(DimSize(W_PossiblePosForRegistration1, 0) + 1, -1) W_PossiblePosForRegistration1, W_PossiblePosForRegistration2
			W_PossiblePosForRegistration1[DimSize(W_PossiblePosForRegistration1, 0) - 1] = combinationPosName1
			W_PossiblePosForRegistration2[DimSize(W_PossiblePosForRegistration2, 0) - 1] = combinationPosName2
			
			DeletePoints /M=0 selectedCombinationRow, 1, M_RegistrationCombinations
			Redimension /N=(-1, 2) M_RegistrationCombinations	// avoids DeletePoints setting the wave to 1D
			
			// check for error
			string error = CheckRegistrationPanelForError()
			TitleBox TBErrorString win=RegistrationMapPanel,title=error
			if (strlen(error) > 0)
				Button BTCreateRegistrationMap win=RegistrationMapPanel, disable=2
			else
				Button BTCreateRegistrationMap win=RegistrationMapPanel, disable=0
			endif
			
			return 1
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function PMRegistrationAlgorithmChanged(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			string windowName = pa.win
			
			if (StringMatch(popStr, "Local Weighted Mean") == 1)
				SetVariable SVReferencesInLocalRegion,win=$windowName, disable=0
			else
				SetVariable SVReferencesInLocalRegion,win=$windowName, disable=2
			endif
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function /S CheckRegistrationPanelForError()
	DFREF packageFolder = root:Packages:Localizer
	wave /T W_PossiblePosForRegistration1 = packageFolder:W_PossiblePosForRegistration1
	wave /T W_PossiblePosForRegistration2 = packageFolder:W_PossiblePosForRegistration2
	wave /T M_RegistrationCombinations = packageFolder:M_RegistrationCombinations
	
	variable nCombinationsRequested = DimSize(M_RegistrationCombinations, 0)
	if (nCombinationsRequested == 0)
		return "add one or more sets of calibration positions using the '>>' button"
	endif
	
	variable i, j
	variable localizationMethod = -1
	for (i = 0; i < nCombinationsRequested; i+=1)
		for (j = 0; j < 2; j += 1)
			wave thisPos = GetPositionsWaveReference(M_RegistrationCombinations[i][j])
			if (WaveExists(thisPos) == 0)
				return M_RegistrationCombinations[i][j] + " does to appear to be a valid positions wave"
			endif
			
			// all waves should have been localized in the same way
			if (NumType(NumberByKey("LOCALIZATION METHOD", note(thisPos))) == 2)
				// no localization method
				// flag as error
				return M_RegistrationCombinations[i][j] + " does to appear to be a valid positions wave"
			endif
			
			if (localizationMethod == -1)	// first run
				localizationMethod = NumberByKey("LOCALIZATION METHOD", note(thisPos))
			else
				if (localizationMethod != NumberByKey("LOCALIZATION METHOD", note(thisPos)))
					return "the positions appear to have been analyzed using different localization algorithms (e.g. Gaussian Fitting)"
				endif
			endif
		endfor
	endfor
	
	return ""
End

Function BTCreateRegistrationMap(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			wave /T M_RegistrationCombinations = packageFolder:M_RegistrationCombinations
			
			NVAR registrationMaxJump = packageFolder:V_RegistrationMaxJump
			NVAR nReferencesInLocalMean = packageFolder:V_nReferencesInLocalMean
			SVAR gRegistrationAlgorithm = packageFolder:S_RegistrationAlgorithm
			NVAR gSystematicRegShiftX = packageFolder:V_SystematicRegShiftX
			NVAR gSystematicRegShiftY = packageFolder:V_SystematicRegShiftY
			
			ControlInfo /W=RegistrationMapPanel BTAlgorithmToUse
			gRegistrationAlgorithm = S_Value
			
			variable nCombinations = DimSize(M_RegistrationCombinations, 0)
			if (nCombinations == 0)
				return 1
			endif
			
			 // store the positions waves in wave waves
			 Make /O/N=(nCombinations)/WAVE packageFolder:W_Channel1Positions, packageFolder:W_Channel2Positions
			 wave /WAVE W_Channel1Positions = packageFolder:W_Channel1Positions
			 wave /WAVE W_Channel2Positions = packageFolder:W_Channel2Positions
			 variable i, nFramesCol
			 for (i = 0; i < nCombinations; i+=1)
			 	W_Channel1Positions[i] = GetPositionsWaveReference(M_RegistrationCombinations[i][0])
			 	W_Channel2Positions[i] = GetPositionsWaveReference(M_RegistrationCombinations[i][1])
			 	
			 	// make sure that the positions have not been consolidated
			 	wave pos1 = W_Channel1Positions[i]
			 	GetColumnForNFramesPresent(pos1, nFramesCol)
			 	MatrixOP /FREE W_ColSum = sum(col(pos1, nFramesCol))
			 	if (W_ColSum[0] != DimSize(pos1, 0))
			 		Abort NameOfWave(pos1) + " appears to have been consolidated. Only non-consolidated positions can be used"
			 	endif
			 	
			 	wave pos2 = W_Channel2Positions[i]
			 	GetColumnForNFramesPresent(pos2, nFramesCol)
			 	MatrixOP /FREE W_ColSum = sum(col(pos2, nFramesCol))
			 	if (W_ColSum[0] != DimSize(pos2, 0))
			 		Abort NameOfWave(pos2) + " appears to have been consolidated. Only non-consolidated positions can be used"
			 	endif
			 endfor
			 
			 CreateRegistrationMap(W_Channel1Positions, W_Channel2Positions, registrationMaxJump, gRegistrationAlgorithm, nReferencesInLocalMean, gSystematicRegShiftX, gSystematicRegShiftY)
			return 1
			break
	EndSwitch
	
	return 0
End

Function FinishRegistrationMap(M_RegistrationMap)
	wave M_RegistrationMap
	
	DFREF packageFolder = root:Packages:Localizer
	
	NVAR gRegistrationMaxJump = packageFolder:V_RegistrationMaxJump
	
	wave /T M_RegistrationCombinations = packageFolder:M_RegistrationCombinations
	wave /WAVE W_Channel1Positions = packageFolder:W_Channel1Positions
	wave /WAVE W_Channel2Positions = packageFolder:W_Channel2Positions
	
	// get an output wave name
	NewDataFolder /O root:'Registration Maps'
	DFREF outputDF = root:'Registration Maps'
	string outputName = GetNewPositionsWaveName("", outputDF, "Select an output name")
	if (strlen(outputName) == 0)
		return 1	// user cancel
	endif

	Duplicate /O M_RegistrationMap, outputDF:$outputName
	 wave output = outputDF:$outputName
	 
	 // add a wave note
	 variable i, nCombinations = DimSize(M_RegistrationCombinations, 0)
	 string listOfPos1 = "", listOfPos2 = ""
	 for (i = 0; i < nCombinations; i+=1)
	 	listOfPos1 += M_RegistrationCombinations[i][0] + ","
	 	listOfPos2 += M_RegistrationCombinations[i][1] + ","
	 endfor
	 // trim added comma
	 listOfPos1 = listOfPos1[0, strlen(listOfPos1) - 2]
	 listOfPos2 = listOfPos2[0, strlen(listOfPos2) - 2]
	 string waveNote = "TYPE:RegistrationMap;MAX SHIFT:" + num2str(gRegistrationMaxJump) + ";"
	 waveNote += "POSITIONS1:" + listOfPos1 + ";"
	 waveNote += "POSITIONS2:" + listOfPos2 + ";"
	 
	 Note /K output, waveNote
	 
	 DoWindow /K RegistrationMapPanel
	 
	 // provide a graph for the user to check whether the registration makes sense
	 Make /O/N=0/D packageFolder:W_Channel1X, packageFolder:W_Channel1Y
	 Make /O/N=0/D packageFolder:W_Channel2X, packageFolder:W_Channel2Y
	 Make /O/N=0/D packageFolder:W_CorrectedChannel1X, packageFolder:W_CorrectedChannel1Y
	 wave W_Channel1X = packageFolder:W_Channel1X; wave W_Channel1Y = packageFolder:W_Channel1Y
	 wave W_Channel2X = packageFolder:W_Channel2X; wave W_Channel2Y = packageFolder:W_Channel2Y
	 wave W_CorrectedChannel1X = packageFolder:W_CorrectedChannel1X; wave W_CorrectedChannel1Y = packageFolder:W_CorrectedChannel1Y
	 
	  variable xCol, yCol, zCol, offset
	 for (i = 0; i < nCombinations; i+=1)
	 	wave pos1 = W_Channel1Positions[i]
	 	wave pos2 = W_Channel2Positions[i]
	 	if (DimSize(output, 1) == 17)
	 		wave Corrected = ApplyRegistrationMap(pos1, output)
	 	else
	 		wave Corrected = ApplyRegistrationMap_Poly(pos1, output)
	 	endif
	 	
	 	offset = DimSize(W_Channel1X, 0)
	 	Redimension /N=(DimSize(W_Channel1X, 0) + DimSize(pos1, 0)) W_Channel1X, W_Channel1Y
	 	GetColumnsForEmitterPositions(pos1, xCol, yCol, zCol)
	 	W_Channel1X[offset, ] = pos1[p - offset][xCol]
	 	W_Channel1Y[offset, ] = pos1[p - offset][yCol]
	 	
	 	offset = DimSize(W_Channel2X, 0)
	 	Redimension /N=(DimSize(W_Channel2X, 0) + DimSize(pos2, 0)) W_Channel2X, W_Channel2Y
	 	GetColumnsForEmitterPositions(pos2, xCol, yCol, zCol)
	 	W_Channel2X[offset, ] = pos2[p - offset][xCol]
	 	W_Channel2Y[offset, ] = pos2[p - offset][yCol]
	 	
	 	offset = DimSize(W_CorrectedChannel1X, 0)
	 	Redimension /N=(DimSize(W_CorrectedChannel1X, 0) + DimSize(Corrected, 0)) W_CorrectedChannel1X, W_CorrectedChannel1Y
		GetColumnsForEmitterPositions(Corrected, xCol, yCol, zCol)
		W_CorrectedChannel1X[offset, ] = Corrected[p - offset][xCol]
		W_CorrectedChannel1Y[offset, ] = Corrected[p - offset][yCol]
	 endfor
	 
	// and plot the result
	DoWindow /F CreateRegistrationResultViewer
	if (V_flag != 1)
		Display /K=1/N=CreateRegistrationResultViewer W_Channel1Y vs W_Channel1X as "Registration result"
		AppendToGraph W_Channel2Y vs W_Channel2X
		AppendToGraph W_CorrectedChannel1Y vs W_CorrectedChannel1X
		ModifyGraph mode=3,rgb(W_Channel2Y)=(0,0,65535),marker(W_CorrectedChannel1Y)=8,rgb(W_CorrectedChannel1Y)=(1,26214,0)
		Legend/C/N=text0/J "\\s(W_Channel1Y) Channel 1\r\\s(W_Channel2Y) Channel 2\r\\s(W_CorrectedChannel1Y) Ch1/Corrected"
		Label left, "Distance (pixels)"
		Label bottom, "Distance (pixels)"
	endif
	
	// also show the user a visualization of the registration shift
	variable xStart, xEnd, yStart, yEnd
	// arbitrarily fetch the image size from the first positions in channel 1 (should be the same for all but we don't check)
	wave pos = W_Channel1Positions[0]
	xStart = 0
	yStart = 0
	xEnd = NumberByKey("X SIZE", note(pos)) - 1
	if (NumType(xEnd) == 2)
		xEnd = 511 // guess
	endif
	yEnd = NumberByKey("Y SIZE", note(pos)) - 1
	if (NumType(yEnd) == 2)
		yEnd = 511 // guess
	endif
	VisualizeRegistrationMap(M_RegistrationMap, xStart, xEnd, yStart, yEnd)

	return 0
End

Function CreateRegistrationMap(W_Channel1Positions, W_Channel2Positions, maxShiftDistance, algorithmToUse, nReferencesInLocalMean, systematicShiftX, systematicShiftY)
	wave /WAVE W_Channel1Positions, W_Channel2Positions
	variable maxShiftDistance
	string algorithmToUse
	variable nReferencesInLocalMean, systematicShiftX, systematicShiftY
	
	// given sets of positions measured in different detector channels, create a registration map
	// that will contain the shifts needed to overlap the emitters in pos1 onto those in pos2
	
	// The emitters need to be close enough together to qualify as the same particle,
	// where 'close enough' is determined by maxShiftDistance
	
	// The high level overview is at follows: the x and y positions in every combination
	// are combined using GroupEmitters() to yield x-y pairs that appear to belong together.
	
	// see http://www.cs.wright.edu/~agoshtas/IVC88.pdf
	
	if (DimSize(W_Channel1Positions, 0) != DimSize(W_Channel2Positions, 0))
		Abort "Need to have same number of positions waves for both channels in CreateRegistrationMap"
	endif
	variable nPositionsWaves = DimSize(W_Channel1Positions, 0)
	
	Make /N=(0, 4)/D/FREE M_CombinedEmitterPositions
	
	variable n, localizationMethod
	for (n = 0; n < nPositionsWaves; n+=1)
		wave pos1 = W_Channel1Positions[n]
		wave pos2 = W_Channel2Positions[n]
		
		// verify that the positions aren't empty
		if ((DimSize(pos1, 0) == 0) || (DimSize(pos2, 0) == 0))
			Abort "One or more of the positions waves does not contain any positions"
		endif
		
		// check if all emitters were localized in the same way
		if (n == 0)
			localizationMethod = NumberByKey("LOCALIZATION METHOD", note(pos1))
			if (NumType(localizationMethod) == 2)	// NaN
				Abort "The input wave named " + NameOfWave(pos1) + " does not appear to be a regular set of positions"
			endif
		else
			if (localizationMethod != NumberByKey("LOCALIZATION METHOD", note(pos1)))
				Abort "The positions appear to have been localized using different algorithms. Please use the same algorithm for all the images"
			endif
		endif
		
		if (localizationMethod != NumberByKey("LOCALIZATION METHOD", note(pos2)))
			Abort "The positions appear to have been localized using different algorithms. Please use the same algorithm for all the images"
		endif
		
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
		
		// try to estimate whether there is a systematic shift between the particles
		// but only if it is significant
		GuessSystematicShiftInPositions(extractedPositions1, extractedPositions2, maxShiftDistance, systematicShiftX, systematicShiftY)
		systematicShiftX = (abs(systematicShiftX) > 0.25) ? systematicShiftX : 0
		systematicShiftY = (abs(systematicShiftY) > 0.25) ? systematicShiftY : 0
		
		// subract the systematic shift from the positions in channel2. It will be added back in
		// after the emitters have been combined
		extractedPositions2[][xCol] -= systematicShiftX
		extractedPositions2[][yCol] -= systematicShiftY
		
		// the grouping functions that already exist are designed for consolidation or particle tracking
		// so make a new positions wave, in which pos2 is appended to pos1
		variable nEmittersIn1 = DimSize(extractedPositions1, 0)
		variable nEmittersIn2 = DimSize(extractedPositions2, 0)
		Duplicate /O /FREE extractedPositions1, M_AppendedEmitters	// duplicate takes care of wave note
		M_AppendedEmitters[][0] = 0
		Redimension /N=(nEmittersIn1 + nEmittersIn2, -1) M_AppendedEmitters
		M_AppendedEmitters[nEmittersIn1, ] = extractedPositions2[p - nEmittersIn1][q]
		M_AppendedEmitters[nEmittersIn1, ][0] = 1
		
		wave /WAVE W_GroupedEmitters = GroupEmitters(M_AppendedEmitters, CombinedPosEstimator_LastPos, maxShiftDistance, 0)
		
		// In W_GroupedEmitters, we're only interested in those subwaves that contain 2 rows (i.e. for which a match was successfully found). We selectively
		// extract those and convert them to a format suitable for passing to CreateRegistration().
		variable i, nMatchedEmitters = 0
		for (i = 0; i < DimSize(W_GroupedEmitters, 0); i+=1)
			wave thisGroup = W_GroupedEmitters[i]
			if (DimSize(thisGroup, 0) == 2)
				nMatchedEmitters += 1
			endif
		endfor
		
		variable offset = DimSize(M_CombinedEmitterPositions, 0)
		Redimension /N=(DimSize(M_CombinedEmitterPositions, 0) + nMatchedEmitters, -1) M_CombinedEmitterPositions
		
		for (i = 0; i < DimSize(W_GroupedEmitters, 0); i+=1)
			wave thisGroup = W_GroupedEmitters[i]
			if (DimSize(thisGroup, 0) == 2)
				M_CombinedEmitterPositions[offset][0] = thisGroup[0][xCol]
				M_CombinedEmitterPositions[offset][1] = thisGroup[0][yCol]
				M_CombinedEmitterPositions[offset][2] = thisGroup[1][xCol]
				M_CombinedEmitterPositions[offset][3] = thisGroup[1][yCol]
				
				// add the systematic offset for channel 2 back in
				M_CombinedEmitterPositions[offset][2] += systematicShiftX
				M_CombinedEmitterPositions[offset][3] += systematicShiftY
				
				offset += 1
			endif
		endfor
	endfor
	
	if (DimSize(M_CombinedEmitterPositions, 0) < 6)
		Abort "Insufficient number of reference points in the positions"
	endif
	
	// allow the user to clean up the combined positions
	CleanCombinedRegistrationPosDlg(M_CombinedEmitterPositions)
	return 0
End

Function CleanCombinedRegistrationPosDlg(M_CombinedPositions)
	wave M_CombinedPositions
	
	DFREF packageFolder = root:Packages:Localizer
	
	variable nMatchingPositions = DimSize(M_CombinedPositions, 0)
	
	Make /O/N=(3 * nMatchingPositions) packageFolder:W_RegistrationMatchX, packageFolder:W_RegistrationMatchY
	Make /O/N=(3 * nMatchingPositions, 3)/W/U packageFolder:M_RegistrationMatch_ColorsAngle, packageFolder:M_RegistrationMatch_ColorsDist
	wave W_RegistrationMatchX = packageFolder:W_RegistrationMatchX
	wave W_RegistrationMatchY = packageFolder:W_RegistrationMatchY
	wave M_RegistrationMatch_ColorsAngle = packageFolder:M_RegistrationMatch_ColorsAngle
	wave M_RegistrationMatch_ColorsDist = packageFolder:M_RegistrationMatch_ColorsDist
	
	variable i
	Make /FREE/N=(nMatchingPositions) /D W_Angles, W_Distances
	for (i = 0; i < nMatchingPositions; i+=1)
		// x coordinates
		W_RegistrationMatchX[3 * i] = M_CombinedPositions[i][0]
		W_RegistrationMatchX[3 * i + 1] = M_CombinedPositions[i][2]
		W_RegistrationMatchX[3 * i + 2] = NaN
		
		// y coordinates
		W_RegistrationMatchY[3 * i] = M_CombinedPositions[i][1]
		W_RegistrationMatchY[3 * i + 1] = M_CombinedPositions[i][3]
		W_RegistrationMatchY[3 * i + 2] = NaN
		
		// color depending on angle between them
		variable /C asComplex = r2polar(cmplx(M_CombinedPositions[i][2] - M_CombinedPositions[i][0], M_CombinedPositions[i][3] - M_CombinedPositions[i][1]))
		W_Angles[i] = imag(asComplex)
		W_Distances[i] = real(asComplex)
	endfor
	ColorTab2Wave Rainbow
	wave M_Colors
	variable nColors = DimSize(M_Colors, 0)
	WaveStats /Q/M=1 W_Angles
	for (i = 0; i < nMatchingPositions; i+=1)
		M_RegistrationMatch_ColorsAngle[3 * i, 3 * i + 1][] = M_Colors[floor((W_Angles[i] - V_min) / (V_max - V_min + 1e-9) * nColors)][q]
	endfor
	WaveStats /Q/M=1 W_Distances
	for (i = 0; i < nMatchingPositions; i+=1)
		M_RegistrationMatch_ColorsDist[3 * i, 3 * i + 1][] = M_Colors[floor((W_Distances[i] - V_min) / (V_max - V_min + 1e-9) * nColors)][q]
	endfor
	KillWaves /Z M_Colors
	
	// delete the graph or panel in case they already exist
	DoWindow /K CleanRegistrationGraph
	DoWindow /K CleanRegistrationPanel
	
	// display a graph showing the matches
	NewPanel /K=1/N=CleanRegistrationPanel /W=(150,50,1148,650)
	DefineGuide GraphLeft={FL,15},GraphTop={FT,15},GraphRight={FR,-200},GraphBottom={FB,-15}
	DefineGuide PanelLeft={GraphRight,15},PanelTop={GraphTop,0},PanelRight={FR,-15},PanelBottom={GraphBottom,0}
	Display /HOST=CleanRegistrationPanel /FG=(GraphLeft,GraphTop,GraphRight,GraphBottom) /N=CleanRegistrationGraph W_RegistrationMatchY vs W_RegistrationMatchX
	ModifyGraph /W=CleanRegistrationPanel#CleanRegistrationGraph zColor(W_RegistrationMatchY)={M_RegistrationMatch_ColorsAngle,*,*,directRGB,0}
	ModifyGraph /W=CleanRegistrationPanel#CleanRegistrationGraph mode=4,marker=8//,gbRGB=(128,128,128)
	SetWindow CleanRegistrationPanel, hook(cleanupHook)=CleanRegistration_WindowHook
	
	// create a panel with instructions
	NewPanel/N=CleanRegistrationControls /HOST=CleanRegistrationPanel /FG=(PanelLeft,PanelTop,PanelRight,PanelBottom)
	SetDrawLayer UserBack
	DrawText /W=CleanRegistrationPanel#CleanRegistrationControls 8,90,"Remove spurious references\rby clicking near\rone of the emitters,\rthen press continue"
	Button BTContinue,win=CleanRegistrationPanel#CleanRegistrationControls,pos={55,214},size={75,20},title="Continue",proc=BTContinueRegistrationCleanProc,fColor=(1,34817,52428)
	Button BTShowPrevious,win=CleanRegistrationPanel#CleanRegistrationControls,pos={31,116},size={125,20},title="Show previous",disable=2,proc=BTShowPreviousRegistrationCorr
	Button BTHidePrevious,win=CleanRegistrationPanel#CleanRegistrationControls,pos={31,145},size={125,20},title="Hide previous",disable=2,proc=BTHidePreviousRegistrationCorr
	Button BTUsePrevious,win=CleanRegistrationPanel#CleanRegistrationControls,pos={31,174},size={125,20},title="Use previous",disable=2,proc=BTUsePreviousRegistrationCorr
	CheckBox CBColorShowsAngle,win=CleanRegistrationPanel#CleanRegistrationControls,pos={13,249},size={104,14},title="Color shows angle"
	CheckBox CBColorShowsAngle,win=CleanRegistrationPanel#CleanRegistrationControls,value= 1,mode=1,proc=RBRegistrationColorProc
	CheckBox CBColorShowsDistance,win=CleanRegistrationPanel#CleanRegistrationControls,pos={13,268},size={104,14},title="Color shows distance"
	CheckBox CBColorShowsDistance,win=CleanRegistrationPanel#CleanRegistrationControls,value= 0,mode=1,proc=RBRegistrationColorProc
	
	// allow the user to show or use previous corrections if needed
	wave /Z W_PreviousRegistrationMatchX = packageFolder:W_PreviousRegistrationMatchX
	wave /Z W_PreviousRegistrationMatchY = packageFolder:W_PreviousRegistrationMatchY
	if (WaveExists(W_PreviousRegistrationMatchX) && WaveExists(W_PreviousRegistrationMatchY))
		Button BTShowPrevious win=CleanRegistrationPanel#CleanRegistrationControls, disable=0
		Button BTUsePrevious win=CleanRegistrationPanel#CleanRegistrationControls, disable=0
	endif
End

Function RBRegistrationColorProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			DFREF packageFolder = root:Packages:Localizer
			wave M_RegistrationMatch_ColorsAngle = packageFolder:M_RegistrationMatch_ColorsAngle
			wave M_RegistrationMatch_ColorsDist = packageFolder:M_RegistrationMatch_ColorsDist
			Variable checked = cba.checked
			CheckBox CBColorShowsAngle,win=CleanRegistrationPanel#CleanRegistrationControls,value=0
			CheckBox CBColorShowsDistance,win=CleanRegistrationPanel#CleanRegistrationControls,value=0
			
			StrSwitch(cba.ctrlName)
				case "CBColorShowsAngle":
					CheckBox CBColorShowsAngle,win=CleanRegistrationPanel#CleanRegistrationControls,value=1
					ModifyGraph /W=CleanRegistrationPanel#CleanRegistrationGraph zColor(W_RegistrationMatchY)={M_RegistrationMatch_ColorsAngle,*,*,directRGB,0}
					break
				case "CBColorShowsDistance":
					CheckBox CBColorShowsDistance,win=CleanRegistrationPanel#CleanRegistrationControls,value=1
					ModifyGraph /W=CleanRegistrationPanel#CleanRegistrationGraph zColor(W_RegistrationMatchY)={M_RegistrationMatch_ColorsDist,*,*,directRGB,0}
					break
				default:
					Abort "Unknown control name in CBRegistrationColorProc"
					break
			EndSwitch
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function BTContinueRegistrationCleanProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			SVAR gRegistrationAlgorithm = packageFolder:S_RegistrationAlgorithm
			NVAR gNReferencesInLocalMean = packageFolder:V_nReferencesInLocalMean
			
			// save the corrected results for potential reuse
			wave W_RegistrationMatchX = packageFolder:W_RegistrationMatchX
			wave W_RegistrationMatchY = packageFolder:W_RegistrationMatchY
			Duplicate /O W_RegistrationMatchX, packageFolder:W_PreviousRegistrationMatchX
			Duplicate /O W_RegistrationMatchY, packageFolder:W_PreviousRegistrationMatchY
			
			// now piece together those matches that remain
			assert(mod(DimSize(W_RegistrationMatchY, 0), 3) == 0)
			variable nMatchesRemaining = DimSize(W_RegistrationMatchY, 0) / 3, i
			Make /N=(nMatchesRemaining, 4) /D/O packageFolder:M_RemainingCombinedPositions
			wave M_RemainingCombinedPositions = packageFolder:M_RemainingCombinedPositions
			for (i = 0; i < nMatchesRemaining; i+=1)
				M_RemainingCombinedPositions[i][0] = W_RegistrationMatchX[3 * i]
				M_RemainingCombinedPositions[i][1] = W_RegistrationMatchY[3 * i]
				M_RemainingCombinedPositions[i][2] = W_RegistrationMatchX[3 * i + 1]
				M_RemainingCombinedPositions[i][3] = W_RegistrationMatchY[3 * i + 1]
			endfor
			
			Duplicate /O M_RemainingCombinedPositions, packageFolder:M_CleanedCombinedPositionsSav
	
			if (StringMatch(gRegistrationAlgorithm, "Local Weighted Mean") == 1)
				wave M_RegistrationMap = CreateRegistration(M_RemainingCombinedPositions, gNReferencesInLocalMean)
			else
				wave M_RegistrationMap = CreateRegistration_Poly(M_RemainingCombinedPositions)
			endif
			
			FinishRegistrationMap(M_RegistrationMap)
			
			DoWindow /K CleanRegistrationPanel
			
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function BTShowPreviousRegistrationCorr(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			string baseWindowName = StringFromList(0, ba.win, "#")
			string graphName = baseWindowName + "#CleanRegistrationGraph"
			
			wave /Z W_PreviousRegistrationMatchX = packageFolder:W_PreviousRegistrationMatchX
			wave /Z W_PreviousRegistrationMatchY = packageFolder:W_PreviousRegistrationMatchY
			assert(WaveExists(W_PreviousRegistrationMatchX) && WaveExists(W_PreviousRegistrationMatchY))
			
			AppendToGraph /W=$graphName W_PreviousRegistrationMatchY vs W_PreviousRegistrationMatchX
			ModifyGraph /W=$graphName mode(W_PreviousRegistrationMatchY)=4,marker(W_PreviousRegistrationMatchY)=6,rgb(W_PreviousRegistrationMatchY)=(0,0,65535)
			Button BTShowPrevious win=$baseWindowName#CleanRegistrationControls, disable=2
			Button BTHidePrevious win=$baseWindowName#CleanRegistrationControls, disable=0
			break
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function BTHidePreviousRegistrationCorr(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			string baseWindowName = StringFromList(0, ba.win, "#")
			string graphName = baseWindowName + "#CleanRegistrationGraph"
			
			wave /Z W_PreviousRegistrationMatchX = packageFolder:W_PreviousRegistrationMatchX
			wave /Z W_PreviousRegistrationMatchY = packageFolder:W_PreviousRegistrationMatchY
			assert(WaveExists(W_PreviousRegistrationMatchX) && WaveExists(W_PreviousRegistrationMatchY))
			
			RemoveFromGraph /W=$graphName W_PreviousRegistrationMatchY
			Button BTShowPrevious win=$baseWindowName#CleanRegistrationControls, disable=0
			Button BTHidePrevious win=$baseWindowName#CleanRegistrationControls, disable=2
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function BTUsePreviousRegistrationCorr(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			wave /Z W_PreviousRegistrationMatchX = packageFolder:W_PreviousRegistrationMatchX
			wave /Z W_PreviousRegistrationMatchY = packageFolder:W_PreviousRegistrationMatchY
			assert(WaveExists(W_PreviousRegistrationMatchX) && WaveExists(W_PreviousRegistrationMatchY))
			
			wave W_RegistrationMatchX = packageFolder:W_RegistrationMatchX
			wave W_RegistrationMatchY = packageFolder:W_RegistrationMatchY
			
			Duplicate /O W_PreviousRegistrationMatchX, W_RegistrationMatchX
			Duplicate /O W_PreviousRegistrationMatchY, W_RegistrationMatchY
			
			BTContinueRegistrationCleanProc(ba)
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function CleanRegistration_WindowHook(s)
	STRUCT WMWinHookStruct &s
	
	switch (s.eventcode)
		case 5:	// mouseup
			
			DFREF packageFolder = root:Packages:Localizer
			string baseWindowName = StringFromList(0, s.winName, "#")
			string graphName = baseWindowName + "#CleanRegistrationGraph"
			
			// get the coordinates of the user click in terms of the graph axes
			variable minX, minY, maxX, maxY, xVal, yVal
			GetAxis /Q/W=$graphName Bottom
			minX = V_min
			maxX = V_max
			GetAxis /Q/W=$graphName Left
			minY = V_min
			maxY = V_max
			
			xVal = AxisValFromPixel(graphName, "Bottom", s.mouseLoc.h)
			yVal = AxisValFromPixel(graphName, "Left", s.mouseLoc.v)
			
			if ((xVal < minX) || (xVal > maxX) || (yVal < minY) || (yVal > maxY))
				return 1
			endif
			
			wave W_RegistrationMatchX = packageFolder:W_RegistrationMatchX
			wave W_RegistrationMatchY = packageFolder:W_RegistrationMatchY
			wave M_RegistrationMatch_ColorsAngle = packageFolder:M_RegistrationMatch_ColorsAngle
			wave M_RegistrationMatch_ColorsDist = packageFolder:M_RegistrationMatch_ColorsDist
			
			// now look for the point that is closest to that clicked by the user
			variable closestDistance = inf, closestIndex
			variable distance, i
			for (i = 0; i < DimSize(W_RegistrationMatchX, 0); i+=1)
				distance = sqrt((W_RegistrationMatchX[i] - xVal)^2 + (W_RegistrationMatchY[i] - yVal)^2)
				if (distance < closestDistance)
					closestDistance = distance
					closestIndex = i
				endif
			endfor
			
			// remove the point if it is with 10 pixels of the mouseclick (arbitrary limit)
			if (closestDistance < 10)
				DeletePoints floor(closestIndex / 3) * 3, 3, W_RegistrationMatchX, W_RegistrationMatchY, M_RegistrationMatch_ColorsAngle, M_RegistrationMatch_ColorsDist
			endif
			
			DoUpdate
			
			break
		default:
			return 0
	endswitch
	
	return 0
End

Function ApplyRegistrationMap_Menu()
	DFREF mapFolder = root:'Registration Maps'
	DFREF positionsFolder = root:'Localized Positions'
	
	// start by getting a list of positions waves and registration maps
	string positionsWaves = GetPossiblePositionsWaves()
	string registrationMaps = GetPossibleRegistrationMaps()
	
	variable nPositionsWaves = ItemsInList(positionsWaves)
	variable nRegistrationMaps = ItemsInList(registrationMaps)
	
	if (nPositionsWaves == 0)
		Abort "No positions waves seem to be present"
	endif
	if (nRegistrationMaps == 0)
		Abort "No registration maps seem to be present"
	endif
	
	string posName, registrationMapName
	Prompt posName, "Positions to correct (from the first channel!):", popup, positionsWaves
	Prompt registrationMapName, "Registration map to use", popup, registrationMaps
	DoPrompt "Select calculation parameters", posName, registrationMapName
	if (V_flag == 1)
		return 0
	endif
	
	wave posWave = GetPositionsWaveReference(posName)
	wave mapWave = mapFolder:$registrationMapName
	
	string newOutputName = GetNewPositionsWaveName(posName + "_mapped", positionsFolder, "Enter the name of the output wave")
	if (strlen(newOutputName) == 0)
		return 0
	endif
	
	variable ccdXSize, ccdYSize, pixelSize
	GetCCDDimensionsFromPositions(posWave, ccdXSize, ccdYSize, pixelSize)
	
	// algorithm to use is determined by the number of columns in the registration map
	if (DimSize(mapWave, 1) == 17)
		wave M_Mapped = ApplyRegistrationMap(posWave, mapWave)
	else
		wave M_Mapped = ApplyRegistrationMap_Poly(posWave, mapWave)
	endif
	
	// remove all positions for which no shift could be determined, or positions that were shifted excessively
	variable xCol, yCol, zCol
	GetColumnsForEmitterPositions(M_Mapped, xCol, yCol, zCol)
	
	variable minX = -ccdXSize / 2
	variable minY = -ccdYSize / 2
	variable maxX = 1.5 * ccdXSize
	variable maxY = 1.5 * ccdYSize
	
	variable nCorrectedPositions = DimSize(M_Mapped, 0)
	
	Duplicate /FREE M_Mapped, M_CleanedMappedPositions
	
	variable offset = 0, i
	for (i = 0; i < nCorrectedPositions; i+=1)
		if ((NumType(M_Mapped[i][xCol]) == 2) || !Within(M_CleanedMappedPositions[i][xCol], minX, maxX) || !Within(M_CleanedMappedPositions[i][yCol], minY, maxY))
			continue
		endif
		M_CleanedMappedPositions[offset][] = M_Mapped[i][q]
		offset += 1
	endfor
	
	Redimension /N=(offset, -1) M_CleanedMappedPositions
	Duplicate /O M_CleanedMappedPositions, positionsFolder:$newOutputName
End

Function VisualizeRegistrationMap(RegistrationMap, xStart, xEnd, yStart, yEnd)
	wave RegistrationMap
	variable xStart, xEnd, yStart, yEnd	// must be integers - not checked for errors
	
	DFREF packageFolder = root:Packages:Localizer
	
	variable nPixelsX = round(xEnd - xStart + 1)
	variable nPixelsY = round(yEnd - yStart + 1)
	
	// make a set of fake positions to correct, with one position at every pixel
	Make /FREE/D/N=(nPixelsX * nPixelsY, 12) M_FakePositions
	Note M_FakePositions, "LOCALIZATION METHOD:" + num2str(LOCALIZATION_GAUSS_FITTING)
	variable i, j, offset = 0
	for (i = xStart; i <= xEnd; i+=1)
		for (j = yStart; j <= yEnd; j+=1)
			M_FakePositions[offset][3] = i
			M_FakePositions[offset][4] = j
			offset += 1
		endfor
	endfor
	
	if (DimSize(RegistrationMap, 1) == 17)
		wave CorrectedPositions = ApplyRegistrationMap(M_FakePositions, RegistrationMap)
	else
		wave CorrectedPositions = ApplyRegistrationMap_Poly(M_FakePositions, RegistrationMap)
	endif
	
	// create output images to mark the shift
	Make /O/N=(nPixelsX, nPixelsY) /D packageFolder:M_CorrectionX, packageFolder:M_CorrectionY
	wave M_CorrectionX = packageFolder:M_CorrectionX
	wave M_CorrectionY = packageFolder:M_CorrectionY
	SetScale /P x, xStart, 1, M_CorrectionX, M_CorrectionY
	SetScale /P y, yStart, 1, M_CorrectionX, M_CorrectionY
	
	variable dx, dy
	for (i = 0; i < DimSize(M_FakePositions, 0); i+=1)
		dx = CorrectedPositions[i][3] - M_FakePositions[i][3]
		dy = CorrectedPositions[i][4] - M_FakePositions[i][4]
		
		M_CorrectionX[M_FakePositions[i][3]][M_FakePositions[i][4]] = dx
		M_CorrectionY[M_FakePositions[i][3]][M_FakePositions[i][4]] = dy
	endfor
	
	DoWindow /F RegistrationVisualizationX
	if (V_flag != 1)
		Display /K=1/N=RegistrationVisualizationX as "Registration shift in X"
		AppendImage /W=$S_Name M_CorrectionX
	endif
	
	DoWindow /F RegistrationVisualizationY
	if (V_flag != 1)
		Display /K=1/N=RegistrationVisualizationY as "Registration shift in Y"
		AppendImage /W=$S_Name M_CorrectionY
		AutoPositionWindow /R=RegistrationVisualizationX RegistrationVisualizationY
	endif
End

Function Create3DCalibration_Menu()
	variable method
	Prompt method, "Format:", popup, "Separate positions per focus depth;One set of positions for all depths;"
	DoPrompt "Select data format", method
	if (V_flag == 1)
		return 0
	endif
	
	switch (method)
		case 1:	// separate per depth
			Make3DCalibration_ManyPos()
			break
		case 2:
			Make3DCalibration_SinglePos()
			break
		default:
			Abort "Unknown method in Create3DCalibration_Menu()"
			break
	endswitch
End

Function Make3DCalibration_ManyPos()
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	// make list and selection waves for the listboxes
	Make /N=(0) /O /T packageFolder:W_EligibleAstigmPosListWave
	Make /N=(0) /O /D packageFolder:W_EligibleAstigmPosSelWave
	Make /N=(0, 2) /O /T packageFolder:M_Selected3DPairsListWave
	Make /N=(0, 2) /O /D packageFolder:M_Selected3DPairsSelWave
	wave/T W_EligibleAstigmPosListWave = packageFolder:W_EligibleAstigmPosListWave
	wave W_EligibleAstigmPosSelWave = packageFolder:W_EligibleAstigmPosSelWave
	wave/T M_Selected3DPairsListWave = packageFolder:M_Selected3DPairsListWave
	wave M_Selected3DPairsSelWave = packageFolder:M_Selected3DPairsSelWave
	
	// title waves
	Make /N=2 /O /T packageFolder:W_Selected3DPairsTitleWave = {"Selected positions", "z position (nm)"}
	Make /N=1 /O /T packageFolder:W_EligibleAstigmPosTitleWave = {"Eligible positions"}
	wave/T W_Selected3DPairsTitleWave = packageFolder:W_Selected3DPairsTitleWave
	wave/T W_EligibleAstigmPosTitleWave = packageFolder:W_EligibleAstigmPosTitleWave
	
	// only positions fitted using astigmatism are eligible
	string eligiblePositionsNames = GetGaussAstigPositionsWaves()
	
	Redimension /N=(ItemsInList(eligiblePositionsNames)) W_EligibleAstigmPosListWave, W_EligibleAstigmPosSelWave
	variable i
	for (i = 0; i < DimSize(W_EligibleAstigmPosListWave, 0); i+=1)
		W_EligibleAstigmPosListWave[i] = StringFromList(i, eligiblePositionsNames)
	endfor
	
	DoWindow /F Make3DCalibrationPanel
	if (V_flag == 0)
		NewPanel /K=1/N=Make3DCalibrationPanel/W=(56,168,728,434)
		SetDrawLayer UserBack
		SetDrawEnv fstyle= 2
		DrawText 10,22,"Select positions fitted on calibration data using ellipsoidal 2D Gauss (astigmatism)"
		ListBox LBEligible3DCalibrationPos,pos={9,29},size={255,181},listWave=W_EligibleAstigmPosListWave,selWave=W_EligibleAstigmPosSelWave,titleWave=W_EligibleAstigmPosTitleWave,mode=4
		ListBox LBSelected3DCalibrationData,pos={335,29},size={327,181},listWave=M_Selected3DPairsListWave,selWave=M_Selected3DPairsSelWave,titleWave=W_Selected3DPairsTitleWave,mode=4
		Button BTAdd3DCalibrationPos,pos={275,85},size={50,20},title=">>",proc=BTAdd3DAstigmPosProc
		Button BTRemove3DCalibrationPos,pos={275,123},size={50,20},title="<<",proc=BTRemove3DAstigmPosProc
		PopupMenu PM3DCalibrationDataFormat,pos={8,218},size={331,20},bodyWidth=225,title="Calibration data format:"
		PopupMenu PM3DCalibrationDataFormat,mode=1,value= #"\"Beads in a single z-plane;Matrix of beads at different depths\"",disable=1
		SetVariable SVCameraPixelSize,pos={464,219},size={200,15},title="Camera pixel size (nm, required):",value=_NUM:100
		Button BTRegularZSpacing,pos={363,241},size={175,20},title="Regular z spacing",proc=BTRegularZCalibrationDepthProc
		Button BTMake3DCalibration,pos={586,241},size={75,20},title="Do it",proc=BTDo3DCalibrationProc
	endif
End

Function BTAdd3DAstigmPosProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			wave/T W_EligibleAstigmPosListWave = packageFolder:W_EligibleAstigmPosListWave
			wave W_EligibleAstigmPosSelWave = packageFolder:W_EligibleAstigmPosSelWave
			wave/T M_Selected3DPairsListWave = packageFolder:M_Selected3DPairsListWave
			wave M_Selected3DPairsSelWave = packageFolder:M_Selected3DPairsSelWave
			
			// find out which positions are selected for adding
			string selectedPositionsNames = ""
			
			variable i
			for (i = 0; i < DimSize(W_EligibleAstigmPosListWave, 0); i+=1)
				if (W_EligibleAstigmPosSelWave[i] & 2^0)
					selectedPositionsNames += W_EligibleAstigmPosListWave[i] + ";"
				endif
			endfor
			
			variable nSelectedPositionsWaves = ItemsInList(selectedPositionsNames)
			if (nSelectedPositionsWaves == 0)
				return 0
			endif
			
			// only add the ones that are not already scheduled for inclusion
			string thisPosName, newSelectedPosNames = ""
			for (i = 0; i < nSelectedPositionsWaves; i += 1)
				thisPosName = StringFromList(i, selectedPositionsNames)
				FindValue /TEXT=thisPosName M_Selected3DPairsListWave
				if (V_Value == -1)
					newSelectedPosNames += thisPosName + ";"
				endif
			endfor
			
			variable nNewPosNames = ItemsInList(newSelectedPosNames)
			if (nNewPosNames == 0)
				return 0
			endif
			
			variable nCurrentCalibrationPos = DimSize(M_Selected3DPairsListWave, 0)
			Redimension /N=(nCurrentCalibrationPos + nNewPosNames, -1) M_Selected3DPairsListWave, M_Selected3DPairsSelWave
			M_Selected3DPairsListWave[nCurrentCalibrationPos,][0] = StringFromList(p - nCurrentCalibrationPos, newSelectedPosNames)
			M_Selected3DPairsListWave[nCurrentCalibrationPos,][1] = ""
			M_Selected3DPairsSelWave[nCurrentCalibrationPos,][0] = 0
			M_Selected3DPairsSelWave[nCurrentCalibrationPos,][1] = 2^1
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTRemove3DAstigmPosProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			wave/T M_Selected3DPairsListWave = packageFolder:M_Selected3DPairsListWave
			wave M_Selected3DPairsSelWave = packageFolder:M_Selected3DPairsSelWave
			
			// loop over all positions, starting from the last, and delete the ones that are selected
			variable i
			for (i = DimSize(M_Selected3DPairsListWave, 0) - 1; i >= 0; i-=1)
				if (M_Selected3DPairsSelWave[i] & 2^0)
					DeletePoints /M=1 i, 1, M_Selected3DPairsSelWave, M_Selected3DPairsListWave
				endif
			endfor
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTRegularZCalibrationDepthProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			string windowName = ba.win
			
			variable firstDepth, depthDelta
			Prompt firstDepth, "Depth of first plane (nm):"
			Prompt depthDelta, "Distance between planes (nm):"
			DoPrompt "Enter spacing params", firstDepth, depthDelta
			if (V_flag != 0)
				return 0
			endif
			
			wave/T M_ListWave = packageFolder:M_Selected3DPairsListWave
			
			variable i, thisDepth
			for (i = 0; i < DimSize(M_ListWave, 0); i+=1)
				thisDepth = firstDepth + i * depthDelta
				M_ListWave[i][1] = num2str(thisDepth)
			endfor
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTDo3DCalibrationProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			DFREF positionsFolder = root:'Localized Positions'
			DFREF packageFolder = root:Packages:Localizer
			
			wave/T M_Selected3DPairsListWave = packageFolder:M_Selected3DPairsListWave
			
			// no calibration positions? Do nothing
			if (DimSize(M_Selected3DPairsListWave, 0) == 0)
				return 0
			endif
			
			// pixel size
			ControlInfo /W=$ba.win SVCameraPixelSize
			variable pixelSize = V_value
			
			// single-plane or matrix calibration?
			variable isSinglePlaneCalibration
			ControlInfo /W=$ba.win PM3DCalibrationDataFormat
			StrSwitch (S_Value)
				case "Beads in a single z-plane":
					isSinglePlaneCalibration = 1
					break
				case "Matrix of beads at different depths":
					isSinglePlaneCalibration = 0
					break
				default:
					Abort "Unknown calibration type string"
					break
			EndSwitch
			
			// check for errors and format data appropriately
			variable nCalibrationPairs = DimSize(M_Selected3DPairsListWave, 0)
			Make /FREE/N=(nCalibrationPairs) /WAVE W_Positions
			Make /FREE /N=(nCalibrationPairs) /D W_FocusDepths
			
			variable i
			for (i = 0; i < nCalibrationPairs; i+=1)
				wave thisPos = GetPositionsWaveReference(M_Selected3DPairsListWave[i][0])
				// if the positions have no pixel size info associated with them, add it based on the user's input
				if (NumType(NumberByKey("X PIXEL SIZE", note(thisPos)) == 2))
					Note thisPos, "X PIXEL SIZE:" + num2str(pixelSize) + ";"
				endif
				variable thisDepth = str2num(M_Selected3DPairsListWave[i][1])
				if (NumType(thisDepth) != 0)
					Abort "Invalid z position entered for one or more pairs"
				endif
				
				W_Positions[i] = thisPos
				W_FocusDepths[i] = thisDepth
			endfor
			
			// do the actual work (also sets up for visualization of the results)
			variable err = Do3DCalibration(W_Positions, W_FocusDepths, pixelSize, isSinglePlaneCalibration)
			if (err == 0)
				DoWindow /K Make3DCalibrationPanel
			endif
			break
	endswitch

	return 0
End

Function Make3DCalibration_SinglePos()

	// only positions fitted using astigmatism are eligible
	string eligiblePositionsNames = GetGaussAstigPositionsWaves()
	string positionsName
	Prompt positionsName, "Calibration positions:", popup, eligiblePositionsNames
	DoPrompt "Select positions to use", positionsName
	if (V_flag == 1)
		return 0
	endif
	wave posWave = GetPositionsWaveReference(positionsName)
	
	NewDataFolder /O root:Packages
	NewDataFolder /O root:Packages:Localizer
	DFREF packageFolder = root:Packages:Localizer
	
	string /G packageFolder:Cal3DSinglePosName
	SVAR Cal3DSinglePosName = packageFolder:Cal3DSinglePosName
	Cal3DSinglePosName = positionsName
	
	variable nPositions = DimSize(posWave, 0)
	variable firstFrame = posWave[0][0]
	variable lastFrame = posWave[nPositions - 1][0]
	
	// make list and selection waves for the listbox
	Make /N=(0, 2) /O /T packageFolder:M_SinglePos3DListWave
	Make /N=(0, 2) /O /D packageFolder:M_SinglePos3DSelWave
	Make /N=2 /O /T packageFolder:W_SinglePos3DTitleWave = {"Frame", "z position (nm)"}
	wave/T M_SinglePos3DListWave = packageFolder:M_SinglePos3DListWave
	wave M_SinglePos3DSelWave = packageFolder:M_SinglePos3DSelWave
	wave/T W_SinglePos3DTitleWave = packageFolder:W_SinglePos3DTitleWave
	
	// extract all frames for which we have positions available
	MatrixOP /FREE W_FrameIndices = col(posWave, 0)
	variable i, framePos, nFramesFound = 0
	for (i = firstFrame; i <= lastFrame; i+=1)
		framePos = BinarySearch(W_FrameIndices, i)
		if (W_FrameIndices[framePos] == i)	// have positions for this frame?
			Redimension /N=(nFramesFound + 1, -1) M_SinglePos3DListWave, M_SinglePos3DSelWave
			M_SinglePos3DListWave[nFramesFound][0] = num2str(i)
			M_SinglePos3DListWave[nFramesFound][1] = ""
			M_SinglePos3DSelWave[nFramesFound][0] = 0
			M_SinglePos3DSelWave[nFramesFound][1] = 0x02
			nFramesFound += 1
		endif
	endfor
	
	DoWindow /F Make3DCalibration_SinglePanel
	if (V_flag != 1)
		NewPanel /K=1 /W=(515,138,866,394) /N=Make3DCalibration_SinglePanel
		ListBox LBSelected3DCalibrationData,win=$S_name,pos={8,9},size={327,181},mode=0
		ListBox LBSelected3DCalibrationData,win=$S_name,titleWave=W_SinglePos3DTitleWave, listWave=M_SinglePos3DListWave, selWave=M_SinglePos3DSelWave
		SetVariable SVCameraPixelSize,win=$S_name,pos={137,199},size={200,15},title="Camera pixel size (nm, required):"
		SetVariable SVCameraPixelSize,win=$S_name,value= _NUM:100
		Button BTRegularZSpacing,win=$S_name,pos={9,221},size={175,20},proc=BTRegularZDepth_SinglePosProc,title="Regular z spacing"
		Button BTMake3DCalibration,win=$S_name,pos={259,221},size={75,20},proc=BTDo3DCalibrationProc_SinglePos,title="Do it"
	endif

End

Function BTRegularZDepth_SinglePosProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			DFREF packageFolder = root:Packages:Localizer
			
			variable firstDepth, depthDelta
			Prompt firstDepth, "Depth of first plane (nm):"
			Prompt depthDelta, "Distance between planes (nm):"
			DoPrompt "Enter spacing params", firstDepth, depthDelta
			if (V_flag != 0)
				return 0
			endif
			
			wave/T M_ListWave = packageFolder:M_SinglePos3DListWave
			
			variable i, thisDepth, thisFrameIndex
			for (i = 0; i < DimSize(M_ListWave, 0); i+=1)
				thisFrameIndex = str2num(M_ListWave[i][0])
				thisDepth = firstDepth +thisFrameIndex * depthDelta
				M_ListWave[i][1] = num2str(thisDepth)
			endfor
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTDo3DCalibrationProc_SinglePos(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			DFREF positionsFolder = root:'Localized Positions'
			DFREF packageFolder = root:Packages:Localizer
			
			wave/T M_SinglePos3DListWave = packageFolder:M_SinglePos3DListWave
			SVAR Cal3DSinglePosName = packageFolder:Cal3DSinglePosName
			
			wave positions = GetPositionsWaveReference(Cal3DSinglePosName)
			
			// pixel size
			ControlInfo /W=$ba.win SVCameraPixelSize
			variable pixelSize = V_value
			
			// check for errors and format data appropriately
			variable nCalibrationPairs = DimSize(M_SinglePos3DListWave, 0)
			Make /FREE/N=(nCalibrationPairs) /WAVE W_Positions
			Make /FREE /N=(nCalibrationPairs) /D W_FocusDepths
			
			variable i
			for (i = 0; i < nCalibrationPairs; i+=1)
				variable frameNumber = str2num(M_SinglePos3DListWave[i][0])
				variable thisDepth = str2num(M_SinglePos3DListWave[i][1])
				if (NumType(thisDepth) != 0)
					Abort "Invalid z position entered for one or more pairs"
				endif
				
				wave /wave W_ExtractedPos = ExtractPositionsInFrame(positions, frameNumber)
				if (NumType(thisDepth) != 0)
					Abort "Invalid z position entered for one or more pairs"
				endif
				
				W_Positions[i] = W_ExtractedPos[0]
				W_FocusDepths[i] = thisDepth
			endfor
			
			// if the positions have no pixel size info associated with them, add it based on the user's input
			if (NumType(NumberByKey("X PIXEL SIZE", note(positions)) == 2))
				Note positions, "X PIXEL SIZE:" + num2str(pixelSize) + ";"
			endif
			
			// do the actual work (also sets up for visualization of the results)
			variable err = Do3DCalibration(W_Positions, W_FocusDepths, pixelSize, 1)
			if (err == 0)
				DoWindow /K Make3DCalibration_SinglePanel
			endif
			break
	endswitch

	return 0
End

Function Do3DCalibration(W_Positions, W_FocusDepths, pixelSize, isSinglePlane)
	wave /WAVE W_Positions
	wave W_FocusDepths
	variable pixelSize, isSinglePlane
	
	DFREF packageFolder = root:Packages:Localizer
	// get an output wave name
	NewDataFolder /O root:'Calibration 3D'
	DFREF outputDF = root:'Calibration 3D'
	string outputName = GetNew3DCalibrationWaveName("", outputDF, "Select an output name")
	if (strlen(outputName) == 0)
		return kUserAbort	// user cancel
	endif
	
	if (isSinglePlane)
		wave W_Calibration = Do3DCalibrationWork_SinglePlane(W_Positions, W_FocusDepths)
	else
		wave W_Calibration = Do3DCalibrationWork_MatrixBeads(W_Positions, W_FocusDepths)
	endif
	Note /K W_Calibration, "TYPE:Astigmatic3DCalibrationMap;X PIXEL SIZE:" + num2str(pixelSize) + ";"
	Duplicate /O W_Calibration, outputDF:$outputName
	
	// create a visualization of the calibration data	
	wave W_StdDev1 = packageFolder:W_StdDev1
	wave W_StdDev2 = packageFolder:W_StdDev2
	wave W_FocusPositions = packageFolder:W_FocusPositions
	wave W_RotationAngle = packageFolder:W_RotationAngle
	// and of the fits
	Duplicate /O W_StdDev1, packageFolder:W_StdDev1_Fit
	Duplicate /O W_StdDev2, packageFolder:W_StdDev2_Fit
	wave W_StdDev1_Fit = packageFolder:W_StdDev1_Fit
	wave W_StdDev2_Fit = packageFolder:W_StdDev2_Fit
	Make /N=(5) /D/FREE W_FitParams
	W_FitParams = W_Calibration[p + 3]
	W_StdDev1_Fit = AstigDefocusFitFunction(W_FitParams, W_FocusPositions[p])
	W_FitParams = W_Calibration[p + 8]
	W_StdDev2_Fit = AstigDefocusFitFunction(W_FitParams, W_FocusPositions[p])
	
	DoWindow /F CalibrationAxisRatios
	if (V_flag != 1)
		Display /N=CalibrationAxisRatios /K=1 W_StdDev1 vs W_FocusPositions
		AppendToGraph /W=CalibrationAxisRatios W_StdDev2 vs W_FocusPositions
		Label Bottom, "Focus position (nm)"
		Label Left, "Standard deviation along primary axes"
		ModifyGraph mode=3
		AppendToGraph /W=CalibrationAxisRatios W_StdDev1_Fit vs W_FocusPositions
		AppendToGraph /W=CalibrationAxisRatios W_StdDev2_Fit vs W_FocusPositions
	endif
	DoWindow /F CalibrationRotationAngles
	if (V_flag != 1)
		Display /N=CalibrationRotationAngles /K=1 W_RotationAngle vs W_FocusPositions
		Label Bottom, "Focus position (nm)"
		Label Left, "Angle of rotation"
		ModifyGraph mode=3
		AutoPositionWindow /R=CalibrationAxisRatios /M=0 CalibrationRotationAngles
	endif
	
	return 0
End

Function Apply3DCalibrationMap_Menu()
	NewDataFolder /O root:'Calibration 3D'
	DFREF calibrationFolder = root:'Calibration 3D'
	DFREF positionsFolder = root:'Localized Positions'
	
	// start by getting a list of positions waves and calibration maps
	string positionsWaves = GetPossiblePositionsWaves()
	string calibrationMaps = GetPossibleAstig3DCalibrMaps()
	
	// also check for the "* no positions *" return
	if (StringMatch(positionsWaves, "*No Positions*"))
		positionsWaves = ""
	endif
	
	variable nPositionsWaves = ItemsInList(positionsWaves)
	variable nCalibrationMaps = ItemsInList(calibrationMaps)
	
	// only positions fitted using astigmatism are eligible
	string eligiblePositionsNames = ""
	variable i
	for (i = 0; i < nPositionsWaves; i+=1)
		wave thesePositions = GetPositionsWaveReference(StringFromList(i, positionsWaves))
		if (NumberByKey("LOCALIZATION METHOD", note(thesePositions)) == LOCALIZATION_ELLIPSGAUSS_ASTIG)
			eligiblePositionsNames += StringFromList(i, positionsWaves) + ";"
		endif
	endfor
	
	nPositionsWaves = ItemsInList(eligiblePositionsNames)
	if (nPositionsWaves == 0)
		Abort "No positions fitted using astigmatism seem to be present"
	endif
	if (nCalibrationMaps == 0)
		Abort "No 3D calibration maps seem to be present"
	endif
	
	string posName, calibrationMapName
	Prompt posName, "Positions to process:", popup, eligiblePositionsNames
	Prompt calibrationMapName, "Calibration map to use", popup, calibrationMaps
	DoPrompt "Select calculation parameters", posName, calibrationMapName
	if (V_flag == 1)
		return 0
	endif
	
	wave posWave = GetPositionsWaveReference(posName)
	wave calibrationWave = calibrationFolder:$calibrationMapName
	
	string newOutputName = GetNewPositionsWaveName(posName + "_3D", positionsFolder, "Enter the name of the output wave")
	if (strlen(newOutputName) == 0)
		return 0
	endif
	
	wave CalibratedPositions = Apply3DCalibrationMap(posWave, calibrationWave)
	Duplicate /O CalibratedPositions, positionsFolder:$newOutputName
End

Function BTCorrectDriftProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string baseWindowName = GetBaseWindowName(ba.win)
			DoDriftCorrection(calledFromButton=1, baseWindowName=baseWindowName)
			break
	EndSwitch
End

Function DoDriftCorrection([calledFromButton, baseWindowName])
	variable calledFromButton
	string baseWindowName
	
	if (ParamIsDefault(calledFromButton))
		calledFromButton = 0
	endif
	
	string correctionModeName
	string positionsWaveName
	if (calledFromButton)
		ControlInfo /W=$baseWindowName#OtherControls PMDriftCorrectionMode
		correctionModeName = S_Value
		
		ControlInfo /W=$baseWindowName#OtherControls PMDriftCorrSelectPositionsWave
		positionsWaveName = S_Value	
		if (StringMatch(positionsWaveName, "* No Positions *") == 1)
			Abort "Please select positions to correct first"
		endif
		
		wave positions = GetPositionsWaveReference(positionsWaveName)
		if (WaveExists(positions) != 1)
			Abort "The selected positions do not appear to be valid"
		endif
	else
		Prompt positionsWaveName, "Positions to correct:", popup, GetPossiblePositionsWaves()
		Prompt correctionModeName, "Algorithm:", popup, "subimages;ad-hoc fiducials;average position;"	// no manual fiducials when called from a menu
		DoPrompt "Select parameters", positionsWaveName, correctionModeName
		if (V_flag != 0)
			return 0
		endif
		
		wave positions = GetPositionsWaveReference(positionsWaveName)
		if (WaveExists(positions) != 1)
			Abort "The selected positions do not appear to be valid"
		endif
	endif
	
	// check if the positions have been consolidated, if so then refuse
	variable nFramesPresentCol
	GetColumnForNFramesPresent(positions, nFramesPresentCol)
	ImageTransform /G=(nFramesPresentCol) sumCol, positions
	if (V_value != DimSize(positions, 0))
		Abort "The requested positions have been consolidated, there is no point in performing the drift correction"
	endif
	
	StrSwitch (correctionModeName)
		case "ad-hoc fiducials":
			wave /WAVE /Z W_Result = DoAutomaticDriftCorrection(positions)
			if (WaveExists(W_Result) != 1)
				return 0	// cancel
			endif
			wave M_DriftEstimates = W_Result[0]
			wave W_nMarkersUsed = W_Result[1]
			break
		case "manual fiducials":
			DoManualDriftCorrection(baseWindowName)
			// manual drift correction is special because it is modeless, and needs to call the downstream functionality itself
			return 0
			break
		case "subimages":
			wave /Z M_DriftEstimates = DoSubImageDriftCorrection(positions)
			break
		case "average position":
			wave /Z M_DriftEstimates = DoAveragePositionDriftCorr(positions)
			break
		default:
			Abort "Unknown correction mode in BTCorrectDriftProc()"
			break
	EndSwitch
	
	if (WaveExists(M_DriftEstimates) != 1)	// cancel
		return 0
	endif
	
	if (!WaveExists(W_nMarkersUsed))
		InterpAndApplyDriftCorrection(positions, M_DriftEstimates)
	else
		InterpAndApplyDriftCorrection(positions, M_DriftEstimates, W_nEmittersUsed = W_nMarkersUsed)
	endif
End

Function /WAVE DoManualDriftCorrection(baseWindowName)
	string baseWindowName
	
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	DFREF packageFolder = root:Packages:Localizer
	NVAR /Z driftCorrectionMaxJump = windowDataFolder:V_driftCorrectionMaxJump
	if (!NVAR_Exists(driftCorrectionMaxJump))
		variable /G windowDataFolder:V_driftCorrectionMaxJump
		NVAR driftCorrectionMaxJump = windowDataFolder:V_driftCorrectionMaxJump
		driftCorrectionMaxJump = 1
	endif
	variable msTimer, sec, err
	string positionsWaves, panelName
	
	// set up the waves that will contain the clicked positions and the listbox entries (both updated by the window hook function)
	Make /D/O/N=(0, 2) windowDataFolder:M_ClickedCoordinates
	Make /T/O/N=(0, 2) windowDataFolder:M_ClickedTextCoordinates
	Make /T/O/N=(2) windowDataFolder:M_DriftListBoxTitles = {"X", "Y"}
	
	wave M_ClickedCoordinates = windowDataFolder:M_ClickedCoordinates
	wave /T M_ClickedTextCoordinates = windowDataFolder:M_ClickedTextCoordinates
	wave /T M_DriftListBoxTitles = windowDataFolder:M_DriftListBoxTitles
	
	// draw the main panel
	DoWindow /K DriftCorrection	// in case the user still had a panel floating around
	
	NewPanel /K=1 /N=DriftCorrection /W=(949,50,1197,306)
	panelName = S_name
	SetDrawLayer UserBack
	SetDrawEnv fsize= 11
	DrawText 17,38,"1. Make sure the first 'good' frame\rin the movie is displayed"
	SetDrawEnv fsize= 11
	DrawText 18,65,"2. Click the fiducial markers in the image."
	SetDrawEnv fsize= 11
	DrawText 22,222,"3. Click continue"
	ListBox LBDriftCorrectionMarkerPosition,pos={68,78},size={107,118}, listWave = M_ClickedTextCoordinates, titleWave=M_DriftListBoxTitles
	Button BTDriftCorrectionContinue,pos={74,227},size={100,20},proc=BTDriftMarkerContinue,title="Continue"
	SetWindow $panelName, hook(KillHook) = DriftCorrection_WindowHook
	
	// set the name of the corresponding interface panel in the userdata
	// since it's difficult to retrieve otherwise
	SetWindow $panelName, UserData = baseWindowName
	
	AutoPositionWindow /R=$baseWindowName /M=0 $panelName
	// attach the window hook function to the CCD Viewer
	SetWindow $baseWindowName, hook(DriftMarkerHook) = CaptureClickedCoordinates_Hook
End

Function DriftCorrection_WindowHook(s)
	STRUCT WMWinHookStruct &s
	
	switch (s.eventcode)
		case 2:	// window kill, remove the window hook function from CCDViewer
			string interfacePanelName = GetUserData(s.winName, "", "")
			DoWindow $interfacePanelName
			if (V_flag != 0)
				SetWindow $interfacePanelName hook(DriftMarkerHook) = $""
				RemoveFromGraph /Z /W=$interfacePanelName#CCDViewer M_ClickedCoordinates
			endif
			break
		default:
			return 0
	endswitch
	
	return 0
End

Function CaptureClickedCoordinates_Hook(s) 
	STRUCT WMWinHookStruct &s
	// this function attaches to a window of choice and appends the X- and Y-coordinates where the user clicked the mouse to a specific wave
	// both a numeric and a text version of the wave is made (for use with listbox controls)
	
	switch (s.eventcode)
		case 3:	// mousedown
			variable xVal, yVal
			string baseWindowName = GetBaseWindowName(s.winName)
			DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
			
			wave /Z M_ClickedCoordinates = windowDataFolder:M_ClickedCoordinates
			wave /T/Z M_ClickedTextCoordinates = windowDataFolder:M_ClickedTextCoordinates
			
			// get the axis range of the viewer graph so clicks outside the graph window can be rejected
			variable minX, minY, maxX, maxY
			GetAxis /Q/W=$baseWindowName#CCDViewer Bottom
			minX = V_min
			maxX = V_max
			GetAxis /Q/W=$baseWindowName#CCDViewer Left
			minY = V_min
			maxY = V_max
			
			if (WaveExists(M_ClickedCoordinates) == 0)
				Make /D/N=(0, 2) windowDataFolder:M_ClickedCoordinates
				wave M_ClickedCoordinates = windowDataFolder:M_ClickedCoordinates
			endif
			if (WaveExists(M_ClickedTextCoordinates) == 0)
				Make /T /N=(0, 2) windowDataFolder:M_ClickedTextCoordinates
				wave /T M_ClickedTextCoordinates = windowDataFolder:M_ClickedTextCoordinates
			endif
			
			if (DimSize(M_ClickedCoordinates, 0) == 0)
				// no clicked coordinates yet
				// append the wave to the graph so that the user can see
				// where they clicked
				AppendToGraph /W=$baseWindowName#CCDViewer M_ClickedCoordinates[][1] vs M_ClickedCoordinates[][0]
				ModifyGraph /W=$baseWindowName#CCDViewer mode(M_ClickedCoordinates)=3
			endif
			
			xVal = AxisValFromPixel(baseWindowName + "#CCDViewer", "Bottom", s.mouseLoc.h)
			yVal = AxisValFromPixel(baseWindowName + "#CCDViewer", "Left", s.mouseLoc.v)
			
			if ((xVal < minX) || (xVal > maxX) || (yVal < minY) || (yVal > maxY))
				return 1
			endif
			
			Redimension /N=(DimSize(M_ClickedCoordinates, 0) + 1, 2) M_ClickedCoordinates
			M_ClickedCoordinates[DimSize(M_ClickedCoordinates, 0) - 1][0] = xVal
			M_ClickedCoordinates[DimSize(M_ClickedCoordinates, 0) - 1][1] = yVal
			
			Redimension /N=(DimSize(M_ClickedTextCoordinates, 0) + 1, 2) M_ClickedTextCoordinates
			M_ClickedTextCoordinates[DimSize(M_ClickedTextCoordinates, 0) - 1][0] = num2str(xVal)
			M_ClickedTextCoordinates[DimSize(M_ClickedTextCoordinates, 0) - 1][1] = num2str(yVal)
			
			return 1
			break
		default:
			return 0
			break
	endswitch
End

Function BTDriftMarkerContinue(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			string interfacePanelName = GetUserData(ba.win, "", "")
			DFREF windowDataFolder = GetWindowDataFolder(interfacePanelName)
			DFREF packageFolder = root:Packages:Localizer
			NVAR currentFrame = windowDataFolder:V_currentImage
			NVAR maxJump = windowDataFolder:V_driftCorrectionMaxJump
			
			wave /T M_ClickedTextCoordinates = windowDataFolder:M_ClickedTextCoordinates
			wave M_ClickedCoordinates = windowDataFolder:M_ClickedCoordinates
			
			variable nMarkers = DimSize(M_ClickedTextCoordinates, 0), nPositions
			string positionsWaveName
			variable i, nearestX, nearestY, index
			
			ControlInfo /W=$interfacePanelName#OtherControls PMDriftCorrSelectPositionsWave
			positionsWaveName = S_Value
			
			DoWindow /K DriftCorrection
			
			if (StringMatch(positionsWaveName, "* No Positions *") == 1)
				KillWaves /Z M_ClickedCoordinatesText
				KillWaves /Z windowDataFolder:M_ClickedCoordinates
				return 0
			endif
			
			wave positions = GetPositionsWaveReference(positionsWaveName)
			nPositions = DimSize(positions, 0)
			
			if (nMarkers == 0)	// no markers indicated
				KillWaves /Z M_ClickedCoordinatesText
				KillWaves /Z M_ClickedCoordinates
				return 0
			endif
			
			// check if the positions have been consolidated, if so then refuse
			variable nFramesPresentCol
			GetColumnForNFramesPresent(positions, nFramesPresentCol)
			ImageTransform /G=(nFramesPresentCol) sumCol, positions
			if (V_value != DimSize(positions, 0))
				Abort "The requested positions have been consolidated, there is no point in performing the drift correction"
			endif
			
			// extract the positions that are fitted in the current frame
			wave /WAVE M_ExtractedPositions = ExtractPositionsInFrame(positions, currentFrame)
			wave extractedPositions = M_ExtractedPositions[0]
			wave extractedPositionsIndices = M_ExtractedPositions[1]
			
			if (DimSize(extractedPositions, 0) == 0)	// no positions were fitted in this frame
				KillWaves /Z M_ClickedCoordinatesText, extractedPositions
				KillWaves /Z M_ClickedCoordinates
				Abort "No positions were fitted in the current frame, not even the markers!"
			endif
			
			Make /D/O/N=(nMarkers, 2) windowDataFolder:RecoveredMarkerPositions	// the coordinates that we have recovered for the markers
			wave RecoveredMarkerPositions = windowDataFolder:RecoveredMarkerPositions
			RecoveredMarkerPositions = str2num(M_ClickedTextCoordinates[p][q])
			
			for (i = 0; i < nMarkers; i+=1)
				returnPointNearestFitPositions(extractedPositions, RecoveredMarkerPositions[i][0], RecoveredMarkerPositions[i][1], nearestX, nearestY, index)
				RecoveredMarkerPositions[i][0] = nearestX
				RecoveredMarkerPositions[i][1] = nearestY
			endfor
			
			// the user might have clicked the same marker more than once
			// so delete identical positions
			variable offset = 0, j, k
			for (i = 0; i < DimSize(RecoveredMarkerPositions, 0); i+=1)
				for (j = i + 1; j < DimSize(RecoveredMarkerPositions, 0); j+=1)
					if ((RecoveredMarkerPositions[i][0] == RecoveredMarkerPositions[j][0]) && (RecoveredMarkerPositions[i][1] == RecoveredMarkerPositions[j][1]))
						// remove the markers pointed to by j
						Duplicate /FREE RecoveredMarkerPositions, OldRecoveredMarkerPositions
						Redimension /N=(DimSize(OldRecoveredMarkerPositions, 0) - 1, -1) RecoveredMarkerPositions
						for (k = 0; k < DimSize(OldRecoveredMarkerPositions, 0); k+=1)
							if (offset == j)
								continue
							endif
							RecoveredMarkerPositions[offset][] = OldRecoveredMarkerPositions[k][q]
							offset += 1
						endfor
						j -= 1
					endif
				endfor
			endfor
					
			nMarkers = DimSize(RecoveredMarkerPositions, 0)
			
			Make /D/O/N=(nMarkers) windowDataFolder:W_Markers_X
			Make /D/O/N=(nMarkers) windowDataFolder:W_Markers_Y
			wave W_Markers_X = windowDataFolder:W_Markers_X
			wave W_Markers_Y = windowDataFolder:W_Markers_Y
			
			W_Markers_X = RecoveredMarkerPositions[p][0]
			W_Markers_Y = RecoveredMarkerPositions[p][1]
			
			DoWindow /F $interfacePanelName
			AppendToGraph /W=$interfacePanelName#CCDViewer W_Markers_Y vs W_Markers_X
			ModifyGraph /W=$interfacePanelName#CCDViewer mode(W_Markers_Y)=3
			ModifyGraph /W=$interfacePanelName#CCDViewer rgb(W_Markers_Y)=(0,65535,0)
			
			DoUpdate
			
			DoAlert 1, "The green markers on the image highlight the presumed marker positions. Accept?"
			if (V_flag != 1)	// user didn't choose 'yes'
				RemoveFromGraph /W=$interfacePanelName#CCDViewer W_Markers_Y
				KillWaves /Z M_ClickedCoordinatesText, extractedPositions, RecoveredMarkerPositions, W_Markers_X, W_Markers_Y
				KillWaves /Z windowDataFolder:M_ClickedCoordinates
				return 0
			endif
			
			RemoveFromGraph /W=$interfacePanelName#CCDViewer W_Markers_Y
			KillWaves /Z W_Markers_X, W_Markers_Y
			
			// do the real work
			wave /WAVE W_Result = CalculateDrift_ExplicitMarkers(positions, RecoveredMarkerPositions, maxJump)
			wave M_DriftParameters = W_Result[0]
			wave W_nMarkersUsed = W_Result[1]
			InterpAndApplyDriftCorrection(positions, M_DriftParameters, W_nEmittersUsed = W_nMarkersUsed)
	endswitch
	return 0
End

Function /WAVE DoAutomaticDriftCorrection(positions)
	wave positions
	
	DFREF packageFolder = root:Packages:Localizer
	NVAR /Z maxJump = packageFolder:V_driftCorrectionMaxJump
	if (!NVAR_Exists(maxJump))
		variable /G packageFolder:V_driftCorrectionMaxJump = 1
		NVAR maxJump = packageFolder:V_driftCorrectionMaxJump
	endif
	NVAR /Z autoDriftCorrectionMinLength = packageFolder:V_autoDriftCorrectionMinLength
	if (!NVAR_Exists(autoDriftCorrectionMinLength))
		variable /G packageFolder:V_autoDriftCorrectionMinLength = 5
		NVAR autoDriftCorrectionMinLength = packageFolder:V_autoDriftCorrectionMinLength
	endif
	NVAR /Z autoDriftCorrectionNBlinking = packageFolder:V_autoDriftCorrectionNBlinking
	if (!NVAR_Exists(autoDriftCorrectionNBlinking))
		variable /G packageFolder:V_autoDriftCorrectionNBlinking = 0
		NVAR autoDriftCorrectionNBlinking = packageFolder:V_autoDriftCorrectionNBlinking
	endif
	
	variable maxJumpLocal = maxJump
	variable trackLengthLocal = autoDriftCorrectionMinLength
	variable nBlinkingLocal = autoDriftCorrectionNBlinking
	
	Prompt maxJumpLocal, "Max positions jump:"
	Prompt trackLengthLocal, "Min number of frames emitter must be visible:"
	Prompt nBlinkingLocal, "Emitter may disappear for ? frames:"
	DoPrompt "Enter calculation parameters", maxJumpLocal, trackLengthLocal, nBlinkingLocal
	if (V_flag == 1)
		return $""
	endif
	
	maxJump = maxJumpLocal
	autoDriftCorrectionMinLength = trackLengthLocal
	autoDriftCorrectionNBlinking = nBlinkingLocal
	
	wave /WAVE W_Result = DriftCorrection_automatic(positions, maxJumpLocal, trackLengthLocal, nBlinkingLocal)
	return W_Result
End

Function /WAVE DoSubImageDriftCorrection(positions)
	wave positions
	
	DFREF packageFolder = root:Packages:Localizer
	NVAR /Z minPosPerSubImage = packageFolder:V_driftCorrMinPosPerSubImage
	if (!NVAR_Exists(minPosPerSubImage))
		variable /G packageFolder:V_driftCorrMinPosPerSubImage = 5000
		NVAR minPosPerSubImage = packageFolder:V_driftCorrMinPosPerSubImage
	endif
	NVAR /Z subImagePixelSize = packageFolder:V_driftCorrSubImagePixelSize
	if (!NVAR_Exists(subImagePixelSize))
		variable /G packageFolder:V_driftCorrSubImagePixelSize = 0.25
		NVAR subImagePixelSize = packageFolder:V_driftCorrSubImagePixelSize
	endif
	NVAR /Z subImageUseThreads = packageFolder:V_driftCorrSubImageUseThreads
	if (!NVAR_Exists(subImageUseThreads))
		variable /G packageFolder:V_driftCorrSubImageUseThreads = 0
		NVAR /Z subImageUseThreads = packageFolder:V_driftCorrSubImageUseThreads
	endif
	
	variable minPositionsPerImageLocal = minPosPerSubImage
	variable subImagePixelSizeLocal = subImagePixelSize
	variable subImageUseThreadsLocal = subImageUseThreads
	
	Prompt minPositionsPerImageLocal, "Min positions per subimage:"
	Prompt subImagePixelSizeLocal, "Subimage pixel size (pixels):"
	Prompt subImageUseThreadsLocal, "Use multiple threads (faster but needs more memory)"
	DoPrompt "Enter calculation parameters", minPositionsPerImageLocal, subImagePixelSizeLocal, subImageUseThreadsLocal
	if (V_flag == 1)
		return $""
	endif
	
	if ((minPositionsPerImageLocal <= 0) || (subImagePixelSizeLocal <= 0))
		Abort "Invalid parameters"
	endif
	
	minPosPerSubImage = minPositionsPerImageLocal
	subImagePixelSize = subImagePixelSizeLocal
	subImageUseThreads = subImageUseThreadsLocal
	
	wave M_DriftParameters = DriftCorrection_SubImage(positions, minPositionsPerImageLocal, subImagePixelSizeLocal, subImageUseThreadsLocal)
	return M_DriftParameters
End


Function /WAVE DoAveragePositionDriftCorr(positions)
	wave positions
	
	DFREF packageFolder = root:Packages:Localizer
	
	NVAR /Z driftCorrNAveragePos = packageFolder:V_driftCorrNAveragePos
	if (!NVAR_Exists(driftCorrNAveragePos))
		variable /G packageFolder:V_driftCorrNAveragePos = 2000
		NVAR driftCorrNAveragePos = packageFolder:V_driftCorrNAveragePos
	endif
	
	variable nPositionsToGroup = driftCorrNAveragePos
	
	Prompt nPositionsToGroup, "Number of positions to group:"
	DoPrompt "Enter calculation parameters", nPositionsToGroup
	if (V_flag == 1)
		return $""
	endif
	
	if ((nPositionsToGroup <= 0) || (nPositionsToGroup >= DimSize(positions, 0)))
		Abort "Invalid parameters"
	endif
	
	driftCorrNAveragePos = nPositionsToGroup
	
	wave M_DriftParameters = DriftCorrection_AvgImage(positions, nPositionsToGroup)
	return M_DriftParameters
End

Function InterpAndApplyDriftCorrection(positions, M_DriftEstimates, [W_nEmittersUsed])
	wave positions
	wave M_DriftEstimates, W_nEmittersUsed
	
	DFREF packageFolder = root:Packages:Localizer:
	NewDataFolder /O packageFolder:DriftInterpolation
	DFREF driftInterpFolder = packageFolder:DriftInterpolation
	
	NVAR /Z testVar = driftInterpFolder:V_DriftInterpPolyFitXOrder
	if (!NVAR_Exists(testVar))
		variable /G driftInterpFolder:V_DriftInterpPolyFitXOrder = 3
		variable /G driftInterpFolder:V_DriftInterpPolyFitYOrder = 3
		variable /G driftInterpFolder:V_DriftInterpSplineSmoothX = 1
		variable /G driftInterpFolder:V_DriftInterpSplineSmoothY = 1
		variable /G driftInterpFolder:V_DriftInterpSplineStdDevX = 0.1
		variable /G driftInterpFolder:V_DriftInterpSplineStdDevY = 0.1
		string /G driftInterpFolder:S_outputPositionsName
	endif
	
	NVAR gDriftInterpPolyFitXOrder = driftInterpFolder:V_DriftInterpPolyFitXOrder
	NVAR gDriftInterpPolyFitYOrder = driftInterpFolder:V_DriftInterpPolyFitYOrder
	NVAR gDriftInterpSplineSmoothX = driftInterpFolder:V_DriftInterpSplineSmoothX
	NVAR gDriftInterpSplineSmoothY = driftInterpFolder:V_DriftInterpSplineSmoothY
	NVAR gDriftInterpSplineStdDevX = driftInterpFolder:V_DriftInterpSplineStdDevX
	NVAR gDriftInterpSplineStdDevY = driftInterpFolder:V_DriftInterpSplineStdDevY
	SVAR gOutputPositionsName = driftInterpFolder:S_outputPositionsName
	
	variable nPositions = DimSize(positions, 0)
	variable firstFrame = positions[0][0]
	variable lastFrame = positions[nPositions - 1][0]
	variable nFrames = lastFrame - firstFrame + 1
	
	// make local copies of the drift estimates for use in the graph
	Make /O/N=(DimSize(M_DriftEstimates, 0)) /D packageFolder:W_EstimatedDriftLocations, packageFolder:W_EstimatedDriftX, packageFolder:W_EstimatedDriftY
	wave W_EstimatedDriftLocations = packageFolder:W_EstimatedDriftLocations
	wave W_EstimatedDriftX = packageFolder:W_EstimatedDriftX
	wave W_EstimatedDriftY = packageFolder:W_EstimatedDriftY
	W_EstimatedDriftLocations = M_DriftEstimates[p][0]
	W_EstimatedDriftX = M_DriftEstimates[p][1]
	W_EstimatedDriftY = M_DriftEstimates[p][2]
	if (!ParamIsdefault(W_nEmittersUsed))
		Duplicate /O W_nEmittersUsed, packageFolder:W_nEmittersUsed
		wave W_nEmittersUsedLocal = packageFolder:W_nEmittersUsed
	endif
	
	DoWindow /K DriftInterpolatePanel
	
	NewPanel /N=DriftInterpolatePanel /K=1 /W=(327,90,1112,574)
	SetDrawLayer UserBack
	DrawText 509,23,"Interpolation method to use:"
	DrawText 21,21,"Raw and interpolated drift estimates:"
	DefineGuide GraphLeft={FL,19},TopGraphTop={FT,29},BottomGraphBottom={FB,-23},GraphRight={FR,-303}
	DefineGuide TopGraphBottom={FT,0.48,FB}, BottomGraphTop={FT,0.52,FB}
	DefineGuide PanelLeft={FR,-281},PanelRight={FR,-22}
	
	Display/N=DriftInterpolateGraph /FG=(GraphLeft,TopGraphTop,GraphRight,TopGraphBottom)/HOST=DriftInterpolatePanel
	Display/N=AppliedCorrectionPreviewGraph /FG=(GraphLeft,BottomGraphTop,GraphRight,BottomGraphBottom)/HOST=DriftInterpolatePanel
	
	NewPanel/N=DriftInterpolateControlPanel /FG=(PanelLeft,TopGraphTop,PanelRight,BottomGraphBottom)/HOST=DriftInterpolatePanel
	TabControl TBInterpolationMode,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={9,18},size={241,24},tabLabel(0)="Raw"
	TabControl TBInterpolationMode,win=DriftInterpolatePanel#DriftInterpolateControlPanel,tabLabel(1)="Poly Fit",tabLabel(2)="Spline"
	TabControl TBInterpolationMode,win=DriftInterpolatePanel#DriftInterpolateControlPanel,value= 0,proc=TBDriftInterpolatePanelProc
	TitleBox TBNoSettingsForLinear,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={93,64},size={69,20},title="<Nothing to modify>",disable=0
	SetVariable SVPolyFitOrderX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={12,57},size={150,15},title="Order for X:",value=gDriftInterpPolyFitXOrder
	SetVariable SVPolyFitOrderX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={1,20,1},proc=SVRecalculateDriftInterpProc,disable=1
	SetVariable SVPolyFitOrderY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={13,86},size={150,15},title="Order for Y:",value=gDriftInterpPolyFitYOrder
	SetVariable SVPolyFitOrderY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={1,20,1},proc=SVRecalculateDriftInterpProc,disable=1
	SetVariable SVStandardDeviationX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={12,140},size={200,15},title="Standard deviation for X:",value=gDriftInterpSplineStdDevX,disable=1,proc=SVRecalculateDriftInterpProc
	SetVariable SVStandardDeviationX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={0,inf,0.1}
	SetVariable SVStandardDeviationY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={14,198},size={200,15},title="Standard deviation for Y:",value=gDriftInterpSplineStdDevY,disable=1,proc=SVRecalculateDriftInterpProc
	SetVariable SVStandardDeviationY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={0,inf,0.1}
	SetVariable SVSplineSmoothingFactorX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={14,116},size={200,15},title="Smoothing factor for X (~1.0):",value=gDriftInterpSplineSmoothX,disable=1,proc=SVRecalculateDriftInterpProc
	SetVariable SVSplineSmoothingFactorX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={0,inf, 0.05}
	SetVariable SVSplineSmoothingFactorY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={12,177},size={200,15},title="Smoothing factor for Y (~1.0):",value=gDriftInterpSplineSmoothY,disable=1,proc=SVRecalculateDriftInterpProc
	SetVariable SVSplineSmoothingFactorY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={0,inf, 0.05}
	Button BTCancelDriftInterpolation,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={26,300},size={75,20},title="Cancel",proc=BTCancelDriftInterpolationProc
	Button BTContinueDriftInterpolation,win=DriftInterpolatePanel#DriftInterpolateControlPanel,pos={144,300},size={75,20},title="Continue",proc=BTContinueDriftInterpProc
	
	// create the waves needed to show the output
	Make /O/N=(nFrames) /D driftInterpFolder:W_InterpolatedDriftX, driftInterpFolder:W_InterpolatedDriftY
	wave W_InterpolatedDriftX = driftInterpFolder:W_InterpolatedDriftX
	wave W_InterpolatedDriftY = driftInterpFolder:W_InterpolatedDriftY
	SetScale /P x, firstFrame, 1, W_InterpolatedDriftX, W_InterpolatedDriftY
	wave CorrectedPositions = ApplyDriftCorrection(positions, W_InterpolatedDriftX, W_InterpolatedDriftY)
	Duplicate /O CorrectedPositions, driftInterpFolder:M_TentativelyCorrectedPositions
	wave M_TentativelyCorrectedPositions = driftInterpFolder:M_TentativelyCorrectedPositions
	Make /O/N=(DimSize(M_TentativelyCorrectedPositions, 0)) /D driftInterpFolder:W_TentativelyCorrectedFrames
	wave W_TentativelyCorrectedFrames = driftInterpFolder:W_TentativelyCorrectedFrames
	W_TentativelyCorrectedFrames = M_TentativelyCorrectedPositions[p][0]
	
	
	// and show them in the graph
	AppendToGraph /W=DriftInterpolatePanel#DriftInterpolateGraph W_EstimatedDriftX vs W_EstimatedDriftLocations	// if the order of these appends changes then the legend needs to be changed as well
	AppendToGraph /W=DriftInterpolatePanel#DriftInterpolateGraph W_EstimatedDriftY vs W_EstimatedDriftLocations
	AppendToGraph /W=DriftInterpolatePanel#DriftInterpolateGraph W_InterpolatedDriftX
	AppendToGraph /W=DriftInterpolatePanel#DriftInterpolateGraph W_InterpolatedDriftY
	
	variable xCol, yCol, zCol
	getColumnsForEmitterPositions(positions, xCol, yCol, zCol)
	AppendToGraph /W=DriftInterpolatePanel#AppliedCorrectionPreviewGraph M_TentativelyCorrectedPositions[][yCol] vs M_TentativelyCorrectedPositions[][xCol]
	
	// fancify the layout
	Label /W=DriftInterpolatePanel#DriftInterpolateGraph bottom, "Acquisitions"
	Label /W=DriftInterpolatePanel#DriftInterpolateGraph left, "Estimated drift (pixels)"
	ModifyGraph /W=DriftInterpolatePanel#DriftInterpolateGraph rgb($NameOfWave(W_EstimatedDriftX))=(0,0,65535),rgb($NameOfWave(W_InterpolatedDriftX))=(0,0,65535)
	ModifyGraph /W=DriftInterpolatePanel#DriftInterpolateGraph mode($NameOfWave(W_EstimatedDriftX))=3,marker($NameOfWave(W_EstimatedDriftX))=8
	ModifyGraph /W=DriftInterpolatePanel#DriftInterpolateGraph mode($NameOfWave(W_EstimatedDriftY))=3,marker($NameOfWave(W_EstimatedDriftY))=8
	Legend /W=DriftInterpolatePanel#DriftInterpolateGraph /C/N=text0/J/F=2 "\\s(#0) Estimated X\r\\s(#1) Estimated Y\r\\s(#2) Interpolated X\r\\s(#3) Interpolated Y"
	ModifyGraph /W=DriftInterpolatePanel#AppliedCorrectionPreviewGraph mode=2, zColor(M_TentativelyCorrectedPositions)={W_TentativelyCorrectedFrames,*,*,YellowHot,0}
	
	
	// optional wave showing the number of markers used in each frame
	if (!ParamIsDefault(W_nEmittersUsed))
		AppendToGraph /W=DriftInterpolatePanel#DriftInterpolateGraph /R W_nEmittersUsedLocal
		AppendText /W=DriftInterpolatePanel#DriftInterpolateGraph /N=text0 "\\s(#4) N markers used per frame"
		ModifyGraph /W=DriftInterpolatePanel#DriftInterpolateGraph rgb($NameOfWave(W_nEmittersUsedLocal))=(3,52428,1)
		Label /W=DriftInterpolatePanel#DriftInterpolateGraph right, "Number of markers used per frame"
	endif
	
	// and format the waves so that RecalculateDriftInterpolation() will find them
	Make /N=(4) /WAVE /O driftInterpFolder:W_DriftInterpolationWaveRefs	// holds wave refs to positions and drift waves
	wave /WAVE W_DriftInterpolationWaveRefs = driftInterpFolder:W_DriftInterpolationWaveRefs
	W_DriftInterpolationWaveRefs[0] = positions
	W_DriftInterpolationWaveRefs[1] = W_EstimatedDriftLocations
	W_DriftInterpolationWaveRefs[2] = W_EstimatedDriftX
	W_DriftInterpolationWaveRefs[3] = W_EstimatedDriftY
	
	// limit the maximum order of polynomial fit
	variable maxFitOrder = min(DimSize(W_EstimatedDriftX, 0) - 1, 10)
	SetVariable SVPolyFitOrderX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={1,maxFitOrder,1}
	SetVariable SVPolyFitOrderY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,limits={1,maxFitOrder,1}
	
	// actualize the view
	RecalculateDriftInterpolation()
End

constant kDriftCorrLineTab=0
constant kDriftCorrPolyFitTab=1
constant kDriftCorrSplineTab=2

Function TBDriftInterpolatePanelProc(tca) : TabControl
	STRUCT WMTabControlAction &tca

	switch( tca.eventCode )
		case 2: // mouse up
			Variable tab = tca.tab
			
			// hide all controls, then show the controls that are needed
			TitleBox TBNoSettingsForLinear,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=1
			SetVariable SVPolyFitOrderX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=1
			SetVariable SVPolyFitOrderY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=1
			SetVariable SVStandardDeviationX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=1
			SetVariable SVStandardDeviationY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=1
			SetVariable SVSplineSmoothingFactorY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=1
			SetVariable SVSplineSmoothingFactorX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=1
			
			switch (tab)
				case kDriftCorrLineTab:
					TitleBox TBNoSettingsForLinear,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=0
					break
				case kDriftCorrPolyFitTab:
					SetVariable SVPolyFitOrderX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=0
					SetVariable SVPolyFitOrderY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=0
					break
				case kDriftCorrSplineTab:
					SetVariable SVStandardDeviationX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=0
					SetVariable SVStandardDeviationY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=0
					SetVariable SVSplineSmoothingFactorY,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=0
					SetVariable SVSplineSmoothingFactorX,win=DriftInterpolatePanel#DriftInterpolateControlPanel,disable=0
					break
				default:
					Abort "Unknown tab in TBDriftInterpolatePanelProc()"
			EndSwitch
			
			RecalculateDriftInterpolation()
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SVRecalculateDriftInterpProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	
	DFREF packageFolder = root:Packages:Localizer
	DFREF driftInterpFolder = packageFolder:DriftInterpolation

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			// poly orders need to be real number
			NVAR gDriftInterpPolyFitXOrder = driftInterpFolder:V_DriftInterpPolyFitXOrder
			NVAR gDriftInterpPolyFitYOrder = driftInterpFolder:V_DriftInterpPolyFitYOrder
			gDriftInterpPolyFitXOrder = round(gDriftInterpPolyFitXOrder)
			gDriftInterpPolyFitYOrder = round(gDriftInterpPolyFitYOrder)
			
			RecalculateDriftInterpolation()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTCancelDriftInterpolationProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			DFREF packageFolder = root:Packages:Localizer
			DFREF driftInterpFolder = packageFolder:DriftInterpolation
			
			DoWindow /K DriftInterpolatePanel
			
			DFREF savDF = GetDataFolderDFR()
			SetDataFolder driftInterpFolder
			KillWaves /A/Z
			SetDataFolder savDF
			KillDataFolder /Z driftInterpFolder
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function BTContinueDriftInterpProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			DFREF packageFolder = root:Packages:Localizer
			DFREF driftInterpFolder = packageFolder:DriftInterpolation
			
			wave /WAVE W_DriftInterpolationWaveRefs = driftInterpFolder:W_DriftInterpolationWaveRefs	// holds wave refs to positions and drift waves
			wave positions = W_DriftInterpolationWaveRefs[0]
			
			wave W_InterpolatedDriftX = driftInterpFolder:W_InterpolatedDriftX
			wave W_InterpolatedDriftY = driftInterpFolder:W_InterpolatedDriftY
			
			// ask for an output name
			DFREF positionsFolder = root:'Localized Positions'
			string outputName = CleanupName(NameOfWave(positions) + "_driftCorr", 0)
			outputName = GetNewPositionsWaveName(outputName, positionsFolder, "Enter a name for the corrected positions wave")
			if (strlen(outputName) == 0)	// the user canceled the dialog
				return 0
			endif
			
			wave CorrectedPositions = ApplyDriftCorrection(positions, W_InterpolatedDriftX, W_InterpolatedDriftY)
			Duplicate /O CorrectedPositions, root:'Localized Positions':$outputName
			
			DoWindow /K DriftInterpolatePanel
			
			DFREF savDF = GetDataFolderDFR()
			SetDataFolder driftInterpFolder
			KillWaves /A/Z
			SetDataFolder savDF
			KillDataFolder /Z driftInterpFolder
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function RecalculateDriftInterpolation()
	
	DFREF packageFolder = root:Packages:Localizer
	DFREF driftInterpFolder = packageFolder:DriftInterpolation
	
	wave /WAVE W_DriftInterpolationWaveRefs = driftInterpFolder:W_DriftInterpolationWaveRefs	// holds wave refs to positions and drift waves
	wave positions = W_DriftInterpolationWaveRefs[0]
	wave W_EstimatedDriftLocations = W_DriftInterpolationWaveRefs[1]
	wave W_EstimatedDriftX = W_DriftInterpolationWaveRefs[2]
	wave W_EstimatedDriftY = W_DriftInterpolationWaveRefs[3]
	
	Assert((DimSize(W_EstimatedDriftLocations, 0) == DimSize(W_EstimatedDriftX, 0)) && (DimSize(W_EstimatedDriftLocations, 0) == DimSize(W_EstimatedDriftY, 0)))
	
	variable nPositions = DimSize(positions, 0)
	variable firstFrame = positions[0][0]
	variable lastFrame = positions[nPositions - 1][0]
	variable nFrames = lastFrame - firstFrame + 1
	
	variable nDriftEstimates = DimSize(W_EstimatedDriftLocations, 0)
	variable firstFrameWithDriftEstimate = W_EstimatedDriftLocations[0]
	variable lastFrameWithDriftEstimate = W_EstimatedDriftLocations[nDriftEstimates - 1]
	
	Make /O/N=(nFrames) /D driftInterpFolder:W_InterpolatedDriftX, driftInterpFolder:W_InterpolatedDriftX
	wave W_InterpolatedDriftX = driftInterpFolder:W_InterpolatedDriftX
	wave W_InterpolatedDriftY = driftInterpFolder:W_InterpolatedDriftY
	SetScale /P x, firstFrame, 1, W_InterpolatedDriftX, W_InterpolatedDriftY
	
	// get the selected method to use
	ControlInfo /W=DriftInterpolatePanel#DriftInterpolateControlPanel TBInterpolationMode
	variable currentTab = V_value
	
	Switch (currentTab)
		case kDriftCorrLineTab:
			W_InterpolatedDriftX = Interp(clip(x, firstFrameWithDriftEstimate, lastFrameWithDriftEstimate), W_EstimatedDriftLocations, W_EstimatedDriftX)
			W_InterpolatedDriftY = Interp(clip(x, firstFrameWithDriftEstimate, lastFrameWithDriftEstimate), W_EstimatedDriftLocations, W_EstimatedDriftY)
			break
			
		case kDriftCorrPolyFitTab:
			NVAR gDriftInterpPolyFitXOrder = driftInterpFolder:V_DriftInterpPolyFitXOrder
			NVAR gDriftInterpPolyFitYOrder = driftInterpFolder:V_DriftInterpPolyFitYOrder
			
			variable polyFitXDegree = gDriftInterpPolyFitXOrder + 1	// because CurveFit expects e.g. 2 for a line
			variable polyFitYDegree = gDriftInterpPolyFitYOrder + 1	// because CurveFit expects e.g. 2 for a line
			
			// line needs to be handled separately
			if (polyFitXDegree == 2)
				CurveFit /Q/W=2 /N=1 line, W_EstimatedDriftX /X=W_EstimatedDriftLocations
			else
				CurveFit /Q/W=2 /N=1 poly polyFitXDegree, W_EstimatedDriftX /X=W_EstimatedDriftLocations
			endif
			wave W_Coef
			Duplicate /FREE W_Coef, W_PolyFit_X
			
			if (polyFitYDegree == 2)
				CurveFit /Q/W=2 /N=1 line, W_EstimatedDriftY /X=W_EstimatedDriftLocations
			else
				CurveFit /Q/W=2 /N=1 poly polyFitYDegree, W_EstimatedDriftY /X=W_EstimatedDriftLocations
			endif
			wave W_Coef
			Duplicate /FREE W_Coef, W_PolyFit_Y
			
			W_InterpolatedDriftX = Poly(W_PolyFit_X, x)
			W_InterpolatedDriftY = Poly(W_PolyFit_Y, x)
			
			break
			
		case kDriftCorrSplineTab:
			NVAR gDriftInterpSplineSmoothX = driftInterpFolder:V_DriftInterpSplineSmoothX
			NVAR gDriftInterpSplineSmoothY = driftInterpFolder:V_DriftInterpSplineSmoothY
			NVAR gDriftInterpSplineStdDevX = driftInterpFolder:V_DriftInterpSplineStdDevX
			NVAR gDriftInterpSplineStdDevY = driftInterpFolder:V_DriftInterpSplineStdDevY
			
			Interpolate2 /T=3/I=3 /F=(gDriftInterpSplineSmoothX)/S=(gDriftInterpSplineStdDevX) /Y=W_InterpolatedDriftX W_EstimatedDriftLocations, W_EstimatedDriftX
			Interpolate2 /T=3/I=3 /F=(gDriftInterpSplineSmoothY)/S=(gDriftInterpSplineStdDevY) /Y=W_InterpolatedDriftY W_EstimatedDriftLocations, W_EstimatedDriftY
			break
		default:
			Abort "Unknown tab in TBDriftInterpolatePanelProc()"
	EndSwitch
	
	// update the positions to correct with the new estimates
	wave M_TentativelyCorrectedPositions = driftInterpFolder:M_TentativelyCorrectedPositions
	wave correctedPositions = ApplyDriftCorrection(positions, W_InterpolatedDriftX, W_InterpolatedDriftY)
	Duplicate /O correctedPositions, M_TentativelyCorrectedPositions
	
	DoUpdate
End

Function RecoverEmittersByBleachingSteps(maxJumpDistance, windowName)
	variable maxJumpDistance
	string windowName
	
	string baseWindowName = GetBaseWindowName(windowName)
	DFREF windowDataFolder = GetWindowDataFolder(baseWindowName)
	
	NewDataFolder /O root:'Localized Positions'
	DFREF savDF = GetDataFolderDFR()
	
	Struct FitData fp
	// get all the fit options and parameters the user defined in the GUI
	GetGUIFitSettingsAndOptions(fp, baseWindowName)
	
	// Currently only symmetric 2D Gaussian fitting is supported
	if ((fp.localizationMethod != LOCALIZATION_GAUSS_FITTING) && (fp.localizationMethod != LOCALIZATION_GAUSS_FITTING_FIX))
		Abort "Currently only 2D symmetric Gaussian fitting is supported for the bleaching analysis"
	endif
	
	variable previousLoadedImage = fp.currentImage
	variable i, j, nPositionsFound, nPositionsOnWaitingList, nearestX, nearestY, closestIndex
	variable xCol, yCol, zCol, intensityCol, widthCol, dummyCol
	variable firstRun = 1
	variable offset, err
	
	for (i = fp.numberOfImages - 1; i >= 0; i-=1)
		// provide a progress window
		LocalizerProgressWindow("Bleaching", fp.numberOfImages - 1 - i, fp.numberOfImages - 1)
		
		ChangeCurrentImage(i, baseWindowName)
		wave Viewer_Temp = windowDataFolder:M_CCDImage
		
		 if (firstRun == 1)
			// set up a wave that will contain the subtracted frame
			Duplicate /O Viewer_Temp, windowDataFolder:M_Subtracted
			wave M_Subtracted = windowDataFolder:M_Subtracted
			Redimension /D M_Subtracted
			
			// set up a wave that will contain the contributions of all emitters that have been identified
			Duplicate /O Viewer_Temp, windowDataFolder:M_SummedEmitters
			Wave M_SummedEmitters = windowDataFolder:M_SummedEmitters
			Redimension /D M_SummedEmitters
			FastOP M_SummedEmitters = 0
		endif
		
		SetDataFolder windowDataFolder
		
		// subtract all emitters that have been localized in previous frames
		MatrixOP /O M_Subtracted = Viewer_Temp - M_SummedEmitters
		
		LocalizationAnalysis /M=(fp.localizationMethod) /D=(fp.thresholdMethod) /Y=(CAMERA_TYPE_IGOR_WAVE) /G={fp.preprocessing, fp.postprocessing} /F=(fp.particlefinder) /PFA=(fp.PFA) /T=(fp.directThresholdLevel) /R=(fp.minDistanceBetweenParticles) /W=(fp.PSFWidth) /PVER={fp.PVerSymm, fp.PVerEllipseSymm, fp.PVerOverlap} /Z /Q /DEST=windowDataFolder:POS_out "M_Subtracted"
		SetDataFolder savDF
		
		if (V_flag != 0)	// an error occurred during the fitting
			string errorMessage = "An error occurred while fitting frame " + num2str(i)
			Abort errorMessage
		endif
		
		wave POS_out = windowDataFolder:POS_out
		 nPositionsFound = DimSize(POS_out, 0)
		 
		 if (firstRun == 1)
			// make the output wave for the positions
			Duplicate /O POS_out, root:'Localized Positions':BleachingAnalysis
			wave M_BleachingAnalysis = root:'Localized Positions':BleachingAnalysis
			Redimension /N=(0, DimSize(POS_out, 1)) M_BleachingAnalysis
			
			Duplicate /O POS_out, windowDataFolder:M_WaitingList			// a wave that will contain a list of points encountered in the previous frame
			wave M_WaitingList = windowDataFolder:M_WaitingList			// the reason for its existence is that the frame in which we first encounter an emitter may not show its correct intensity since the emitter bleaches somewhere in that frame
			Redimension /N=(0, DimSize(POS_out, 1)) M_WaitingList			// so when a new point is found, add it to the waiting list first, and then only include it in the analysis if the next frame shows it again
			
			// get some indices into the columns of the positions wave
			GetColumnsForEmitterPositions(POS_out, xCol, yCol, zCol)
			GetColumnForIntegratedIntensity(POS_out, intensityCol)
			GetColumnsForFittedWidth(POS_out, widthCol, dummyCol)
			
			// if no width or amplitude is provided then abort
			if ((intensityCol < 0) || (widthCol < 0))
				Abort "The bleaching analysis requires a fitting method that provides amplitude and width information (e.g. Gaussian fitting)"
			endif
			
			// set up some viewer windows
			DoWindow /F BleachingAnalysis_Original
			if (V_flag == 0)
				Display /K=1 /N=BleachingAnalysis_Original ; AppendImage Viewer_Temp
				AppendToGraph /W=BleachingAnalysis_Original M_BleachingAnalysis[][yCol] vs M_BleachingAnalysis[][xCol]
				ModifyGraph /W=BleachingAnalysis_Original mode=3
			endif
			DoWindow /F BleachingAnalysis_Subtracted
			if (V_flag == 0)
				Display /K=1 /N=BleachingAnalysis_Subtracted ; AppendImage M_Subtracted
				AutoPositionWindow /M=0 /R=BleachingAnalysis_Original BleachingAnalysis_Subtracted
			endif
			
			ModifyGraph /W=BleachingAnalysis_Original width={Plan,1,bottom,left}
			ModifyGraph /W=BleachingAnalysis_Subtracted width={Plan,1,bottom,left}
			
			firstRun = 0
		endif
		 
		 if (i != 0)	// for all frames except the first, work with a waiting list
			 // did we localize any emitters that are on the waiting list?
			 for (j = 0; j < DimSize(M_WaitingList, 0); j+=1)
			 	// did we localize this position again?
			 	err = returnPointNearestFitPositions(POS_Out, M_WaitingList[j][xCol], M_WaitingList[j][yCol], nearestX, nearestY, closestIndex)
			 	
			 	if ((err != 0) || sqrt((nearestX - M_WaitingList[j][xCol])^2 + (nearestY - M_WaitingList[j][yCol])^2) > maxJumpDistance * fp.PSFWidth)	// is it close enough to count as the same emitter?
			 		// delete this point from the waiting list, we can't reproduce it
			 		DeletePoints j, 1, M_WaitingList
			 		j -= 1
			 		continue
			 		
			 	else		// the point corresponds to the same emitter as the one on the waiting list
			 			// add it to the list of recovered emitters, and remove it from the fitted positions list (otherwise it will end up on the waiting list again)
			 		Redimension /N=(DimSize(M_BleachingAnalysis, 0) + 1, DimSize(M_BleachingAnalysis, 1)) M_BleachingAnalysis
			 		M_BleachingAnalysis[DimSize(M_BleachingAnalysis, 0) - 1][] = POS_out[closestIndex][q]
			 		M_BleachingAnalysis[DimSize(M_BleachingAnalysis, 0) - 1][0] = i
			 		
			 		// add its contribution
			 		AddEmitter(M_SummedEmitters, POS_out[closestIndex][intensityCol], POS_out[closestIndex][widthCol], POS_out[closestIndex][xCol], POS_out[closestIndex][yCol])
			 		
			 		DeletePoints closestIndex, 1, POS_out
			 		DeletePoints j, 1, M_WaitingList
			 		j -= 1
			 	endif
			 endfor
		else	// for the very first frame, don't use the waiting list but just add all fitted points to the output waves
			// requested by Rob 30102009
			
			offset = DimSize(M_BleachingAnalysis, 0)
			Redimension /N=(DimSize(M_BleachingAnalysis, 0) + DimSize(POS_Out, 0), DimSize(M_BleachingAnalysis, 1)) M_BleachingAnalysis
	 		
	 		for (j = 0; j < DimSize(POS_Out, 0); j+=1)
	 			M_BleachingAnalysis[j + offset][] = POS_Out[j][q]
	 			M_BleachingAnalysis[j + offset][0] = i
	 			AddEmitter(M_SummedEmitters, POS_out[j][intensityCol], POS_out[j][widthCol], POS_out[j][xCol], POS_out[j][yCol])
	 		endfor
			
		endif
		 
		 // all points that are in the fitting list now count as new potential candidates
		 // add them to the waiting list
		 if (DimSize(POS_out, 0) > 0)
			Redimension /N=(DimSize(POS_out, 0), DimSize(POS_out, 1)) M_WaitingList
			M_WaitingList = POS_out
		endif
		 
		 DoUpdate
	endfor
	
	// now reverse the order of the fitted positions so the positions are located according to incrementing frame number
	Duplicate /FREE/O M_BleachingAnalysis, ShiftedBleaching
	ShiftedBleaching = M_BleachingAnalysis[DimSize(M_BleachingAnalysis, 0) - 1 - p][q]
	M_BleachingAnalysis = ShiftedBleaching
	
	KillWaves /Z POS_out
	
	ChangeCurrentImage(previousLoadedImage, baseWindowName)
	
End

Function LBDeleteDoubleClickedRowProc(lba) : ListBoxControl
	STRUCT WMListboxAction &lba

	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave

	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 3: // double click
			
			DeletePoints /M=0 row, 1, listWave
			if (WaveExists(selWave))
				DeletePoints /M=0 row, 1, selWave
			endif
			break
		case 4: // cell selection
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
	endswitch

	return 0
End
