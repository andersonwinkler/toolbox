SetZoomLevel 1
RedrawScreen

# Hide cursor
SetDisplayFlag 3 0

# Hide axes
SetDisplayFlag 23 0 

set outpth "$env(SUBJECTS_DIR)/$env(SUBJECT_NAME)/shots"

## Surfaces =================================================
LoadMainSurface 0 lh.white
LoadMainSurface 1 rh.white
SetSurfaceLineWidth 0 0 2
SetSurfaceLineWidth 0 2 2
SetSurfaceLineWidth 1 0 2
SetSurfaceLineWidth 1 2 2
SetDisplayFlag 5 0

# Coronal
SetOrientation 0
SetCursor 0 128 128 068
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor068.tif"
SaveTIFF $tiff
SetCursor 0 128 128 068
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor068.tif"
SaveTIFF $tiff
SetCursor 0 128 128 088
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor088.tif"
SaveTIFF $tiff
SetCursor 0 128 128 098
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor098.tif"
SaveTIFF $tiff
SetCursor 0 128 128 108
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor108.tif"
SaveTIFF $tiff
SetCursor 0 128 128 118
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor118.tif"
SaveTIFF $tiff
SetCursor 0 128 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor128.tif"
SaveTIFF $tiff
SetCursor 0 128 128 138
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor138.tif"
SaveTIFF $tiff
SetCursor 0 128 128 148
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor148.tif"
SaveTIFF $tiff
SetCursor 0 128 128 158
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor158.tif"
SaveTIFF $tiff
SetCursor 0 128 128 168
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor168.tif"
SaveTIFF $tiff
SetCursor 0 128 128 188
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_cor188.tif"
SaveTIFF $tiff

# Transversal
SetOrientation 1
SetCursor 0 128 118 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_tra118.tif"
SaveTIFF $tiff
SetCursor 0 128 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_tra128.tif"
SaveTIFF $tiff
SetCursor 0 128 138 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_tra138.tif"
SaveTIFF $tiff

# Sagittal
SetOrientation 2
SetCursor 0 88 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_sag088.tif"
SaveTIFF $tiff
SetCursor 0 168 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_surfaces_sag168.tif"
SaveTIFF $tiff
UnloadAllSurfaces

## Subcortical =================================================
LoadSegmentationVolume 0 aseg.mgz $env(FREESURFER_HOME)/FreeSurferColorLUT.txt

# Coronal
SetOrientation 0
SetCursor 0 128 128 068
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor068.tif"
SaveTIFF $tiff
SetCursor 0 128 128 088
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor088.tif"
SaveTIFF $tiff
SetCursor 0 128 128 098
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor098.tif"
SaveTIFF $tiff
SetCursor 0 128 128 108
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor108.tif"
SaveTIFF $tiff
SetCursor 0 128 128 118
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor118.tif"
SaveTIFF $tiff
SetCursor 0 128 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor128.tif"
SaveTIFF $tiff
SetCursor 0 128 128 138
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor138.tif"
SaveTIFF $tiff
SetCursor 0 128 128 148
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor148.tif"
SaveTIFF $tiff
SetCursor 0 128 128 158
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor158.tif"
SaveTIFF $tiff
SetCursor 0 128 128 168
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor168.tif"
SaveTIFF $tiff
SetCursor 0 128 128 188
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_cor188.tif"
SaveTIFF $tiff

# Transversal
SetOrientation 1
SetCursor 0 128 118 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_tra118.tif"
SaveTIFF $tiff
SetCursor 0 128 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_tra128.tif"
SaveTIFF $tiff
SetCursor 0 128 138 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_tra138.tif"
SaveTIFF $tiff

# Sagittal
SetOrientation 2
SetCursor 0 88 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_sag088.tif"
SaveTIFF $tiff
SetCursor 0 168 128 128
RedrawScreen
set tiff "$outpth/$env(SUBJECT_NAME)_aseg_sag168.tif"
SaveTIFF $tiff

exit
