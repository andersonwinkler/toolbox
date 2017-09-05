# This TCL script is supposed to be called from within tksurfer, by the
# script called 'takeshots'.
# The environment variables that should be present are:
# - SUBJECTS_DIR: regular FS variable.
# - SUBJECT_NAME: name of the current subject.
# - SURF: surface to be loaded, e.g. pial, white, inflated, etc.
# - WHAT: parcellation scheme, e.g. aparc, a2005s, a2009s.
#
# _____________________________________
# Anderson M. Winkler
# Yale University / Institute of Living
# Feb/2010
# http://brainder.org

open_window

set lablpth "$env(SUBJECTS_DIR)/$env(SUBJECT_NAME)/label"
set outpth "$env(SUBJECTS_DIR)/$env(SUBJECT_NAME)/shots"

labl_import_annotation "$lablpth/$hemi.$env(WHAT).annot"
redraw
UpdateAndRedraw

puts "Taking Snapshots..."

make_lateral_view
redraw
set tiff "$outpth/$env(SUBJECT_NAME)_${hemi}_$env(SURF)_$env(WHAT)_lat.tif"
save_tiff $tiff

rotate_brain_x 90
redraw
set tiff "$outpth/$env(SUBJECT_NAME)_${hemi}_$env(SURF)_$env(WHAT)_inf.tif"
save_tiff $tiff

rotate_brain_x 180
redraw
set tiff "$outpth/$env(SUBJECT_NAME)_${hemi}_$env(SURF)_$env(WHAT)_sup.tif"
save_tiff $tiff

make_lateral_view
rotate_brain_y 180
redraw
set tiff "$outpth/$env(SUBJECT_NAME)_${hemi}_$env(SURF)_$env(WHAT)_med.tif"
save_tiff $tiff

exit
