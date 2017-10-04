# this file to be read after running something like 
# tksurfer JS_071026js rh inflated -curv -tcl scriptStartLHAll.tcl

set hemiSphere "---HEMI---"
set subject "---SJ---"
set fs_subject "---FS_SJ---"
set polecc "---POLECC---"
set condition "---CONDITION---"
set name "---NAME---"
set exit_when_ready "---EXIT---"
set mask "---MASK---"

if {$hemiSphere == "rh"} {
	# for right hemisphere we need the following parameters
	set yDirection -1
} else {		
# for left hemisphere we need the following parameters
	set yDirection 1
}

#set polardir [format "%s/surf" $condition]
set polardir $subject

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

### for backward compatibility (old script-specific mechanism)
set floatstem $polecc                   ;# float file stem
set realname real                      ;# analyse infix
set complexname imag                   ;# analyse infix
set rgbname polar                   ;# name of rgbfiles

#### parm defaults: can reset in csh script with setenv
puts "tksurfer: [file tail $script]: read and smooth complex Fourier comp"
set overlayflag 1       ;# overlay data on gray brain
set surfcolor 1         ;# draw the curvature under data
set avgflag 1           ;# make half convex/concave
set complexvalflag 1    ;# two-component data
set colscale 0          ;# 0=wheel,1=heat,2=BR,3=BGR,4=twocondGR,5=gray
set angle_offset -.125   ;# phase offset (-0.25 for up semicircle start)
set angle_cycles 1.0    ;# adjust range
set fthresh 0.05       ;# val/curv sigmoid zero (neg=>0)
set fslope 1         	;# contast (was fsquash 2.5)
set fmid   0.1         ;# set linear region
set smoothsteps 0
set offset 0.20    		;# default lighting offset
# smooth the curvature and surface before doing anything else
set invphaseflag 0
set revphaseflag 0

if { [info exists revpolarflag] } { 
  set revphaseflag $revpolarflag 
}

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
read_binary_curv

#### setenv polardir overrides setenv dir
if [info exists polardir] { set dir $polardir }

#### read and smooth complex component MRI Fourier transform of data
puts "tksurfer: [file tail $script]: read and smooth complex Fourier comp"
setfile val */$dir/surf/${condition}/---POLECC_S---_${complexname}_${mask}-$hemi.$fs_subject.mgz     ;# polarangle
read_binary_values
smooth_val $smoothsteps 
shift_values     ;# shift complex component out of way

#### read and smooth real component MRI Fourier transform of data
puts "tksurfer: [file tail $script]: read and smooth real Fourier comp"
setfile val */$dir/surf/${condition}/---POLECC_S---_${realname}_${mask}-$hemi.$fs_subject.mgz    ;# polarangle
read_binary_values
smooth_val $smoothsteps

#### scale and position brain
puts "tksurfer: [file tail $script]: scale, position brain"
open_window
make_lateral_view       ;# rotate either hemisphere
do_lighting_model -1 -1 -1 -1 $offset ;# -1 => nochange; diffuse curv (def=0.15)

# setup overlay characteristics
set gaLinkedVars(fthresh) $fthresh
SendLinkedVarGroup overlay
set gaLinkedVars(fmid) $fmid
SendLinkedVarGroup overlay
set gaLinkedVars(fslope) $fslope
SendLinkedVarGroup overlay

scale_brain 1.8
set nrimages 1
set rotation_gain 360.0
set rot [ expr { $rotation_gain / $nrimages } ]

for {set base_y_rotation 50} {$base_y_rotation < 140} {incr base_y_rotation 40} {
	make_lateral_view
	rotate_brain_y [ expr { ($base_y_rotation * $yDirection) } ]
	rotate_brain_x [ expr { -$rotation_gain / 2.0 } ]
	redraw
	for {set i 0} {$i < $nrimages} {incr i} {
		rotate_brain_x [ expr { $rot * $i } ]
		redraw
		
		set l [string length "$i"] 
		if {$l == 1} {
		    set label "00$i"
		} elseif {$l == 2} {
		    set label "0$i"
		} elseif {$l == 3} {
		    set label "$i"
		}
		
		set fN [format "---FIGPATH---/%s_%s_---POLECC_S---_%s-%s-%s-%s-%s.tiff" $mask $subject $condition $hemiSphere $label $base_y_rotation $fs_subject]
		save_tiff $fN
		rotate_brain_x [ expr { -$rot * $i } ]
	}
}
if {$exit_when_ready == 1} {
    exit
}
