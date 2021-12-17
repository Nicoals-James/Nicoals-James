/*
 * File: plot_model.tcl
 * Created Date: 2021-12-14 22:52:43
 * Author: Nicolas zhan
 * Contact: 1059291645@qq.com
 * -----
 * Last Modified: 2021-12-14 22:52:45
 */





###################
#Author: yaoyl
#Last modified: 2021\12\13
#Plot the model of the model xcm and write as txt file
###################

namespace import ::tcl::mathfunc::cos
namespace import ::tcl::mathfunc::sqrt
query y
#lmod relxill ~/relxill2/
lmod relxill /data/relxill2/

#Input of para names and para values

# set name [lindex $argv 0]
# set values [lindex $argv 1]

#Iput of test_num and xcm_names
set test_num [lindex $argv 0]
set xcm_names [lindex $argv 1]
puts $test_num
puts $xcm_names

#At the model xcm
#As an example constants*tbabs*diskbb*relxillcp*xillvercp
#@ctdrx.xcm

@$xcm_names


#plot the model
# tclout modpar
# set num_pars $xspec_tclout
# #The order of the input paras must be in order
# set j 0
# for { set i 1}  {$i < $num_pars+1} {incr i} {
#    tclout pinfo $i
#    if { $xspec_tclout == $name }{
#    	newpar i $value -1 $value $value $value $value
#    }
# }


#Here choose the type of x axis
#set plot energy
#set plot channel

cpd $test_num\_model_plot.gif/GIF
#plot eemodel
#Get the x values
tclout plot model x
set x_list $xspec_tclout
#Get the y values
tclout plot model y
set y_list $xspec_tclout
cpd none

set fp [open "$test_num\_model.txt" w+]
for {set i 0} {$i < [llength $x_list]} {incr i} {
	set x_value [lindex $x_list $i]
	set y_value [lindex $y_list $i]
	puts $fp "$x_value $y_value"
}
close $fp

rm $test_num\_model_plot.gif
file rename $test_num\_model_plot.gif_3 $test_num\_model_plot.gif
rm *gif_*
