// build_lstexcess.do -- compile the Mata library for stexcess.
// Run from the repo root:
//   /Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp -q -b do ado/build_lstexcess.do
// Produces ado/lstexcess.mlib (shipped alongside the .ado files).

clear all
do "ado/stexcess_mata.do"
mata: mata mlib create lstexcess, dir("ado") replace
mata: mata mlib add lstexcess _stx_*(), dir("ado")
mata: mata mlib index
di as txt "built ado/lstexcess.mlib"
