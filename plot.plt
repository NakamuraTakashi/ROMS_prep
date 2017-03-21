# plot.plt
#
#$ gnuplot plot.plt
#   or
#$ gnuplot
#gnuplot> load 'plot.plt'
#


#set term png
set view map
#set xrange[0:500]
#set yrange[0:500]
#set cbrange[-1:1]
set size ratio 1
set palette rgbformulae 22,13,-31

#do for [i=0:100]{
#  set label 1 sprintf("t = %d", i*10) at 80,105 left # ƒ‰ƒxƒ‹1‚ğİ’è
#  infile = sprintf('output01/u%03d.dat',i)
#  outfile = sprintf('figs01/u%03d.png',i)
#  set output outfile
#  splot infile matrix with pm3d
#}
splot 'check4.dat' matrix with pm3d

