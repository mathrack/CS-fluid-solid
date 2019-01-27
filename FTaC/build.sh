set -e
#
cd plot_valid
for i in all_diss_K2.py all_diss_K2_fluid.py all_tt_K2.py all_tt_K2_fluid.py 395_diss.py 395_diss_fluid.py 395_tt.py 395_tt_fluid.py 1020_tt.py
do
  echo $i
  python $i &
done
cd ..
#
cd plot_higher
for i in tt_chti2.py tt_cht2j.py diss_chti2.py diss_cht2j.py contour_tt.py contour_tt_rel_error.py contour_ratio_eps.py contour_rel_error.py diss_g1k1.py diss_g01k02.py
do
  echo $i
  python $i &
done
cd ..
#
cd appendix_re_1020/plots
for i in diss.py rxx.py
do
  echo $i
  python $i &
done
cd ../..
#
wait
#
pdflatex LES_fluid_solid.tex
pdflatex LES_fluid_solid.tex
bibtex LES_fluid_solid
bibtex LES_fluid_solid
pdflatex LES_fluid_solid.tex
pdflatex LES_fluid_solid.tex
pdflatex LES_fluid_solid.tex
evince LES_fluid_solid.pdf
