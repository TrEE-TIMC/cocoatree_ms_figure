python plot_subfigure_panel.py halabi -m SCA -c none
python plot_subfigure_panel.py halabi -m MI -c none
python plot_subfigure_panel.py halabi -m NMI -c none
python plot_subfigure_panel.py halabi -m MI -c APC

python plot_subfigure_panel.py rhomboid -m SCA -c none
python plot_subfigure_panel.py rhomboid -m MI -c none
python plot_subfigure_panel.py rhomboid -m NMI -c none
python plot_subfigure_panel.py rhomboid -m MI -c APC

python plot_subfigure_panel.py DHFR -m SCA -c none
python plot_subfigure_panel.py DHFR -m MI -c none
python plot_subfigure_panel.py DHFR -m NMI -c none
python plot_subfigure_panel.py DHFR -m MI -c APC

python plot_figure_conservation_vs_cumscore.py halabi
python plot_figure_conservation_vs_cumscore.py rhomboid
python plot_figure_conservation_vs_cumscore.py DHFR
