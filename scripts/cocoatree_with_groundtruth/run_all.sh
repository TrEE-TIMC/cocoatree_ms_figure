python run_cocoatree_with_gt.py halabi -o results/cocoatree_gt/halabi
python run_cocoatree_with_gt.py rivoire -o results/cocoatree_gt/rivoire
python run_cocoatree_with_gt.py rhomboid -o results/cocoatree_gt/rhomboid

# Now plot results
python plot_sectors.py results/cocoatree_gt/rhomboid/cocoatree_SCA_none.csv \
    -o images/cocoatree_with_gt/rhomboid_sectors.png
python plot_sectors.py results/cocoatree_gt/halabi/cocoatree_SCA_none.csv \
    -o images/cocoatree_with_gt/halabi_sectors.png
python plot_sectors.py results/cocoatree_gt/rivoire/cocoatree_SCA_none.csv \
    -o images/cocoatree_with_gt/rivoire_sectors.png
