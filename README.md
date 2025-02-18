# README

This git repository has been built using the `template-research-project`
template or is the `template-research project`.

## Project organization

The project structure follows the recommendation of [quick guide to organize
your computation biology project](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

1. `data` contains the raw data of your project. If the data is too large, do
   not commit it (and ideally, create a small script to automatically download
   the data). If the data is very very large, and you cannot even download it
   on your computer, try to have a small subsample that you can commit here,
   and symlink (on the cluster) the large data in this folder. *e.g.*, if you
   work on the whole set of complete genomes of the NCBI, create a small
   subsample:

   `data/NCBI-partial-dataset`

   And add, *on the server*, a symlink to the whole dataset such that the
   structure of the folder is:

   ```
   data/
	NCBI-partial-dataset
	NCBI-full-dataset
   ```

   Code such that your scripts work transparently locally on the small
   subsample data, and remotely on the full dataset (by using snakemake for
   instance).

2. `scripts` contains a series of folders, each corresponding to a step in the
   data analysis process. In those folders, add all scripts necessary to this
   step. That includes visualization. A concrete example: working on Hi-C
   data, the first step is always to normalize the data. I thus create a
   folder `normalization` in `scripts`. In it, I add a script called
   `normalize_counts.py` and another one called `plot_contact_counts.py`. The
   first is the script I use to normalize the data. The second is a script I
   use to visualize the raw or the normalize contact counts: it outputs a png
   file for each of my datasets which I can then add to my labnotebook (more
   below).

3. `results` contains the labnotebook. This notebook should either by an html
   file or something that can be compiled in an HTML file so that it is easily
   shareable with collaborators.o

4. `doc` contains any communication document you may be creating relating to
   this project.


Check out the original publication to understand the additional folder

## The electronic labnotebook

Keep track of all and any experiment your run, whether it fails or not! Do
this in an electronic labnotebook, so that you can easily share it with
collaborators. This is particularly useful if you are working remotely and/or
have video confcall meeting regurlarly.

This labnotebook will be built automatically, and the resulting HTML will be
pushed on the `gh-pages`.

- This repository's labnotebook is done using
  [sphinx](http://www.sphinx-doc.org/en/stable/), but anything that builds in
  HTML in `results/_build/html` using `make all` will work with the current
  github workflow setup.
- Try to write something everyday.
