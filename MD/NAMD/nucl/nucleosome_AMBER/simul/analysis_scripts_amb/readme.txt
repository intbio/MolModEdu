Here is more information on using scripts to analyze nucleosme trajectory.




Step 1. Preparation
Nucleosome orientation and trajectory preparation
1prepare.sh

Step 2. Visualization
2fast_visual.sh - will produce a movie, and an image with overlayed states.
tweek fast_render/movie.tcl to change it

High quality visualization
has to be done separatly using scripts in render_scripts

you can run vmd for viewing different parts of trj by launching
scripts
view_nucl_md.sh
view_nucl_md_sol.sh
view_nucl_min_eq.sh
vmd -e view_nucl_md_hq.tcl

Step 3. Analysis part 1. Movement analysis.
Here we get RMSD, RMSF, B-factor, PCA analysis

3gen_mov.sh - see details in script, many comments.

Step 4. Ion distribution analisys

We need to characterize ion distribution.
RDF with DNA and nucleosome, volume map,
see what DNA grooves are occupied, residence times,
condensation.

4ion.sh

Step 5. Water

Step 6. DNA

Step 7. Interactions

Energetic analysis


