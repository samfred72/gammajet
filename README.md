Step 1: replace all file paths with samson72 to the top of this git repo
  find . -type f -exec sed -i 's|/home/samson72/sphnx/gammajet|/your/path/to/this/directory|g' {} +
Step 2: unzip ttrees and put in 'trees' directory
Step 3: go to src/ and open 'make.sh', edit the file path to your install dir
  then: bash make.sh
Step 4: Make the histograms
  cd macros
  bash runall.sh
Step 5: Process the multijet
  (in macros dir)
  bash run_process_multijet.sh
Step 6: Run the grid scan
  (in macros dir)
  bash run_grid.sh
Step 7: Draw things
  cd ~/drawing
  root draw_all.C
  root draw_insitu_fit.C
  root draw_global_systematics.C
  root draw_multijet.C
