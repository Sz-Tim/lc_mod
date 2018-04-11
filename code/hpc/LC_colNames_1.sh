#!/bin/bash

FILES=out/full_Clim*.csv
orig_csv=($FILES)
orig_csv_full=($(basename -a $FILES))

for f in {0..23}
do
	orig_full_f="${orig_csv_full[$f]}"
	orig_f_name="${orig_full_f%.csv}"
	outf=$orig_f_name
	outf+=".csv"
	sed '/^#/ d' out/$orig_full_f | head -1 > full_colNames/$outf
done