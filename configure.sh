#!/bin/bash -l
##### Date:21/05/2017
##### configure.sh: set the current path for loading R functions
istr="/path/to"
ostr="$PWD/R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/FuSeq.R"
echo "Done!"