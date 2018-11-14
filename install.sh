#!/bin/bash
##############################
##### Date:21/05/2017
##### This file is used to install FuSeq from source codes. It will automatically Sailfish and its dependences.
##############################
CMAKECMD="cmake"
if [ "$DFETCH_BOOST" != "" ]; then 
    #echo "DFETCH_BOOST"
    CMAKECMD=$(echo $CMAKECMD" -DFETCH_BOOST="$DFETCH_BOOST)
fi;
if [ "$DBOOST_ROOT" != "" ]; then 
    #echo "DBOOST_ROOT" 
    CMAKECMD=$(echo $CMAKECMD" -DBOOST_ROOT="$DBOOST_ROOT)
fi;
if [ "$DTBB_INSTALL_DIR" != "" ]; then 
    #echo "DTBB_INSTALL_DIR"
    CMAKECMD=$(echo $CMAKECMD" -DTBB_INSTALL_DIR="$DTBB_INSTALL_DIR)
fi;
if [ "$DCMAKE_INSTALL_PREFIX" != "" ]; then 
    #echo "DCMAKE_INSTALL_PREFIX" 
    if [[ "$DCMAKE_INSTALL_PREFIX" != /* ]]; then
    	DCMAKE_INSTALL_PREFIX=$(echo $PWD"/"$DCMAKE_INSTALL_PREFIX)
	fi;
    CMAKECMD=$(echo $CMAKECMD" -DCMAKE_INSTALL_PREFIX="$DCMAKE_INSTALL_PREFIX)
fi;

#echo $CMAKECMD
#download sailfish
wget https://github.com/kingsfordgroup/sailfish/archive/v0.10.0.tar.gz
tar -xzvf v0.10.0.tar.gz
#copy our files 
cp -r include sailfish-0.10.0/
cp -r src sailfish-0.10.0/
#move to sailfish directory
cd sailfish-0.10.0
#do cmake
eval $CMAKECMD
make
make install
cd ..
#Configuration: set the current path for loading R functions
istr="/path/to"
ostr="$PWD/R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/FuSeq.R"
echo "Done!"
##### done



