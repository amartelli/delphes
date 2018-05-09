./configure
sed -i -e 's/c++0x/c++1y/g' Makefile

export PYTHIA8=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/pythia8/223-mlhled2/
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH

make -j 4 HAS_PYTHIA8=true DelphesPythia8  
