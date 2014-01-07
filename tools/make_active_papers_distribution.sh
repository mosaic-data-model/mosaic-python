rm -rf temp
mkdir temp
cd temp
mkdir code
mkdir code/python-packages
mkdir code/tests
mkdir documentation

cp -rp ../../lib/mosaic code/python-packages
cp -rp ../../lib/mosaic_pdb code/python-packages
cp -rp ../../tests/*_tests.py code/tests
rm code/tests/all_tests.py  # useless
rm code/tests/pdb_test.py   # requires an external file
for fn in `ls code/tests/*.py`; do
    python ../fix_test_runner.py $fn
done
find . -name \*.py[co] -exec rm \{\} \; 
find . -name __pycache__ -exec rm -r \{\} \; 

cp -rp ../../doc/source/*.rst documentation
cp -p ../../LICENSE.txt documentation
cp -p ../..//README.txt documentation
cat <<EOF >>documentation/README.txt

=============================================
Note about the ActivePapers version of Mosaic
=============================================

ActivePapers libraries require no installation
and profit from automatic dependency handling.
Using this ActivePaper requires no more than
a working ActivePapers installation.

This ActivePaper contains only the Mosaic Python
library and the documentation. The full pyMosaic
distribution also contains a command-line tool for
file conversion.  This tool make no sense in
the ActivePapers framework and was thus removed.

EOF

aptool -p pyMosaic.ap create
aptool checkin -t text documentation/*
aptool checkin -t module code/python-packages/*
aptool checkin -t importlet code/tests/*
aptool ln doi:10.6084/m9.figshare.692144: code/python-packages/immutable

mkdir -p ../../build
mv pyMosaic.ap ../../build

cd ..
rm -rf temp
