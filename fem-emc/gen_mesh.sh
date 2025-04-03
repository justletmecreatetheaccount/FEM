cd genmesh/build
cmake ..
make
./genmesh plot
cd ../../data
python3 fixMesh.py
cd ..