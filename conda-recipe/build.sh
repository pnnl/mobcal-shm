export C_INCLUDE_PATH=/usr/include/

make || echo

mkdir -p $PREFIX/bin
cp mobcal $PREFIX/bin/mobcal
cp mobcal_shm $PREFIX/bin/mobcal_shm
