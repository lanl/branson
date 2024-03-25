```
CALIPER_ROOT_DIR=<path>/Caliper-Adiak/install CC=cc CXX=CC cmake   -Dadiak_DIR=<path>/Adiak/install/lib/cmake/adiak/ -Dcaliper_DIR=<path>/Caliper-Adiak/install/share/cmake/caliper/  -DMETIS_ROOT_DIR==<path>/new/metis/install -DCMAKE_BUILD_TYPE=Release  -DWITH_ADIAK=ON  ../src 
```
