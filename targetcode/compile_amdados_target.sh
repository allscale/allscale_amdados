#/bin/bash

# cross-compile the amdados generated code from the AllScale compiler 
# for Power system

CC=g++-7
AllScale_directory="/home/fearghal/Fearghal/Allscale/"
HPX_build_dir="/home/fearghal/Fearghal/Allscale/hpx_release/hpx/"
HPX_install_dir="/home/fearghal/local/hpx_release/"
Boost_install="/home/fearghal/local/boost_1_67/"
echo $AllScale_directory"allscale_runtime"
target="amdados_generated"
echo $Boost_install"lib"  
$CC $target".cpp" -o amdados_cc -Bdynamic -I$AllScale_directory"allscale_runtime" -I$HPX_install_dir -I$HPX_install_dir"include" -I$HPX_build_dir -I$AllScale_directory"allscale_api/code/api/include" -I$AllScale_directory"allscale_api/code/utils/include" -I$Boost_install"include" -I/home/fearghal/hwloc/include --std=c++14 -L$Boost_install"lib" -L$AllScale_directory"allscale_runtime/build/src" -L$HPX_build_dir"my_hpx_build/lib" -L/home/fearghal/local/hwloc/lib -lm -lpthread -O3 -DALLSCALE_WITH_HPX -DHPX_DISABLE_ASSERTS -DBOOST_DISABLE_ASSERTS -lhpx_allscale -lhpx_init -lhpx -lboost_chrono -lboost_date_time -lboost_filesystem -lboost_program_options -lboost_regex -lboost_system -lboost_thread -lboost_atomic -lrt -pipe -fshow-column -fpermissive
