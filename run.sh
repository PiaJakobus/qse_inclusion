#! /bin/sh
prefix="${prefix:-/usr/local/hdf5}"
H5TOOL="h5fc"               # The tool name
H5TOOL_BIN="${prefix}/bin/${H5TOOL}"   # The path of the tool binary
${H5TOOL_BIN} hdf5_utils.f95 linspace.f95 read_qse.f95 one_zone.f95 

rm *.o
rm *.mod 

