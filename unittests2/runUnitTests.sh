#! /bin/bash
mkdir unittest_results_surface_core
./UnitedAtom --config-file=unittest_surface_core.config
mkdir unittest_results_surface
./UnitedAtom --config-file=unittest_surface.config
rm unittest_error_out.txt
for i in {1..20}
do
    ./UnitedAtom --config-file=unittest_error_estimate.config | tail -n 1 &>> unittest_error_out.txt
done
python compareUnitTests
