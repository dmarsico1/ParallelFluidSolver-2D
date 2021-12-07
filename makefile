.PHONY: clean

fluid_par: grid.o mpi_interface.o netcdf_write.o therm.o init.o advect.o diffuse.o matrix.o poisson.o rk_time_step.o step.o fluid_par.o
	mpifort -o fluid_par grid.o mpi_interface.o netcdf_write.o therm.o init.o advect.o diffuse.o matrix.o poisson.o rk_time_step.o step.o fluid_par.o \
	-I/usr/include -L/usr/lib/x86_64-gnu-linux/ -lnetcdff -lnetcdf

grid.o: grid.f95
	mpifort -c grid.f95

mpi_interface.o: mpi_interface.f95
	mpifort -c mpi_interface.f95

netcdf_write.o: netcdf_write.f95
	mpifort -c netcdf_write.f95 -I /usr/include -L/usr/lib/x86_64-gnu-linux/ -lnetcdff -lnetcdf

therm.o: therm.f95
	mpifort -c therm.f95

init.o: init.f95
	mpifort -c init.f95

advect.o: advect.f95
	mpifort -c advect.f95

diffuse.o: diffuse.f95
	mpifort -c diffuse.f95

matrix.o: matrix.f95
	mpifort -c matrix.f95

poisson.o: poisson.f95
	mpifort -c poisson.f95

rk_time_step.o: rk_time_step.f95
	mpifort -c rk_time_step.f95

step.o: step.f95
	mpifort -c step.f95

fluid_par.o: fluid_par.f95
	mpifort -c fluid_par.f95


clean:
		rm -f grid.o mpi_interface.o netcdf_write.o therm.o init.o advect.o diffuse.o matrix.o poisson.o rk_time_step.o step.o fluid_par.o fluid_par \
			grid.mod mpi_interface.mod netcdf_write.mod therm.mod init.mod advect.mod diffuse.mod matrix.mod poisson.mod rk_time_step.mod step.mod

