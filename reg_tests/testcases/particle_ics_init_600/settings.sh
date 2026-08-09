# --- General settings
jorekmodel="600"
description="kinetic_main with ICS density-proportional impurity initialisation (B@10%) using model$jorekmodel"
options="with_vpar=.true. with_TiTe=.false. with_neutrals=.false. with_impurities=.false. with_refluid=.false."
mpitasks=2
num_threads=2
binaries="kinetic_main"
binaries_initial="jorek_model${jorekmodel}_1"
extra_restart="part_restart.h5"
requiredfiles="input acd96_h.dat ccd96_h.dat plt96_h.dat prb96_h.dat prc96_h.dat scd96_h.dat acd89_b.dat ccd89_b.dat plt89_b.dat prb89_b.dat prc89_b.dat scd89_b.dat"
extra_remote_files="acd96_h.dat ccd96_h.dat plt96_h.dat prb96_h.dat prc96_h.dat scd96_h.dat acd89_b.dat ccd89_b.dat plt89_b.dat prb89_b.dat prc89_b.dat scd89_b.dat part_restart.h5"
particle_example="kinetic_main"
particle_example_dir="particles/examples"


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_coord_tor=1 l_pol_domm=0 n_plane=1 n_period=1 n_coord_period=1  $options  || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                                                                || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                                                                || exit 1
    make cleanall                                                                                                         || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=1 n_coord_tor=1 l_pol_domm=0 n_plane=1 n_period=1 n_coord_period=1   $options  || exit 1
  make $compilopt $debugoptions kinetic_main                                                                              || exit 1
}


# --- Initial run: build equilibrium then test particle initialisation
function initial_run () {
  ${codedir}/util/setinput.sh input restart_particles=.f. restart=.f. nstep_n=10,10,10 tstep_n=10.,100.,1000.             || exit 1
  $MPIRUN 1 ./jorek_model${jorekmodel}_1 < input | tee logfile_initial                                                    || exit 1
  ${codedir}/util/setinput.sh input restart_particles=.f. restart=.t. nstep_n=10 tstep_n=10.                              || exit 1
  export OMP_NUM_THREADS=$num_threads                                                                                     || exit 1
  $MPIRUN $mpitasks ./kinetic_main < input | tee logfile_initial2                                                         || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart_particles=.t. restart=.t. nstep_n=2 tstep_n=10.                               || exit 1
  export OMP_NUM_THREADS=$num_threads                                                                                     || exit 1
  $MPIRUN $mpitasks ./kinetic_main < input | tee logfile                                                                  || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-5                                                                                           || exit 1
}
