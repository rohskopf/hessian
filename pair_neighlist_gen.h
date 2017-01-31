/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(neighlist_gen,PairNeighlist_Gen)

#else

#ifndef LMP_PAIR_NEIGHLIST_GEN_H
#define LMP_PAIR_NEIGHLIST_GEN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNeighlist_Gen : public Pair {
 public:
  PairNeighlist_Gen(class LAMMPS *);
  virtual ~PairNeighlist_Gen();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_global;
  double **cut;
  double **cutmin;
  double **k,**r0;
  double **offset;
  double *iilist;                // i-i force constants from i-i list file
  double *ijlist;                // i-j force constants from i-j list file
  double **poseq;           // register the initial position of atoms -> consider as equilibrium position

  char **elements;              // names of unique elements
  int *map;                     // mapping from atom types to elements
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  int maxshort;                 // size of short neighbor list array
  int *neighshort;              // short neighbor list array

  void allocate();
  //virtual void read_second_order(char *);
  virtual void build_poseq(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
