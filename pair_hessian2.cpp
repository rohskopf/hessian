/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_hessian2.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

#include <stdlib.h>
#include <string.h>
#include "create_box.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "region_prism.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
#define MAXLINE 1024

FILE * debug_read_files = fopen ("DEBUGHESS_READ","w");
FILE * debug_compute = fopen ("DEBUGHESS_COMPUTE","w");
//FILE * fh_neighlist = fopen ("NEIGHLIST","w");
//FILE * debug_poseq = fopen ("DEBUGHESS","w");

/* ---------------------------------------------------------------------- */

PairHessian2::PairHessian2(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  iilist = NULL;
  ijlist = NULL;
  poseq = NULL;
}

/* ---------------------------------------------------------------------- */

PairHessian2::~PairHessian2()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cutmin);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(offset);

    memory->destroy(iilist);
    memory->destroy(ijlist);
    memory->destroy(poseq);
  }
}

/* ---------------------------------------------------------------------- */

void PairHessian2::compute(int eflag, int vflag)
{
  //fprintf(debug_compute, "Beginning compute\n");
  int i,j,ii,jj,inum,jnum,itype,jtype;
  tagint itag,jtag,ktag;
  imageint *image = atom->image;
  double xtmp,ytmp,ztmp,delx,dely,delz,delxij,delyij,delzij,delxik,delyik,delzik,evdwl,fpair;
  double rsq,r,dr,rk,factor_lj;
  double poseqi[3],poseqj[3],deli[3],delj[3],fzeros[3], fpair2[3];//,fpair[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xhi, yhi, zhi;
  int *jlist_test;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int natoms = atom->natoms;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  xhi = domain->boxhi[0];
  yhi = domain->boxhi[1];
  zhi = domain->boxhi[2];
  int iipair_indx = 0;
  int ijpair_indx = 0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // unwrap contains the position of each atom in the supercell, useful to know the displacement for neighbor atoms
  double **unwrap;
  memory->create(unwrap,natoms,3,"domain:unwrap");
  for (i = 0; i < natoms; i++) {
    for (j = 0; j < 3; j++) {
      unwrap[i][j] = 0.0;  
    }
  }
 
  for (i = 0; i < nlocal; i++)
    domain->unmap(x[i],image[i],unwrap[tag[i]-1]);

  double unwloc[natoms];
  double unwall[natoms];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < natoms; j++) {
      unwloc[j] = unwrap[j][i];
    }
    MPI_Allreduce(&unwloc,&unwall,natoms,MPI_DOUBLE,MPI_SUM,world);
    for (j = 0; j < natoms; j++) {
      unwrap[j][i] = unwall[j];
    }
  }
/*  if (comm->me == 0) {
    for (i = 0; i < atom->natoms; i++) {
      printf("itag = %d %f %f %f\n",i+1,unwrap[i][0],unwrap[i][1],unwrap[i][2]);
    }
  }
*/

  //fprintf(debug_compute, "test\n");

  /*if (comm->me == 0) {
    for (int i = 0; i < natoms; i++) {
      fprintf(debug_compute,"%.12f %.12f %.12f\n",poseq[i][0],poseq[i][1],poseq[i][2]);
    }
  }*/

  // loop over neighbors of my atoms
  //fprintf (debug_comp, "test");
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    poseqi[0] = poseq[itag-1][0];
    poseqi[1] = poseq[itag-1][1];
    poseqi[2] = poseq[itag-1][2];

    deli[0] = unwrap[itag-1][0] - poseqi[0];
    deli[1] = unwrap[itag-1][1] - poseqi[1];
    deli[2] = unwrap[itag-1][2] - poseqi[2];
    /*deli[0] = 0;
    deli[1] = 0;
    deli[2] = 0;*/

    //fprintf(debug_compute, "itag: -------------------- %d\n", itag);
    //fprintf(debug_compute, "%.12f %.12f %.12f\n", deli[0], deli[1], deli[2]);
    // Take care of periodic boundary conditions
    /*if (deli[0] > 1) {
      deli[0] = (unwrap[itag-1][0] - xhi) - poseqi[0];
    } else if (deli[0] < -1) deli[0] = unwrap[itag-1][0] - (poseqi[0] - xhi);
    if (deli[1] > 1) {
      deli[1] = (unwrap[itag-1][1] - yhi) - poseqi[1];
    } else if (deli[1] < -1) deli[1] = unwrap[itag-1][1] - (poseqi[1] - yhi);
    if (deli[2] > 1) {
      deli[2] = (unwrap[itag-1][2] - zhi) - poseqi[2];
    } else if (deli[2] < -1) deli[2] = unwrap[itag-1][2] - (poseqi[2] - zhi);*/
    //fprintf(debug_compute, "%.12f %.12f %.12f\n", deli[0], deli[1], deli[2]);

    //fprintf(debug_compute, "%d\n", itag);
    //fprintf(debug_compute, "%f %f %f\n", poseqi[0], poseqi[1], poseqi[2]);
    //fprintf(debug_compute, "%.12f %.12f %.12f\n", deli[0], deli[1], deli[2]);
    //fprintf(debug_compute, "atom1: %d atom2: %d FCxx: %.12f\n", itag, itag, fc2[itag-1][0][itag-1][0]*(1/13.605698066)*(0.529177249)*(0.529177249));


    /*fprintf(debug_compute, "%f %f %f %f %f %f %f %f %f\n", iilist[0+iipair_indx], iilist[1+iipair_indx], iilist[2+iipair_indx], iilist[3+iipair_indx], iilist[4+iipair_indx], \
                                                           iilist[5+iipair_indx], iilist[6+iipair_indx], iilist[7+iipair_indx], iilist[8+iipair_indx]);*/

    //fprintf(debug_compute, "%f\n", iilist[0]);

    //fprintf(debug_compute, "%f\n", cutmin[0][0]);
    // Compute self interactions if cutmin == 0
    if (cutmin[0][0] < 0.01)
    {
      fpair2[0] = -1.0*(iilist[0+iipair_indx]*deli[0] + iilist[1+iipair_indx]*deli[1] + iilist[2+iipair_indx]*deli[2]);
      fpair2[1] = -1.0*(iilist[3+iipair_indx]*deli[0] + iilist[4+iipair_indx]*deli[1] + iilist[5+iipair_indx]*deli[2]);
      fpair2[2] = -1.0*(iilist[6+iipair_indx]*deli[0] + iilist[7+iipair_indx]*deli[1] + iilist[8+iipair_indx]*deli[2]);
      /*fpair2[0] = 0;
      fpair2[1] = 0;
      fpair2[2] = 0;*/
      //fprintf(debug_compute, "%.12f %.12f %.12f\n", fpair2[0], fpair2[1], fpair2[2]);

      f[i][0] += fpair2[0];
      f[i][1] += fpair2[1];
      f[i][2] += fpair2[2];

      if (eflag) {
        evdwl = 0.5*(iilist[0+iipair_indx]*deli[0]*deli[0] + iilist[1+iipair_indx]*deli[0]*deli[1] + iilist[2+iipair_indx]*deli[0]*deli[2] + \
                     iilist[3+iipair_indx]*deli[1]*deli[0] + iilist[4+iipair_indx]*deli[1]*deli[1] + iilist[5+iipair_indx]*deli[1]*deli[2] + \
                     iilist[6+iipair_indx]*deli[2]*deli[0] + iilist[7+iipair_indx]*deli[2]*deli[1] + iilist[8+iipair_indx]*deli[2]*deli[2]);
        //evdwl = 0;
        ev_tally(i,i,nlocal,newton_pair,evdwl,0.0,0.0,0.0,0.0,0.0);
      }
    }

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    iipair_indx = iipair_indx + 9;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      //jlist_test[jj] = j;
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      jtag = tag[j];
      //fprintf(debug_compute, "i: %d j: %d\n", i, j);
      //fprintf(debug_compute, "itag: %d jtag: %d\n", itag, jtag);

      poseqj[0] = poseq[jtag-1][0];
      poseqj[1] = poseq[jtag-1][1];
      poseqj[2] = poseq[jtag-1][2];

      delj[0] = unwrap[jtag-1][0] - poseqj[0];
      delj[1] = unwrap[jtag-1][1] - poseqj[1];
      delj[2] = unwrap[jtag-1][2] - poseqj[2];
      /*delj[0] = 0;
      delj[1] = 0;
      delj[2] = 0;*/

      // Take care of periodic boundary conditions
      /*if (delj[0] > 1) {
        delj[0] = (unwrap[jtag-1][0] - xhi) - poseqj[0];
      } else if (delj[0] < -1) delj[0] = unwrap[jtag-1][0] - (poseqj[0] - xhi);
      if (delj[1] > 1) {
        delj[1] = (unwrap[jtag-1][1] - yhi) - poseqj[1];
      } else if (delj[1] < -1) delj[1] = unwrap[jtag-1][1] - (poseqj[1] - yhi);
      if (delj[2] > 1) {
        delj[2] = (unwrap[jtag-1][2] - zhi) - poseqj[2];
      } else if (delj[2] < -1) delj[2] = unwrap[jtag-1][2] - (poseqj[2] - zhi);*/

      delxij = (xtmp - deli[0]) - (x[j][0] - delj[0]);
      delyij = (ytmp - deli[1]) - (x[j][1] - delj[1]);
      delzij = (ztmp - deli[2]) - (x[j][2] - delj[2]);

      //fprintf(debug_compute, "%d\n", itag);
      //fprintf(debug_compute, "%f %f %f\n", poseqi[0], poseqi[1], poseqi[2]);
      //fprintf(debug_compute, "%.12f %.12f %.12f\n", delj[0], delj[1], delj[2]);

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq > cutmin[itype][jtype] && rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        //fprintf(debug_compute, "2: %d\n", itag);
        /*dr = r - r0[itype][jtype];
        dexp = exp(-alpha[itype][jtype] * dr);
        fpair = factor_lj * morse1[itype][jtype] * (dexp*dexp - dexp) / r;*/

        //fprintf(debug_compute, "itag: %d jtag: %d\n", itag, jtag);
        //fprintf(debug_compute, "%f\n", fc2_temp2[itag-1][0][jtag-1][0]);
        //fprintf(debug_compute, "atom1: %d atom2: %d FCxx: %.12f\n", itag, jtag, fc2[itag-1][0][jtag-1][0]*(1/13.605698066)*(0.529177249)*(0.529177249));
        //fprintf(fh_neighlist, "%d %d %f %f %f\n", itag, jtag, delx, dely, delz);
        //fprintf(debug_compute, "%f %f %f\n", xtmp, ytmp, ztmp);
        //fprintf(debug_compute, "%f %f %f\n", x[j][0], x[j][1], x[j][2]);

        /*fprintf(debug_compute, "%f %f %f %f %f %f %f %f %f\n", ijlist[0+ijpair_indx], ijlist[1+ijpair_indx], ijlist[2+ijpair_indx], ijlist[3+ijpair_indx], ijlist[4+ijpair_indx], \
                                                               ijlist[5+ijpair_indx], ijlist[6+ijpair_indx], ijlist[7+ijpair_indx], ijlist[8+ijpair_indx]);*/

        /*fpair2[0] = -1.0*(fc2[itag-1][0][jtag-1][0]*delj[0] + fc2[itag-1][0][jtag-1][1]*delj[1] + fc2[itag-1][0][jtag-1][2]*delj[2]);
        fpair2[1] = -1.0*(fc2[itag-1][1][jtag-1][0]*delj[0] + fc2[itag-1][1][jtag-1][1]*delj[1] + fc2[itag-1][1][jtag-1][2]*delj[2]);
        fpair2[2] = -1.0*(fc2[itag-1][2][jtag-1][0]*delj[0] + fc2[itag-1][2][jtag-1][1]*delj[1] + fc2[itag-1][2][jtag-1][2]*delj[2]);*/
        /*fpair2[0] = 0;
        fpair2[1] = 0;
        fpair2[2] = 0;*/
        fpair2[0] = -1.0*(ijlist[0+ijpair_indx]*delj[0] + ijlist[1+ijpair_indx]*delj[1] + ijlist[2+ijpair_indx]*delj[2]);
        fpair2[1] = -1.0*(ijlist[3+ijpair_indx]*delj[0] + ijlist[4+ijpair_indx]*delj[1] + ijlist[5+ijpair_indx]*delj[2]);
        fpair2[2] = -1.0*(ijlist[6+ijpair_indx]*delj[0] + ijlist[7+ijpair_indx]*delj[1] + ijlist[8+ijpair_indx]*delj[2]);
        //fprintf(debug_compute, "%.12f %.12f %.12f\n", fpair2[0], fpair2[1], fpair2[2]);

        f[i][0] += fpair2[0];
        f[i][1] += fpair2[1];
        f[i][2] += fpair2[2];

        if (eflag) {
          evdwl = 0.5*(ijlist[0+ijpair_indx]*deli[0]*delj[0] + ijlist[1+ijpair_indx]*deli[0]*delj[1] + ijlist[2+ijpair_indx]*deli[0]*delj[2] + \
                       ijlist[3+ijpair_indx]*deli[1]*delj[0] + ijlist[4+ijpair_indx]*deli[1]*delj[1] + ijlist[5+ijpair_indx]*deli[1]*delj[2] + \
                       ijlist[6+ijpair_indx]*deli[2]*delj[0] + ijlist[7+ijpair_indx]*deli[2]*delj[1] + ijlist[8+ijpair_indx]*deli[2]*delj[2]);
          //evdwl = 0;
          ev_tally(i,i,nlocal,newton_pair,evdwl,0.0,0.0,0.0,0.0,0.0);
        }

        if (evflag) {
          ev_tally(i,i,nlocal,newton_pair,evdwl,0.0,0.0,delxij,delyij,delzij);
          //v_tally_tensor(i,j,nlocal,newton_pair,delxij*fpair[0], delyij*fpair[1], delzij*fpair[2],delxij*fpair[1],delxij*fpair[2],delyij*fpair[2]);
        }

        /*if (newton_pair || j < nlocal) {
          f[j][0] -= fpair2[0];
          f[j][1] -= fpair2[1];
          f[j][2] -= fpair2[2];
        }*/

        dr = r - r0[itype][jtype];
        rk = k[itype][jtype] * dr;
        fpair = factor_lj * -2.0*rk/r;

        /*f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = rk*dr - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);*/

        ijpair_indx = ijpair_indx + 9;
      }
    }
  }

  //memory->destroy(iilist);
  //memory->destroy(ijlist);
  memory->destroy(unwrap);
  //memory->destroy(poseq);

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairHessian2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  //memory->create(cutminsq,n+1,n+1,"pair:cutminsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cutmin,n+1,n+1,"pair:cutmin");
  memory->create(k,n+1,n+1,"pair:k");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairHessian2::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHessian2::coeff(int narg, char **arg)
{
  if (narg != 7) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_one = cut_global;
  double cut_one_min = 0;
  cut_one_min = force->numeric(FLERR,arg[2]);
  // Square minimum cutoff so that it'll be used in compute
  cut_one_min = cut_one_min*cut_one_min;
  cut_one = force->numeric(FLERR,arg[3]);

  // read i-i list and initialize potential parameters
  read_second_order(arg[4]);

  // read i-j list and initialize potential parameters
  read_second_order(arg[5]);

  // build the table of equilibrium position for the atoms
  build_poseq(arg[6]);

  /*double xlo = domain->boxlo[0];
  double xhi = domain->boxhi[0];
  double ylo = domain->boxlo[1];
  double yhi = domain->boxhi[1];
  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  fprintf(debug_read_hessian, "%f\n", xlo);
  fprintf(debug_read_hessian, "%f\n", xhi);
  fprintf(debug_read_hessian, "%f\n", ylo);
  fprintf(debug_read_hessian, "%f\n", yhi);
  fprintf(debug_read_hessian, "%f\n", zlo);
  fprintf(debug_read_hessian, "%f\n", zhi);*/

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      cutmin[i][j] = cut_one_min;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHessian2::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Tersoff requires atom IDs");
  //if (force->newton_pair == 0)
  //  error->all(FLERR,"Pair style Fc requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHessian2::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  if (offset_flag) {
    //double alpha_dr = -alpha[i][j] * (cut[i][j] - r0[i][j]);
    double dr = cut[i][j] - r0[i][j];
    double rk = k[i][j] * dr;
    //offset[i][j] = d0[i][j] * (exp(2.0*alpha_dr) - 2.0*exp(alpha_dr)); 
    offset[i][j] = rk*dr;
  } else offset[i][j] = 0.0;

  k[j][i] = k[i][j];
  r0[j][i] = r0[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairHessian2::read_second_order(char *file)
{
  int params_per_line = 12;
  int params_neigh = 15;
  char **words_neigh = new char*[params_neigh+1]; 
  char **words = new char*[params_per_line+1];
  int bas_neigh,bas_ref,num_neigh;
  int natoms = atom->natoms;

  // open file on proc 0
  FILE *fp;
  if (comm->me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open second-order force constant file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'

  int n,nwords;
  char line[MAXLINE],*ptr, *listtype_ptr, listtype[MAXLINE];
  int eof = 0;
  int i;
  int j;
  int k;
  int l;
  double m; //force constants
  int fc_count = 0;
  //char ii = "ii\n";
  //char ij = "ij\n";

  listtype_ptr = fgets(listtype,MAXLINE,fp);
  //if (strcmp(listtype, "ii\n") == 0) fprintf(debug_read_files, listtype);
  if (strcmp(listtype, "ii\n") == 0){
    //fprintf(debug_read_files, listtype);
    memory->destroy(iilist);
    memory->create(iilist,natoms*9*8,"pair:iilist");
    for (int it = 0; it < natoms*9; it++){
      iilist[it] = 0;
    }
  }
  if (strcmp(listtype, "ij\n") == 0){
    //fprintf(debug_read_files, listtype);
    memory->destroy(ijlist);
    memory->create(ijlist,158004*8,"pair:ijlist");
    for (int it = 0; it < 158004; it++){
      ijlist[it] = 0;
    }
  }

  /*for (int it = 0; it < natoms; it++){
    for (int jt = 0; jt < 3; jt++) {
      for (int kt = 0; kt < natoms; kt++) {
        for (int lt = 0; lt < 3; lt++) {
          fc2[it][jt][kt][lt] = 0;;
          //fprintf(debug_read_hessian, "%d %d %d %d %f\n", it,jt,kt,lt,fc2_temp2[it][jt][kt][lt]);
        }
      }
    }
  }*/
  //fclose(debug_read_hessian);
  // extract the number of atoms in the unit cell
  //fprintf(debug_read_hessian, "asdf\n");
  while (!eof) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    //fprintf(debugfc2,"%d\n",n); // "n" is number of characters in the line

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    //fprintf(debugfc2, line);
    if (nwords == 0) {
      continue;    
    } else {
      nwords = 0;
      words_neigh[nwords++] = strtok(line," \t\n\r\f");
      while ((words_neigh[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

      bas_neigh = atoi(words_neigh[0]);
      bas_ref = atoi(words_neigh[nwords-2]);  
      num_neigh = atoi(words_neigh[nwords-5]);

      i = atoi(words_neigh[0]);
      j = atoi(words_neigh[1]);
      k = atoi(words_neigh[2]);
      l = atoi(words_neigh[3]);
      m = atof(words_neigh[4]);

      /*fprintf(debugfc2, "%s\n", words_neigh[0]);
      fprintf(debugfc2, "%s\n", words_neigh[1]);
      fprintf(debugfc2, "%s\n", words_neigh[2]);
      fprintf(debugfc2, "%s\n", words_neigh[3]);
      fprintf(debugfc2, "%s\n", words_neigh[4]);
      fprintf(debugfc2, "---------------------\n");*/

      //fprintf(debugfc2, "%f\n", m);
    }

    //fprintf(debugfc2, "%d %d %d %d %f\n", i-1,j-1,k-1,l-1,val);

    // Assign value to array
    if (strcmp(listtype, "ii\n") == 0){
      iilist[fc_count] = m*13.605698066*(1/0.529177249)*(1/0.529177249);
      //fprintf(debug_read_files, "ii\n");
    }
    if (strcmp(listtype, "ij\n") == 0){
      ijlist[fc_count] = m*13.605698066*(1/0.529177249)*(1/0.529177249);
      //fprintf(debug_read_files, "jj\n");
    }

    //double val = fc2_temp2[i-1][j-1][k-1][l-1];
    //fprintf(debugfc2, "%d %d %d %d %f\n", i,j,k,l,val);
    //fprintf(debugfc2, "%d\n", fc_count);
    fc_count = fc_count + 1;
    
    /*fc2i = fc2i + 1;
    if (fc2j == 64){
      fc2j = fc2j + 1;
    }*/

    /*fprintf(debugfc2, "%s\n", words_neigh[0]);
    fprintf(debugfc2, "%s\n", words_neigh[1]);
    fprintf(debugfc2, "%s\n", words_neigh[2]);
    fprintf(debugfc2, "%s\n", words_neigh[3]);
    fprintf(debugfc2, "%s\n", words_neigh[4]);
    fprintf(debugfc2, "---------------------\n");*/
 
    /*if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    //fprintf(debugfc2, line);
    //fprintf(debugfc2, "%d\n", nwords);
    if (nwords == 0) {
      continue;
    } else {
      nwords = 0;
      //fprintf(debugfc2, "%d\n", nwords);
      words[nwords++] = strtok(line," \t\n\r\f");
      words[0] = strtok(line," ");
      words[1] = strtok(line," ");
      words[2] = strtok(NULL," ");
      fprintf(debugfc2, "%s\n", words[0]);
      fprintf(debugfc2, "%s\n", words[1]);
      fprintf(debugfc2, "%s\n", words[2]);
      fprintf(debugfc2, words[3]);
      words[nwords++] = strtok(line," \t\n\r\f");
      nufc2 = atoi(words[0]);
      eof = 1;
    }*/
  }
  //fprintf(debug_read_hessian, "asdf\n");

  // Read iilist for testing
  //int test = strcmp(listtype, "ii\n");
  //fprintf(debug_read_files, "%d\n", test);
  /*if (strcmp(listtype, "ii\n") == 0){
    for (int it = 0; it < natoms*9; it++){
      double val = iilist[it];
      fprintf(debug_read_files, "%f\n",val);
    }
  }*/

  // Read ijlist for testing
  /*if (strcmp(listtype, "ij\n") == 0){
    for (int it = 0; it < 16128; it++){
      double val = ijlist[it];
      fprintf(debug_read_files, "%f\n",val);
    }
  }*/

  //fprintf(debugfc2, "%d\n", n);
  //fprintf(debugfc2, "%d\n", nufc2);
  eof = 0;
  if (fp == NULL) fclose(fp);

  if (comm->me == 0) printf("second-order force constant file read!\n");
  //delete [] words;
}

/*----------------------------------------------------------------------
            build the equilibrium position table          
---------------------------------------------------------------------- */
  
void PairHessian2::build_poseq(char *file)
{
  char line[MAXLINE],*ptr;
  FILE *fp;
  int n,nwords,nuc_poscar;
  char **words = new char*[MAXLINE];
  int natoms = atom->natoms;
  int info = 0;

  memory->create(poseq,natoms,3,"pair:poseq");

  for (int i = 0; i < natoms; i++) {
    for (int j = 0; j < 3; j++) {
      poseq[i][j] = 0.0;
    }
  }

  fp = fopen(file, "r");

  if (fp == NULL) {
    char str[MAXLINE];
    sprintf(str,"Cannot open the equilibrium position file %s", file);
    error->all(FLERR,str);
  }

  int eof = 0;
  int count = 0;
  // extract the equilibrium position
  while (!eof) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) {
      continue;
    } else {
      nwords = 0;
      words[nwords++] = strtok(line," \t\n\r\f");
      while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;
      poseq[count][0] = atof(words[0]);
      poseq[count][1] = atof(words[1]);
      poseq[count][2] = atof(words[2]);
      count++;
    }
  }

  /*if (comm->me == 0) {
    for (int i = 0; i < natoms; i++) {
      fprintf(debug_poseq,"%.12f %.12f %.12f\n",poseq[i][0],poseq[i][1],poseq[i][2]);
    }
  }*/

  if (count != natoms) error->all(FLERR,"Number of atoms in the equilibrium position file different from the total number of atoms");

  if (comm->me == 0) printf("The equilibrium position table has been built!\n");
}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairHessian2::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&k[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairHessian2::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&k[i][j],sizeof(double),1,fp);
          fread(&r0[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&k[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairHessian2::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairHessian2::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairHessian2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i][i],r0[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairHessian2::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",
              i,j,k[i][j],r0[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairHessian2::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{

  double phi;

  double r = sqrt(rsq);
  double dr = r - r0[itype][jtype];
  double rk = k[itype][jtype] * dr;

  fforce = 0;
  if (r > 0.0) fforce = factor_lj * -2.0*rk/r;

  phi = rk*dr;
  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void *PairHessian2::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"k") == 0) return (void *) k;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  return NULL;
}
