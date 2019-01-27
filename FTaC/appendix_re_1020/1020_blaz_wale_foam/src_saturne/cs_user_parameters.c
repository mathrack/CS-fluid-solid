/*============================================================================
 * User subroutines for input of calculation parameters.
 *============================================================================*/

/* Code_Saturne version 4.0.4 */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_gui_util.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
* User function to compute the gradient and dissipation of a field
* 
* This function is used to compute the moments
*
*----------------------------------------------------------------------------*/

static void
_my_grad_u(cs_real_33_t *grad)
{
  // Get the calculation option from the field
  const cs_field_t *fid = CS_F_(u);
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(fid, key_cal_opt_id, &var_cal_opt);
  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;
  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  // Compute the gradient
  cs_field_gradient_vector(fid,
                           false,
                           gradient_type,
                           halo_type,
                           1, // inc
                           grad);
}

static void
_my_grad_and_diss(const void *input,
                  cs_real_t *output)
{
  // Allocate memory
  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  cs_real_33_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_33_t);

  // Compute the gradient
  _my_grad_u(grad);

  // Put it in the output
  const int n_cells =  m->n_cells;
  for (int cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        output[9*cell_id+3*i+j] = grad[cell_id][i][j];
      }
    }
  }

  // Free memory
  BFT_FREE(grad);
}

static void
_my_grad_and_diss2(const void *input,
                  cs_real_t *output)
{
  // Allocate memory
  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  cs_real_33_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_33_t);

  // Compute the gradient
  _my_grad_u(grad);

  // Put it in the output
  const cs_real_t *mut = CS_F_(mu_t)->val;
  const int n_cells =  m->n_cells;
  for (int cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        output[9*cell_id+3*i+j] = mut[cell_id]*grad[cell_id][i][j];
      }
    }
  }

  // Free memory
  BFT_FREE(grad);
}
static void
_my_grad_and_diss3(const void *input,
                  cs_real_t *output)
{
  // Allocate memory
  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  cs_real_33_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_33_t);

  // Compute the gradient
  _my_grad_u(grad);

  // Put it in the output
  const cs_real_t *mut = CS_F_(mu_t)->val;
  const int n_cells =  m->n_cells;
  for (int cell_id = 0; cell_id < n_cells; cell_id++) {
    output[cell_id] = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        output[cell_id] += grad[cell_id][i][j]*grad[cell_id][i][j];
      }
    }
    output[cell_id] *= mut[cell_id];
  }

  // Free memory
  BFT_FREE(grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /* If the GUI is used, user fields should preferably be defined with the GUI,
     so that associated numerical options, boundary conditions, initializations,
     and such may also be defined using the GUI. */

  if (cs_gui_file_is_loaded())
    return;

  /*--------------------------------------------------------------------------*/

  /* Example: add 2 scalar variables ("species" in the GUI nomenclature).
   *
   * Note that at this (very early) stage of the data setup, fields are
   * not defined yet. Associated fields will be defined later (after
   * model-defined fields) in the same order as that used here, and
   * after user-defined variables defined throught the GUI, if used.
   *
   * Currently, only 1d (scalar) fields are handled.
   *
   * parameters for cs_parameters_add_variable():
   *   name             <-- name of variable and associated field
   *   dim              <-- variable dimension
   */

  if (false) {
    cs_parameters_add_variable("species_1", 1);
    cs_parameters_add_variable("tracer", 1);
  }

  /*--------------------------------------------------------------------------*/

  /* Example: add the variance of a user variable.
   *
   * parameters for cs_parameters_add_variable_variance():
   *   name          <-- name of variance and associated field
   *   variable_name <-- name of associated variable
   */

  if (false) {
    cs_parameters_add_variable_variance("variance_1",
                                        "species_1");
  }

  /*--------------------------------------------------------------------------*/

  /* Example: add a user property defined on boundary faces.
   *
   * parameters for cs_parameters_add_property():
   *   name        <-- name of property and associated field
   *   dim         <-- property dimension
   *   location_id <-- id of associated mesh location, which must be one of:
   *                     CS_MESH_LOCATION_CELLS
   *                     CS_MESH_LOCATION_INTERIOR_FACES
   *                     CS_MESH_LOCATION_BOUNDARY_FACES
   *                     CS_MESH_LOCATION_VERTICES
   */

  if (false) {
    cs_parameters_add_property("user_b_property_1",
                               1,
                               CS_MESH_LOCATION_BOUNDARY_FACES);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define linear solver options.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined.
 *
 * Available native iterative linear solvers include conjugate gradient,
 * Jacobi, BiCGStab, BiCGStab2, and GMRES. For symmetric linear systems,
 * an algebraic multigrid solver is available (and recommended).
 *
 * External solvers may also be setup using this function, the cs_sles_t
 * mechanism alowing such through user-define functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define time moments.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined, and before fine control of field output options
 * is defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_time_moments(void)
{
  /*
   * We compute temporal means of the type <f1*f2*f3*...*fn>
   * The fi's are variables defined by fields of a same location
   * (usually cells or boundary faces)

   * The parameters for time_moment_define_by_field_ids are:
   *   name         <--  name of associated moment
   *   n_fields     <--  number of associated fields
   *   field_id     <--  ids of associated fields
   *   component_id <--  ids of matching field components (-1 for all)
   *   type         <--  moment type (CS_TIME_MOMENT_MEAN
   *                     or CS_TIME_MOMENT_VARIANCE)
   *   nt_start     <--  starting time step (or -1 to use t_start)
   *   t_start      <--  starting time
   *   restart_mode <--  behavior in case or restart:
   *                     CS_TIME_MOMENT_RESTART_RESET,
   *                     CS_TIME_MOMENT_RESTART_AUTO, or
   *                     CS_TIME_MOMENT_RESTART_EXACT
   *   restart_name <--  name in previous run, NULL for default
   */

  int nt_start = 50000+1;

  {
    /* 1st order moment <mu_t> and <U>, calculated starting from time step */

    int moment_f_id[] = {CS_F_(mu_t)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("mu_t",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       nt_start, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
    moment_f_id[0] = CS_F_(u)->id;
    cs_time_moment_define_by_field_ids("U_mean",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       nt_start, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
    /* End of 1st order moments */
    //
    /* 2nd order moments rij */
    cs_time_moment_define_by_field_ids("rij",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_VARIANCE,
                                       nt_start, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
    /* End of 2nd order moments rij */
  }

  {
    /* Moments for the dissipation  */
    cs_time_moment_define_by_func("mu_t_u_diss",
                                  CS_MESH_LOCATION_CELLS,
                                  1,
                                  _my_grad_and_diss3,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    cs_time_moment_define_by_func("u_grad",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  _my_grad_and_diss,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    cs_time_moment_define_by_func("u_diss",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  _my_grad_and_diss,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_VARIANCE,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    cs_time_moment_define_by_func("mu_t_u_grad",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  _my_grad_and_diss2,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    /* End of dissipation-related moments */
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
