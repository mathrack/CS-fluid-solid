/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* Code_Saturne version 5.0.4-patch */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_base.h"
#include "cs_convection_diffusion.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gui_util.h"
#include "cs_grid.h"
#include "cs_internal_coupling.h"
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
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_selector.h"
#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
*
* User functions to compute the gradient and the dissipation of fields
*
*----------------------------------------------------------------------------*/

// Compute the gradient of u
static void
_my_grad_u(cs_real_33_t *grad)
{
  // Get the calculation option from the field
  const cs_field_t *fid = CS_F_(u);
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(fid, key_cal_opt_id, &var_cal_opt);
  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;
  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);
  // Compute the gradient
  cs_field_gradient_vector(fid,
                           false, // use_previous_t
                           gradient_type,
                           halo_type,
                           1, // inc
                           grad);
}
// This wrapper outputs the tensor grad(u)
static void
_my_ugrad(const void *input, cs_real_t *output)
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
// This wrapper outputs the tensor mu_t * grad(u)
static void
_my_mut_ugrad(const void *input, cs_real_t *output)
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
// This wrapper outputs the scalar mu_t * sum_{i,j} ( dUi / dXj )**2
static void
_my_udiss(const void *input, cs_real_t *output)
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

// Compute the gradient of T
static void
_my_grad_t(cs_real_3_t *grad, const int fi)
{
  // Get the calculation option from the field
  const cs_field_t *fid = cs_field_by_id(fi);
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(fid, key_cal_opt_id, &var_cal_opt);
  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;
  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);
  // Compute the gradient
  cs_field_gradient_scalar(fid,
                           false, // use_previous_t
                           gradient_type,
                           halo_type,
                           1, // inc
                           true, // recompute_cocg
                           grad);
}
// This wrapper outputs the vector grad(T)
static void
_my_tgrad(const void *input, cs_real_t *output)
{ // Input argument is field name
  const int fid = cs_field_id_by_name( *((const char **)input) );
  // Allocate memory
  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);
  // Compute the gradient
  _my_grad_t(grad, fid);
  // Put it in the output
  const int n_cells =  m->n_cells;
  for (int cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 3; i++) {
        output[3*cell_id+i] = grad[cell_id][i];
    }
  }
  // Free memory
  BFT_FREE(grad);
}
// This wrapper outputs the vector mu_t * grad(T) / Prt
static void
_my_prt_tgrad(const void *input, cs_real_t *output)
{ // Input argument is field name
  const int fid = cs_field_id_by_name( *((const char **)input) );
  // Allocate memory
  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);
  // Compute the gradient
  _my_grad_t(grad, fid);
  // Put it in the output
  cs_real_t prt = cs_field_get_key_double(cs_field_by_id(fid),
                                          cs_field_key_id("turbulent_schmidt"));
  const cs_real_t *mut = CS_F_(mu_t)->val;
  const int n_cells =  m->n_cells;
  for (int cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 3; i++) {
      output[3*cell_id+i] = mut[cell_id]*grad[cell_id][i]/prt;
    }
  }
  // Free memory
  BFT_FREE(grad);
}
// This wrapper outputs the scalar mu_t * sum_{j} ( dT / dXj )**2 / Pr_t
static void
_my_tdiss(const void *input, cs_real_t *output)
{ // Input argument is field name
  const int fid = cs_field_id_by_name( *((const char **)input) );
  // Allocate memory
  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells_ext = m->n_cells_with_ghosts;
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);
  // Compute the gradient
  _my_grad_t(grad, fid);
  // Put it in the output
  cs_real_t prt = cs_field_get_key_double(cs_field_by_id(fid),
                                          cs_field_key_id("turbulent_schmidt"));
  const cs_real_t *mut = CS_F_(mu_t)->val;
  const int n_cells =  m->n_cells;
  for (int cell_id = 0; cell_id < n_cells; cell_id++) {
    output[cell_id] = 0.0;
    for (int i = 0; i < 3; i++) {
      output[cell_id] += grad[cell_id][i]*grad[cell_id][i];
    }
    output[cell_id] *= mut[cell_id] / prt;
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
  const int my_nb_scalar = 11;
  const char *my_scal_name[] = {"diric","neuma",
                                "g05a05","g05a1","g05a2",
                                "g1a05","g1a1","g1a2",
                                "g2a05","g2a1","g2a2"};
  for (int i = 0; i<my_nb_scalar; i++) {
    cs_parameters_add_variable(my_scal_name[i], 1);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(void)
{
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  const cs_field_t *fu = CS_F_(u);
  cs_field_get_key_struct(fu, key_cal_opt_id, &var_cal_opt);
  var_cal_opt.nswrgr = 1;
  cs_field_set_key_struct(fu, key_cal_opt_id, &var_cal_opt);

  const int my_nb_scalar = 11;
  const char *my_scal_name[] = {"diric","neuma",
                                "g05a05","g05a1","g05a2",
                                "g1a05","g1a1","g1a2",
                                "g2a05","g2a1","g2a2"};
  for (int i = 0; i<my_nb_scalar; i++) {
    cs_field_t *f1 = cs_field_by_name(my_scal_name[i]);
    cs_field_get_key_struct(f1, key_cal_opt_id, &var_cal_opt);
    var_cal_opt.isstpc = 1;
    var_cal_opt.nswrgr = 1;
    var_cal_opt.iwgrec = 1;
    cs_field_set_key_struct(f1, key_cal_opt_id, &var_cal_opt);
  }
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
  int nt_start = 500000;
  const int my_nb_scalar = 11;
  static const char *my_scal_name[] = {"diric","neuma",
                                "g05a05","g05a1","g05a2",
                                "g1a05","g1a1","g1a2",
                                "g2a05","g2a1","g2a2"};
  int moment_c_id[] = {-1};
  int n_fields = 1;
  int moment_c_id2[] = {-1, -1};
  int n_fields2 = 2;
  // mu_t, <U>, Rij
  // <T>, <T'**2>, <u'T'>
  {
    // mu_t
    int moment_f_id[] = {CS_F_(mu_t)->id};
    cs_time_moment_define_by_field_ids("mu_t",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       nt_start, // nt_start
                                       -1,   // t_start
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
    // <U>
    moment_f_id[0] = CS_F_(u)->id;
    cs_time_moment_define_by_field_ids("u_mean",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       nt_start, // nt_start
                                       -1,   // t_start
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
    // Rij
    cs_time_moment_define_by_field_ids("rij",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_VARIANCE,
                                       nt_start, // nt_start
                                       -1,   // t_start
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
    const char *name_ext[] = {"_mean", "_variance", "_flux"};
    for (int iii = 0; iii<my_nb_scalar; iii++) {
      char str[80];
      // <T>
      moment_f_id[0] = cs_field_id_by_name(my_scal_name[iii]);
      snprintf(str, sizeof(str), "%s%s", my_scal_name[iii], name_ext[0]);
      cs_time_moment_define_by_field_ids(str,
                                         n_fields,
                                         moment_f_id,
                                         moment_c_id,
                                         CS_TIME_MOMENT_MEAN,
                                         nt_start, // nt_start
                                         -1,   // t_start
                                         CS_TIME_MOMENT_RESTART_AUTO,
                                         NULL);
      // <T'**2>
      snprintf(str, sizeof(str), "%s%s", my_scal_name[iii], name_ext[1]);
      cs_time_moment_define_by_field_ids(str,
                                         n_fields,
                                         moment_f_id,
                                         moment_c_id,
                                         CS_TIME_MOMENT_VARIANCE,
                                         nt_start, // nt_start
                                         -1,   // t_start
                                         CS_TIME_MOMENT_RESTART_AUTO,
                                         NULL);
      // <u'T'>
      int moment_f_id2[] = {moment_f_id[0], CS_F_(u)->id};
      snprintf(str, sizeof(str), "%s%s", my_scal_name[iii], name_ext[2]);
      cs_time_moment_define_by_field_ids(str,
                                         n_fields2,
                                         moment_f_id2,
                                         moment_c_id2,
                                         CS_TIME_MOMENT_MEAN,
                                         nt_start, // nt_start
                                         -1,   // t_start
                                         CS_TIME_MOMENT_RESTART_AUTO,
                                         NULL);
    }
  }
  // Dissipation-related moments
  {
    // Average of grad(U)
    cs_time_moment_define_by_func("u_grad",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  _my_ugrad,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    // Variance of grad(U)
    cs_time_moment_define_by_func("u_diss",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  _my_ugrad,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_VARIANCE,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    // Average of mu_t * grad(U)
    cs_time_moment_define_by_func("mu_t_u_grad",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  _my_mut_ugrad,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    // Average of mu_t * sum_{i,j} [ ( dUi / dXj )**2 ]
    cs_time_moment_define_by_func("mu_t_u_diss",
                                  CS_MESH_LOCATION_CELLS,
                                  1,
                                  _my_udiss,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  nt_start,     // nt_start
                                  -1,   // t_start
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
    const char *name_ext[] = {"_grad", "_diss", "_mu_t_t_grad", "_mu_t_t_diss"};
    for (int iii = 0; iii<my_nb_scalar; iii++) {
      char str[80];
      // Average of grad(T)
      snprintf(str, sizeof(str), "%s%s", my_scal_name[iii], name_ext[0]);
      cs_time_moment_define_by_func(str,
                                    CS_MESH_LOCATION_CELLS,
                                    3,
                                    _my_tgrad,
                                    &(my_scal_name[iii]),
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_MEAN,
                                    nt_start,     // nt_start
                                    -1,   // t_start
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
      // Variance of grad(T)
      snprintf(str, sizeof(str), "%s%s", my_scal_name[iii], name_ext[1]);
      cs_time_moment_define_by_func(str,
                                    CS_MESH_LOCATION_CELLS,
                                    3,
                                    _my_tgrad,
                                    &(my_scal_name[iii]),
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_VARIANCE,
                                    nt_start,     // nt_start
                                    -1,   // t_start
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
      // Average of mu_t * grad(T) / Pr_t
      snprintf(str, sizeof(str), "%s%s", my_scal_name[iii], name_ext[2]);
      cs_time_moment_define_by_func(str,
                                    CS_MESH_LOCATION_CELLS,
                                    3,
                                    _my_prt_tgrad,
                                    &(my_scal_name[iii]),
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_MEAN,
                                    nt_start,     // nt_start
                                    -1,   // t_start
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
      // Average of mu_t * sum_{j} [ ( dT / dXj )**2 ] / Pr_t
      snprintf(str, sizeof(str), "%s%s", my_scal_name[iii], name_ext[3]);
      cs_time_moment_define_by_func(str,
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    _my_tdiss,
                                    &(my_scal_name[iii]),
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_MEAN,
                                    nt_start,     // nt_start
                                    -1,   // t_start
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define internal coupling options.
 *
 * Options are usually defined using cs_internal_coupling_add_entity.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling(void)
{

  const int my_nb_scalar = 11;
  const char *my_scal_name[] = {"diric","neuma",
                                "g05a05","g05a1","g05a2",
                                "g1a05","g1a1","g1a2",
                                "g2a05","g2a1","g2a2"};
  for (int i = 2; i<my_nb_scalar; i++) {
    int fid = cs_field_id_by_name(my_scal_name[i]);
    cs_internal_coupling_add_entity(fid);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volumes as internal coupling zones.
 *
 * These zones will be separated from the rest of the domain using automatically
 * defined thin walls.
 *
 * \param[in, out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling_add_volumes(cs_mesh_t  *mesh)
{
   cs_internal_coupling_add_volume(mesh,"-1.<y<1.");
}
