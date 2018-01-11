/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   1) Manage the exchange of data between Code_Saturne and the pre-processor
 *   2) Define (conforming or non-conforming) mesh joinings.
 *   3) Define (conforming or non-conforming) periodicity.
 *   4) Define thin walls.
 *   5) Modify the geometry and mesh.
 *   6) Smoothe the mesh.
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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_selector.h"

#include "cs_base.h"
#include "cs_mesh_boundary.h"
#include "cs_join.h"
#include "cs_join_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"
#include "cs_mesh_extrude.h"
#include "cs_mesh_smoother.h"
#include "cs_mesh_warping.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_preprocessor_data.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.c
 *
 * \brief Definition and modification of the calculation mesh.
 *
 * See \subpage cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define periodic faces.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_periodicity(void)
{
  int    join_num;
  int    verbosity = 1;     /* per-task dump if > 1, debug level if >= 3 */
  int    visualization = 1; /* debug level if >= 3 */
  float  fraction = 0.10, plane = 25.;

  const double translationx[3] = {12.566, 0.0, 0.0}; /* X-Translation vector */
  const double translationz[3] = {0.0, 0.0, 6.284}; /* Z-Translation vector */

  join_num = cs_join_perio_add_translation("x<=0.2 or x>=12.366",
                                           fraction,
                                           plane,
                                           verbosity,
                                           visualization,
                                           translationx);
  join_num = cs_join_perio_add_translation("z<=0.1 or z>=6.084",
                                           fraction,
                                           plane,
                                           verbosity,
                                           visualization,
                                           translationz);
}

