module optimization_mod

   ! ----------------------------------------------------------------
   ! QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA
   !
   ! QUantum Mechanical VIbratioal Anlalysis or QUMVIA, for short, is
   ! an original implementation of Vibrational Self-consistent Field
   ! (VSCF) and Vibrational Configuration Interaction (VCI) using
   ! distributed gaussian basis set with analytical integrals and a
   ! configuration selection algorithm.
   !
   ! Implementation by Diego Javier Alonso de Armino
   ! Copyright(C) 2015 Diego J. Alonso de Armiño. All rights reserved.
   !
   ! This file is part of QUMVIA.
   !
   ! QUMVIA is free software: you can redistribute it and/or modify
   ! it under the terms of the GNU General Public License as published by
   ! the Free Software Foundation, either version 3 of the License, or
   ! (at your option) any later version.
   !
   ! QUMVIA is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU General Public License for more details.
   !
   ! A copy of the licence can be found in the root directory of
   ! QUMVIA.  If not, see <http://www.gnu.org/licenses/>.
   !
   ! QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA
   ! ----------------------------------------------------------------

private

public :: optimize

contains

   include "optimizer_main.f90"
   include "steepest_descent.f90"
   include "conjugate_gradient.f90"
   include "quasi_newton.f90"
   include "opt_initialize.f90"
   include "opt_check_convergence.f90"

end module optimization_mod
