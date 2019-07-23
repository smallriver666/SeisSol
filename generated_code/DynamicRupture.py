#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016-2018, SeisSol Group
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
#

import numpy as np
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.input import parseJSONMatrixFile
from yateto.memory import CSCMemoryLayout
from multSim import OptionalDimTensor

def addKernels(generator, aderdg, matricesDir, dynamicRuptureMethod):
  if dynamicRuptureMethod == 'quadrature':
    numberOfPoints = (aderdg.order+1)**2
  elif dynamicRuptureMethod == 'cellaverage':
    numberOfPoints = int(4**math.ceil(math.log(aderdg.order*(aderdg.order+1)/2,4)))
  else:
    raise ValueError('Unknown dynamic rupture method.')

  clones = dict()

  # Load matrices
  db = parseJSONMatrixFile('{}/dr_{}_matrices_{}.json'.format(matricesDir, dynamicRuptureMethod, aderdg.order), clones, alignStride=aderdg.alignStride, transpose=aderdg.transpose)

  # Determine matrices
  TinvT = Tensor('TinvT', (aderdg.numberOfQuantities(), aderdg.numberOfQuantities()))
  fluxSolverShape = (aderdg.numberOfQuantities(), aderdg.numberOfExtendedQuantities())
  fluxSolver    = Tensor('fluxSolver', fluxSolverShape)

  gShape = (numberOfPoints, aderdg.numberOfQuantities())
  QInterpolated = OptionalDimTensor('QInterpolated', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)

  fluxScale = Scalar('fluxScale')
  generator.add('rotateFluxMatrix', fluxSolver['qp'] <= fluxScale * aderdg.starMatrix(0)['qk'] * aderdg.T['pk'])

  def interpolateQGenerator(i,h):
    return QInterpolated['kp'] <= db.V3mTo2n[i,h][aderdg.t('kl')] * aderdg.Q['lq'] * TinvT['qp']

  interpolateQPrefetch = lambda i,h: QInterpolated
  generator.addFamily('evaluateAndRotateQAtInterpolationPoints', simpleParameterSpace(4,4), interpolateQGenerator, interpolateQPrefetch)

  nodalFluxGenerator = lambda i,h: aderdg.extendedQTensor()['kp'] <= aderdg.extendedQTensor()['kp'] + db.V3mTo2nTWDivM[i,h][aderdg.t('kl')] * QInterpolated['lq'] * fluxSolver['qp']
  nodalFluxPrefetch = lambda i,h: aderdg.I
  generator.addFamily('nodalFlux', simpleParameterSpace(4,4), nodalFluxGenerator, nodalFluxPrefetch)

  # Energy output
  # Minus and plus refer to the original implementation of Christian Pelties,
  # where the normal points from the plus side to the minus side
  QInterpolatedPlus = OptionalDimTensor('QInterpolatedPlus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  QInterpolatedMinus = OptionalDimTensor('QInterpolatedMinus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  slipRateInterpolated = Tensor('slipRateInterpolated', (numberOfPoints,3))
  slipInterpolated = Tensor('slipInterpolated', (numberOfPoints,3))
  absoluteSlipInterpolated = Tensor('absoluteSlipInterpolated', (numberOfPoints,))
  tractionInterpolated = Tensor('tractionInterpolated', (numberOfPoints,3))
  frictionalEnergy = Tensor('frictionalEnergy', ())
  timeWeight = Scalar('timeWeight')
  spaceWeights = Tensor('spaceWeights', (numberOfPoints,))

  selectTractionSpp = np.zeros((aderdg.numberOfQuantities(), 3), dtype=bool)
  selectTractionSpp[0,0] = True
  selectTractionSpp[3,1] = True
  selectTractionSpp[5,2] = True
  tractionPlusMatrix = Tensor('tractionPlusMatrix', selectTractionSpp.shape, selectTractionSpp, CSCMemoryLayout)
  tractionMinusMatrix = Tensor('tractionMinusMatrix', selectTractionSpp.shape, selectTractionSpp, CSCMemoryLayout)

  computeSlipRateInterpolated = slipRateInterpolated['kp'] <= QInterpolatedMinus['kq'] * aderdg.selectVelocity['qp'] - QInterpolatedPlus['kq'] * aderdg.selectVelocity['qp']
  generator.add('computeSlipRateInterpolated', computeSlipRateInterpolated)

  computeTractionInterpolated = tractionInterpolated['kp'] <= QInterpolatedMinus['kq'] * tractionMinusMatrix['qp'] + QInterpolatedPlus['kq'] * tractionPlusMatrix['qp']
  generator.add('computeTractionInterpolated', computeTractionInterpolated)

  addSlipRateInterpolated = slipInterpolated['kp'] <= slipInterpolated['kp'] + timeWeight * slipRateInterpolated['kp']
  generator.add('addSlipRateInterpolated', addSlipRateInterpolated)

  squareSlipRateInterpolated = absoluteSlipInterpolated['k'] <= slipRateInterpolated['kp'] * slipRateInterpolated['kp']
  generator.add('squareSlipRateInterpolated', squareSlipRateInterpolated)

  computeFrictionalEnergy = frictionalEnergy[''] <= frictionalEnergy[''] + timeWeight * tractionInterpolated['kp'] * slipRateInterpolated['kp'] * spaceWeights['k']
  generator.add('computeFrictionalEnergy', computeFrictionalEnergy)
