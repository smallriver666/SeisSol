//-*-c++-*-
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Gilbert Brietzke (gilbert.brietzke AT lrz.de, http://www.lrz.de)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#ifndef CHECKPOINT_SIONLIB_FAULT_H
#define CHECKPOINT_SIONLIB_FAULT_H

#include "CheckPoint.h"
#include "Checkpoint/Fault.h"
#include <iostream>

namespace seissol
{
  namespace checkpoint
  {
    namespace sionlib
    {
      class Fault : public CheckPoint, virtual public seissol::checkpoint::Fault
      {
      private:
	/** Struct describing the  header information in the file */
	struct Header {
	  unsigned long identifier;
	  int timestepFault;
	};
	
      public:
	Fault()
	  : CheckPoint(0x7A127)
	{}	
	bool init(double* mu, double* slipRate1, double* slipRate2, double* slip, double* slip1, 
		  double* slip2, double* state, double* strength,
		  unsigned int numSides, unsigned int numBndGP);
	
	void initLate()
	{
	  if (numSides() == 0)
	    return; 
	  CheckPoint::initLate();
	}
	/**
	 * @param[out] timestepFault Time step of the fault writer in the checkpoint
	 *  (if the fault writer was active)
	 */
	void load(int &timestepFault);
	void write(int timestepFault);	
	void writeinit();
	string which(){return string("sion.fault");};
	void close() {
	if (numSides() == 0)
	    return; 
	  CheckPoint::close();
	}
      };
    } 
  }  
}

#endif // CHECKPOINT_SIONLIB_FAULT_H
