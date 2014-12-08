/*
Copyright (c) 2014, LEAP Pty Ltd
All rights reserved.

Redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following 
conditions are met:

1. Redistributions of source code must retain the above 
   copyright notice, this list of conditions and the following 
   disclaimer.

2. Redistributions in binary form must reproduce the above 
   copyright notice, this list of conditions and the following 
   disclaimer in the documentation and/or other materials provided 
   with the distribution.

3. Neither the name of the copyright holder nor the names of its 
   contributors may be used to endorse or promote products derived 
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Discription:
UDF showing how to use a scheme variable in the DEEFINE_PROFILE
macro.

This example a real values scheme variable named velocity-inlet is 
defined in scheme. This value is used to do define a profile.  
*/
#include "udf.h"


static real _vel = 0.0;


DEFINE_PROFILE(velocity_inlet, t, i)
{
  face_t f;

  begin_f_loop(f,t)
    F_PROFILE(f,t,i) = _vel * 1e-3;
  end_f_loop(f,t)
}


DEFINE_ADJUST(update_scheme_variable, d)
{
  #if !RP_NODE
  _vel = RP_Get_Real("velocity-inlet");
  #endif
  host_to_node_real_1(_vel);  
}


DEFINE_INIT(init_scheme_variable,d)
{
  #if !RP_NODE
  _vel = RP_Get_Real("velocity-inlet");
  #endif
  host_to_node_real_1(_vel);  
}
