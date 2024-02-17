/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2022, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */



//=============================================================================
//
//  CLASS newClass
//
//=============================================================================


#ifndef OPENMESH_UTILS_GLCONSTASSTRING_HH
#define OPENMESH_UTILS_GLCONSTASSTRING_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/System/config.h>


//== FORWARDDECLARATIONS ======================================================


//== NAMESPACES ===============================================================

namespace OpenMesh {
namespace Utils {

//== CLASS DEFINITION =========================================================

inline
const char *GLenum_as_string( GLenum _m )
{
#define MODE(M) case M:return #M
  switch( _m )
  {
    MODE(GL_POINTS);
    MODE(GL_LINES);
    MODE(GL_LINE_STRIP);
    MODE(GL_LINE_LOOP);
    MODE(GL_TRIANGLES);
    MODE(GL_TRIANGLE_STRIP);
    MODE(GL_TRIANGLE_FAN);
    MODE(GL_QUADS);
    MODE(GL_QUAD_STRIP);
    MODE(GL_POLYGON);
    default: return "<unknown>";
  }
#undef MODE
}

//=============================================================================
} // namespace Utils
} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_UTILS_GLCONSTASSTRING_HH defined
//=============================================================================

