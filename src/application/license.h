/*******************************************************************************
Cubiquity - A micro-voxel engine for games and other interactive applications

Written by David Williams

To the extent possible under law, the author(s) have dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

You should have received a copy of the CC0 Public Domain Dedication along with
this software. If not, see http://creativecommons.org/publicdomain/zero/1.0/.
*******************************************************************************/

#ifndef CUBIQUITY_APP_LICENSE_H
#define CUBIQUITY_APP_LICENSE_H

#include "base/logging.h"

void printLicense() {
    std::string license = R"(
Cubiquity
---------
A micro-voxel engine for games and other interactive applications

Written by David Williams

To the extent possible under law, the author(s) have dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

You should have received a copy of the CC0 Public Domain Dedication along with
this software. If not, see http://creativecommons.org/publicdomain/zero/1.0/.


---------------------- Licensing of third-party libraries ----------------------

Cubiquity uses the following open source libraries:

* CLI11
    - Command line parser for C++11 and beyond that provides a rich feature set
      with a simple and intuitive interface.
    - 3-Clause BSD license
    - https://github.com/CLIUtils/CLI11

* {fmt}
    - A modern formatting library
    - MIT license
    - https://fmt.dev

* Glad
    - Multi-Language GL/GLES/EGL/GLX/WGL Loader-Generator based on the
      official specs
    - Only the output is included (not Glad itself) which is public domain
    - https://github.com/Dav1dde/glad
	
* linalg
    - Single header short vector math library for C++
    - Public domain
    - https://github.com/sgorsten/linalg

* PoolSTL
    - Light, self-contained, thread pool-based implementation of C++17 parallel
      standard library algorithms.
    - MIT license
    - https://github.com/alugowski/poolSTL

* SDL
    - Simple DirectMedia Layer
    - zlib license
    - https://libsdl.org/

* Simplex noise
    - Implementation of Simplex noise
    - Public domain
    - https://github.com/stegu/perlin-noise

* STB libraries
    - Single-file public domain (or MIT licensed) libraries for C/C++
    - Public domain
    - https://github.com/nothings/stb

* Tiny Obj Loader
    - Tiny but powerful single file wavefront obj loader written in C++03
    - MIT license
    - https://github.com/tinyobjloader/tinyobjloader

* toml11
    - A feature-rich TOML language library for C++11/14/17/20.
    - MIT license
    - https://github.com/ToruNiina/toml11

The licenses (where applicable) for these libraries are given below:

3-Clause BSD License
--------------------

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

MIT license
-----------

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the “Software”), to
    deal in the Software without restriction, including without limitation the
    rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
    sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
    IN THE SOFTWARE.


zlib license
------------

    This software is provided 'as-is', without any express or implied
    warranty.  In no event will the authors be held liable for any damages
    arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute it
    freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
       claim that you wrote the original software. If you use this software
       in a product, an acknowledgment in the product documentation would be
       appreciated but is not required.
    2. Altered source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.
    3. This notice may not be removed or altered from any source distribution.

)";
    print("{}", license);
}

#endif // CUBIQUITY_APP_LICENSE_H
