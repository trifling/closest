#  This file is part of closest, a library for k-nearest neighbors (kNN) search
#  in n-dimensions.  
#  Copyright (C) 2011-2016 Daniel Pena <trifling.github@gmail.com>
#
#  Closest is free software: you can redistribute it and/or modify it under the
#  terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later
#  version.
#
#  Closest is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#  details.
#
#  You should have received a copy of the GNU General Public License along with
#  closest. If not, see <http://www.gnu.org/licenses/>.

# use python setup.py to handle the python module when installing

FIND_PACKAGE(PythonInterp)

IF( PYTHONINTERP_FOUND )
   INSTALL(CODE "MESSAGE( \"-- Installing: python bindings\" )")
   INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/setup.py 
                --quiet build --build-base ${CMAKE_CURRENT_BINARY_DIR} install 
                --prefix ${CMAKE_INSTALL_PREFIX} 
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})")
ENDIF()

