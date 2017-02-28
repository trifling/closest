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

#CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

PROJECT( Closest C CXX )

SET(VERSION 2.0.0)

SET(CMAKE_MODULE_PATH ${Closest_SOURCE_DIR}/cmake)

SET(CMAKE_BUILD_TYPE Release)

SET(CMAKE_C_STANDARD_REQUIRED 99)
SET(CMAKE_C_STANDARD 99)
SET(CMAKE_C_FLAGS "-Wall -Wextra")

SET(CMAKE_CXX_STANDARD_REQUIRED 11)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "-Wall -Wextra")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(python)
ADD_SUBDIRECTORY(doc)

#ENABLE_TESTING()

