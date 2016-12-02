# This file is part of closest, a library for cull-based and 
# cell-based spatial search of k-nearest neighbors in n-dimensions.
# Copyright (C) 2011-2016 Daniel Pena <trifling.github@gmail.com>
#
# Closest is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Closest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with closest. If not, see <http://www.gnu.org/licenses/>.
 
set(CMAKE_BUILD_TYPE Release)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable( closest_test test.c )
target_link_libraries( closest_test LINK_PUBLIC closest m )

add_executable( closest_timing timing.c kdtree.c )
target_link_libraries( closest_timing LINK_PUBLIC closest m )

