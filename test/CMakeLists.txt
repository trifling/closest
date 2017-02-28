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
 
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

ADD_EXECUTABLE( closest_test test.c )
TARGET_LINK_LIBRARIES( closest_test LINK_PUBLIC closest m )

ADD_EXECUTABLE( closest_test_array test_array.c )
TARGET_LINK_LIBRARIES( closest_test_array LINK_PUBLIC closest m )

ADD_EXECUTABLE( closest_test_metric test_metric.c )
TARGET_LINK_LIBRARIES( closest_test_metric LINK_PUBLIC closest m )

ADD_EXECUTABLE( closest_test_ln test_ln.c )
TARGET_LINK_LIBRARIES( closest_test_ln LINK_PUBLIC closest m )

ADD_EXECUTABLE( closest_test_cxx test_cxx.cpp )
TARGET_LINK_LIBRARIES( closest_test_cxx LINK_PUBLIC closest_cxx closest m )

ADD_EXECUTABLE( closest_timing timing.c kdtree.c )
TARGET_LINK_LIBRARIES( closest_timing LINK_PUBLIC closest m )


