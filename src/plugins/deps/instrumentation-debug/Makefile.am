#####################################################################################
#      Copyright 2012 Barcelona Supercomputing Center                               #
#                                                                                   #
#      This file is part of the NANOS++ library.                                    #
#                                                                                   #
#      NANOS++ is free software: you can redistribute it and/or modify              #
#      it under the terms of the GNU Lesser General Public License as published by  #
#      the Free Software Foundation, either version 3 of the License, or            #
#      (at your option) any later version.                                          #
#                                                                                   #
#      NANOS++ is distributed in the hope that it will be useful,                   #
#      but WITHOUT ANY WARRANTY; without even the implied warranty of               #
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
#      GNU Lesser General Public License for more details.                          #
#                                                                                   #
#      You should have received a copy of the GNU Lesser General Public License     #
#      along with NANOS++.  If not, see <http://www.gnu.org/licenses/>.             #
#####################################################################################

include $(top_srcdir)/src/plugins/deps/version_commons.am
include $(top_srcdir)/src/common.am

if INSTRUMENTATION_DEBUG_VERSION
instrdebug_LTLIBRARIES = \
        libnanox-deps-plain.la\
        libnanox-deps-perfect-regions.la\
        libnanox-deps-regions.la\
	$(END)

nanos_instrdebug_cxxflags=$(common_includes) $(common_instrdebug_cxxflags)

libnanox_deps_plain_la_CXXFLAGS = $(nanos_instrdebug_cxxflags)
libnanox_deps_plain_la_LDFLAGS = $(ld_plugin_flags)
libnanox_deps_plain_la_SOURCES = $(plain_sources)

libnanox_deps_perfect_regions_la_CXXFLAGS = $(nanos_instrdebug_cxxflags)
libnanox_deps_perfect_regions_la_LDFLAGS = $(ld_plugin_flags)
libnanox_deps_perfect_regions_la_SOURCES = $(perfect_regions_sources)

libnanox_deps_regions_la_CXXFLAGS = $(nanos_instrdebug_cxxflags)
libnanox_deps_regions_la_LDFLAGS = $(ld_plugin_flags)
libnanox_deps_regions_la_SOURCES = $(regions_sources)
endif

