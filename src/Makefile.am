lib_LTLIBRARIES = libgimli.la

otherincludedir = $(includedir)/$(pgkincludedir)/kdtree++/
otherinclude_HEADERS = \
			kdtree++/AUTHORS \
			kdtree++/COPYING \
			kdtree++/NEWS \
			kdtree++/README \
			kdtree++/VERSION \
			kdtree++/allocator.hpp \
			kdtree++/function.hpp \
			kdtree++/iterator.hpp \
			kdtree++/kdtree.hpp \
			kdtree++/node.hpp \
			kdtree++/region.hpp 

pgkincludedir = gimli
pkginclude_HEADERS = \
			baseentity.h \
			blockmatrix.h \
			cholmodWrapper.h \
			curvefitting.h \
			datacontainer.h \
			dc1dmodelling.h \
			elementmatrix.h \
			em1dmodelling.h \
			exitcodes.h \
			expressions.h \
			gimli.h \
			gravimetry.h \
			integration.h \
			interpolate.h \
			inversion.h \
			inversionBase.h \
			inversionFramework.h \
			ipcClient.h\
			kdtreeWrapper.h \
			ldlWrapper.h \
			line.h \
			linSolver.h \
			matrix.h \
			mesh.h \
			meshentities.h \
			meshgenerators.h \
			memwatch.h \
			modellingbase.h \
			multimatrix.h \
			node.h \
			numericbase.h \
			optionmap.h \
			plane.h \
			platform.h \
			polynomial.h \
			pos.h \
			quaternion.h \
			regionManager.h \
			shape.h \
			solver.h \
			solverWrapper.h \
			sparsematrix.h \
			spline.h \
			stopwatch.h \
			trans.h \
			triangleWrapper.h \
			ttdijkstramodelling.h \
			vector.h \
			vectortemplates.h


libgimli_la_SOURCES = \
			cholmodWrapper.cpp \
			curvefitting.cpp \
			datacontainer.cpp \
			dc1dmodelling.cpp \
			elementmatrix.cpp \
			em1dmodelling.cpp \
			gimli.cpp \
			gravimetry.cpp \
			integration.cpp \
			interpolate.cpp \
			inversion.cpp \
			ipcClient.cpp\
			kdtreeWrapper.cpp \
			ldlWrapper.cpp \
			line.cpp \
			linSolver.cpp \
			mesh.cpp \
			mesh_io.cpp \
			meshentities.cpp \
			meshgenerators.cpp \
			memwatch.cpp \
			modellingbase.cpp \
			multimatrix.cpp \
			node.cpp \
			numericbase.cpp \
			optionmap.cpp \
			plane.cpp \
			platform.cpp \
			pos.cpp \
			quaternion.cpp \
			regionManager.cpp \
			shape.cpp \
			solverWrapper.cpp \
			spline.cpp \
			stopwatch.cpp \
			triangleWrapper.cpp \
			ttdijkstramodelling.cpp 

libgimli_la_LDFLAGS= \
	-no-undefined \
	-version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE) \
	-release $(LT_RELEASE)
 
#libgimli_la_LIBADD = ../external/lib/libtriangle.lib
#libgimli_la_LIBADD = -Wl,-ltriangle

#libgimli_la_LIBADD = ../external/lib/libamd.a \
# 			../external/lib/libcamd.a \
# 			../external/lib/libccolamd.a \
# 			../external/lib/libcholmod.a \
# 			../external/lib/libcolamd.a \
# 			../external/lib/libldl.a \
# 			../external/lib/libtriangle.a