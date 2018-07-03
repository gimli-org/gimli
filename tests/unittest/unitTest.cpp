#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif

#ifdef HAVE_UNITTEST
    #include <cppunit/ui/text/TestRunner.h>

    #include <vector.h>
    GIMLI::Placeholder x__;

    #include "testVector.h"
    //#include "testGIMLiMisc.h" // not working on Win32 despite USE_IPC
    #include "testDataContainer.h"
//
    #include "testMesh.h"
    #include "testShape.h"
    #include "testGeometry.h"
    #include "testFEM.h"
    #include "testExternals.h"

#endif // HAVE_UNITTEST

#include <cstdlib>
#include <iostream>

int main( int argc, char **argv){
#ifdef HAVE_UNITTEST
    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
    runner.addTest( registry.makeTest() );
    runner.run();
#else
	std::cout << "cppunit is not installed but if you can read this message, applications can be build." << std::endl<< std::endl;
#endif // HAVE_UNITTEST
    return EXIT_SUCCESS;
}
