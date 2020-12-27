import os

from conans import ConanFile, CMake, tools


class MathterTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"

    def build(self):
        cmake = CMake(self)
        # Current dir is "test_package/build/<build_id>" and CMakeLists.txt is
        # in "test_package"
        cmake.configure()
        cmake.build()

    def imports(self):
        pass;

    def test(self):
        if not tools.cross_building(self):
            os.chdir("bin")
            self.run(".%squick_test" % os.sep)
