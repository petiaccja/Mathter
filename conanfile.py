from conans import ConanFile, CMake, tools


class MathterConan(ConanFile):
    name = "Mathter"
    version = "1.0.0"
    license = "The Unlicence"
    author = "Péter Kardos mathter_library@outlook.com"
    url = "https://github.com/petiaccja/Mathter"
    description = "Powerful 3D math and small-matrix linear algebra library for games and science."
    topics = ("game-dev", "linear-algebra", "vector-math", "matrix-library")
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}
    generators = "cmake"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def source(self):
        self.run("git clone https://github.com/petiaccja/Mathter.git")

    def build(self):
        pass

    def package(self):
        self.copy("*.hpp", dst="include", src="Mathter")

    def package_info(self):
        self.cpp_info.libs = []

