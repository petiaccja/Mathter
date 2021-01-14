from conans import ConanFile, CMake, tools


class MathterConan(ConanFile):
    name = "mathter"
    version = "1.0.0"
    license = "The Unlicence"
    author = "Péter Kardos mathter_library@outlook.com"
    url = "https://github.com/petiaccja/Mathter"
    description = "Powerful 3D math and small-matrix linear algebra library for games and science."
    topics = ("game-dev", "linear-algebra", "vector-math", "matrix-library")
    exports_sources = "Mathter/*"

    def package(self):
        self.copy("*.hpp", dst="include/Mathter", src="Mathter")