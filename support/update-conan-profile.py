import pathlib
import sys
import re

path = pathlib.Path(sys.argv[1])
build_type = sys.argv[2]
c_compiler = sys.argv[3]
cxx_compiler = sys.argv[4]
cxx_standard = sys.argv[5]

text = path.read_text()
text = re.sub("build_type=.*", f"build_type={build_type}", text)
text = re.sub("compiler.cppstd=.*", f"compiler.cppstd={cxx_standard}", text)
text = text + "\n[conf]"
text = text + f"\ntools.build:compiler_executables={{\"cpp\": \"{cxx_compiler}\", \"c\": \"{ c_compiler }\"}}"
text = text + "\ntools.cmake.cmaketoolchain:generator=Ninja"
path.unlink()
path.write_text(text)