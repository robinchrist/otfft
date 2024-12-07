from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.tools.files import copy
import os


class otfftppRecipe(ConanFile):
    name = "otfftpp"
    version = "0.0.1"
    package_type = "library"

    # Optional metadata
    license = "MIT"
    url = "https://github.com/conan-io/conan-center-index"
    description = "Very fast FFT - Heavily modified verison of OTFFT by OK Ojisan(Takuya OKAHISA)"
    topics = ("FFT", "SIMD")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "abi_affecting_cflags": ["ANY"], #comma separated list of c-flags, will be converted to a set for normalization and sorting
        "abi_affecting_cxxflags": ["ANY"], #comma separated list of c-flags, will be converted to a set for normalization and sorting
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "abi_affecting_cflags": "",
        "abi_affecting_cxxflags": "",
    }

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "src/*", "include/*", "external/*"

    def _abi_affecting_cflags(self, info=False):
        options = self.info.options if info else self.options
        return sorted(set(str(options.abi_affecting_cflags).split(",")))

    def _abi_affecting_cxxflags(self, info=False):
        options = self.info.options if info else self.options
        return sorted(set(str(options.abi_affecting_cxxflags).split(",")))

    def config_options(self):
        if self.settings.os == "Windows":
            self.options.rm_safe("fPIC")

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.variables["OTFFT_BUILD_SIMPLE_BENCHMARK"] = False
        tc.variables["OTFFT_BUILD_TESTS"] = False
        tc.extra_cflags.extend(self._abi_affecting_cflags())
        tc.extra_cxxflags.extend(self._abi_affecting_cxxflags())
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE.txt", src=self.source_folder, dst=os.path.join(self.package_folder, "licenses"))
        
        cmake = CMake(self)
        cmake.install()

    def package_id(self):
        # normalize the the extra flags (sorted+comma separated)
        self.info.options.abi_affecting_cflags = ",".join(self._abi_affecting_cflags(info=True))
        self.info.options.abi_affecting_cxxflags = ",".join(self._abi_affecting_cxxflags(info=True))

    def package_info(self):
        self.cpp_info.libs = ["otfftpp"]

    

    

