cmake_minimum_required(VERSION 3.15)

# Assuming pybind11 is found from root or scikit-build-core
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_estate py_estate.cpp)

# Ensure we can find headers from other modules
target_include_directories(_estate PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/source_base
    ${PROJECT_SOURCE_DIR}/source_estate
)

target_link_libraries(_estate PRIVATE estate)

# Install the extension module into the source_estate package
# when using scikit-build-core, we should install to the package directory.
# The package root is `source_estate` (implicitly, or defined in pyproject).
# Actually, our package name is `estate` in pyproject.toml, but the python package directory is `source_estate`.
# Let's check pyproject.toml again.
install(TARGETS _estate DESTINATION op2c)
