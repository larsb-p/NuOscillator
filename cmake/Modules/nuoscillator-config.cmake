SET(NUOSCILLATOR_LIB_LIST "-libOscProbCalcer -libOscillator")

SET(ROOT_INCLUDE_DIRS_SEP ${ROOT_INCLUDE_DIRS})
string(REPLACE ";" " -I" ROOT_INCLUDE_DIRS_SEP "-I${ROOT_INCLUDE_DIRS_SEP}")
SET(ROOT_DEFINITIONS_SEP ${ROOT_DEFINITIONS})
string(REPLACE ";" " " ROOT_DEFINITIONS_SEP "${ROOT_DEFINITIONS_SEP}")
SET(ROOT_LIBRARIES_SEP ${ROOT_LIBRARIES})
string(REPLACE ";" " -l" ROOT_LIBRARIES_SEP "-l${ROOT_LIBRARIES_SEP}")

# Set the creation date
string(TIMESTAMP CREATION_DATE "%d-%m-%Y")

configure_file(${CMAKE_CURRENT_LIST_DIR}/../Templates/nuoscillator-config.in
  "${PROJECT_BINARY_DIR}/nuoscillator-config" @ONLY)
install(PROGRAMS
  "${PROJECT_BINARY_DIR}/nuoscillator-config" DESTINATION
  bin)
