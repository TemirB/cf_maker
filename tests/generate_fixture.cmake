# This directory contains only this test's generated fixtures.
file(MAKE_DIRECTORY "${OUTPUT_DIR}")
file(REMOVE "${OUTPUT_DIR}/input.root" "${OUTPUT_DIR}/config.json")
execute_process(COMMAND "${GENERATOR}" "${OUTPUT_DIR}" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Gaussian fixture generation failed: ${result}")
endif()
