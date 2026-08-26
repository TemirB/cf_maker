# CMake generated Testfile for 
# Source directory: /Users/temir-balinov/diploma/programs/cf_maker
# Build directory: /Users/temir-balinov/diploma/programs/cf_maker/build-asan
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[core]=] "/Users/temir-balinov/diploma/programs/cf_maker/build-asan/test_core")
set_tests_properties([=[core]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/temir-balinov/diploma/programs/cf_maker/CMakeLists.txt;155;add_test;/Users/temir-balinov/diploma/programs/cf_maker/CMakeLists.txt;0;")
add_test([=[config]=] "/Users/temir-balinov/diploma/programs/cf_maker/build-asan/test_config")
set_tests_properties([=[config]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/temir-balinov/diploma/programs/cf_maker/CMakeLists.txt;177;add_test;/Users/temir-balinov/diploma/programs/cf_maker/CMakeLists.txt;0;")
add_test([=[fit]=] "/Users/temir-balinov/diploma/programs/cf_maker/build-asan/test_fit")
set_tests_properties([=[fit]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/temir-balinov/diploma/programs/cf_maker/CMakeLists.txt;202;add_test;/Users/temir-balinov/diploma/programs/cf_maker/CMakeLists.txt;0;")
subdirs("external/nlohmann_json")
