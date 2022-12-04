cmake_minimum_required(VERSION 3.12)

add_library(kaos INTERFACE)

target_sources(kaos INTERFACE
  ${CMAKE_CURRENT_LIST_DIR}/chaos.c



)