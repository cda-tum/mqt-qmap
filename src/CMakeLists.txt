set(MQT_QMAP_TARGET_NAME mqt-qmap)

# exact mapper project library
if(Z3_FOUND)
  # add LogicBlocks library
  add_subdirectory(logicblocks)

  add_subdirectory(cliffordsynthesis)
endif()

add_subdirectory(sc)

add_subdirectory(hybridmap)

add_subdirectory(na)

if(BUILD_MQT_QMAP_BINDINGS)
  add_subdirectory(python)
endif()
