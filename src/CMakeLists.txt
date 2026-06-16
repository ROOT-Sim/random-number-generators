# SPDX-FileCopyrightText: 2008-2026 HPCS Group <rootsim@googlegroups.com>
# SPDX-License-Identifier: GPL-3.0-only

if(WIN32)
    set(EXTRA_LIBS ${EXTRA_LIBS} winmm)
endif()

add_library(rsrng STATIC random.c xxtea.c)
target_include_directories(rsrng PRIVATE . PUBLIC ./include)
target_link_libraries(rsrng ${EXTRA_LIBS})

install(DIRECTORY include/ DESTINATION include)
install(TARGETS rsrng LIBRARY DESTINATION lib)
