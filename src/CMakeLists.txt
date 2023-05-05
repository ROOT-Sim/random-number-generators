# SPDX-FileCopyrightText: 2008-2023 HPDCS Group <rootsim@googlegroups.com>
# SPDX-License-Identifier: GPL-3.0-only

add_library(rsrng STATIC random.c xxtea.c)
target_include_directories(rsrng PRIVATE . ./include)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS rsrng LIBRARY DESTINATION lib)
