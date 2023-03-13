/**
* @file ROOT-Sim.h
*
* @brief ROOT-Sim header for model development
*
* This header defines all the symbols which are needed to develop a model
* to be simulated on top of ROOT-Sim.
*
* This header is the only file which should be included when developing
* a simulation model. All function prototypes exposed to the application
* developer are exposed and defined here.
*
* SPDX-FileCopyrightText: 2008-2022 HPDCS Group <rootsim@googlegroups.com>
* SPDX-License-Identifier: GPL-3.0-only
*/
#pragma once

struct topology;

enum topology_geometry {
	TOPOLOGY_HEXAGON = 1,	//!< hexagonal grid topology
	TOPOLOGY_SQUARE,	//!< square grid topology
	TOPOLOGY_TORUS,		//!< a torus shaped grid topology (a wrapping around square topology)
	TOPOLOGY_RING,		//!< a ring shaped topology walkable in a single direction
	TOPOLOGY_BIDRING,	//!< a ring shaped topology
	TOPOLOGY_STAR,		//!< a star shaped topology
	TOPOLOGY_FCMESH,	//!< a fully-connected mesh
	TOPOLOGY_GRAPH,		//!< a (weighted) graph topology
};

enum topology_direction {
	DIRECTION_E,	//!< East direction
	DIRECTION_W,	//!< West direction
	DIRECTION_N,	//!< North direction
	DIRECTION_S,	//!< South direction
	DIRECTION_NE,	//!< North-east direction
	DIRECTION_SW,	//!< South-west direction
	DIRECTION_NW,	//!< North-west direction
	DIRECTION_SE,	//!< South-east direction

	DIRECTION_RANDOM, //!< Get a random direction, depending on the topology
};

/// An invalid direction, used as error value for the functions which return a LP id
#define INVALID_DIRECTION UINT64_MAX

extern lp_id_t CountRegions(struct topology *topology);
extern lp_id_t CountDirections(lp_id_t from, struct topology *topology);
extern lp_id_t GetReceiver(lp_id_t from, struct topology *topology, enum topology_direction direction);

extern void ReleaseTopology(struct topology *topology);
extern bool AddTopologyLink(struct topology *topology, lp_id_t from, lp_id_t to, double probability);
extern bool IsNeighbor(lp_id_t from, lp_id_t to, struct topology *topology);


// The following trick belongs to Laurent Deniau at CERN.
// https://groups.google.com/g/comp.std.c/c/d-6Mj5Lko_s?pli=1

/// Compute the number of arguments to a variadic macro
#define PP_NARG(...) \
         PP_NARG_(__VA_ARGS__,PP_RSEQ_N())
/// Compute the number of arguments to a variadic macro
#define PP_NARG_(...) \
         PP_128TH_ARG(__VA_ARGS__)
/// Enumerate the arguments' count to a variadic macro
#define PP_128TH_ARG( \
          _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
         _11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
         _21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
         _31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
         _41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
         _51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
         _61,_62,_63,_64,_65,_66,_67,_68,_69,_70, \
         _71,_72,_73,_74,_75,_76,_77,_78,_79,_80, \
         _81,_82,_83,_84,_85,_86,_87,_88,_89,_90, \
         _91,_92,_93,_94,_95,_96,_97,_98,_99,_100, \
         _101,_102,_103,_104,_105,_106,_107,_108,_109,_110, \
         _111,_112,_113,_114,_115,_116,_117,_118,_119,_120, \
         _121,_122,_123,_124,_125,_126,_127,N,...) N
/// Enumerate the arguments' count to a variadic macro in reverse order
#define PP_RSEQ_N() \
         127,126,125,124,123,122,121,120, \
         119,118,117,116,115,114,113,112,111,110, \
         109,108,107,106,105,104,103,102,101,100, \
         99,98,97,96,95,94,93,92,91,90, \
         89,88,87,86,85,84,83,82,81,80, \
         79,78,77,76,75,74,73,72,71,70, \
         69,68,67,66,65,64,63,62,61,60, \
         59,58,57,56,55,54,53,52,51,50, \
         49,48,47,46,45,44,43,42,41,40, \
         39,38,37,36,35,34,33,32,31,30, \
         29,28,27,26,25,24,23,22,21,20, \
         19,18,17,16,15,14,13,12,11,10, \
         9,8,7,6,5,4,3,2,1,0

extern struct topology *vInitializeTopology(enum topology_geometry geometry, int argc, ...);


/**
 * @brief Initialize a topology struct
 *
 * @param geometry The topology struct to initialize
 * @param ... The number of regions and/or the number of links, depending on the topology geometry
 * @return A pointer to the topology structure
 */
#define InitializeTopology(geometry, ...) vInitializeTopology(geometry, PP_NARG(__VA_ARGS__), __VA_ARGS__)
