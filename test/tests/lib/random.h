/**
* @file random.h
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

extern double Random(void);
extern uint64_t RandomU64(void);
extern double Poisson(void);
#define Expent(mean) ((mean) * Poisson())
extern double Normal(void);
extern int RandomRange(int min, int max);
extern int RandomRangeNonUniform(int x, int min, int max);
extern double Gamma(unsigned ia);
extern unsigned Zipf(double skew, unsigned limit);
