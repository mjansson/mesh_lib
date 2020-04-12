/* types.h  -  mesh library  -  Public Domain  -  2020 Mattias Jansson
 *
 * This library provides a cross-platform mesh library in C11 providing
 * mesh handling functionality to applications base on our foundation library.
 *
 * The latest source code maintained by Mattias Jansson is always available at
 *
 * https://github.com/mjansson/mesh_lib
 *
 * This library is put in the public domain; you can redistribute it and/or modify it without any
 * restrictions.
 *
 */

#pragma once

/*! \file types.h
    Mesh data types */

#include <foundation/platform.h>
#include <foundation/types.h>

#include <mesh/build.h>

#if defined(MESH_COMPILE) && MESH_COMPILE
#ifdef __cplusplus
#define MESH_EXTERN extern "C"
#define MESH_API extern "C"
#else
#define MESH_EXTERN extern
#define MESH_API extern
#endif
#else
#ifdef __cplusplus
#define MESH_EXTERN extern "C"
#define MESH_API extern "C"
#else
#define MESH_EXTERN extern
#define MESH_API extern
#endif
#endif

typedef struct mesh_config_t mesh_config_t;

#define MESH_INVALID_VERTEX ((unsigned int)-1)

struct mesh_config_t {
	size_t _unused;
};

struct mesh_coordinate_t {
	real x;
	real y;
	real z;
};

struct mesh_normal_t {
	real nx;
	real ny;
	real nz;
};

struct mesh_uv_t {
	real u;
	real v;
};

struct mesh_color_t {
	real r;
	real g;
	real b;
	real a;
};

struct mesh_vertex_t {
	//! Coordinate index
	unsigned int coordinate;
	//! Normal index
	unsigned int normal;
	//! UV index for each uv layer
	unsigned int uv[8];
	//! Color index for each color layer
	unsigned int color[2];
	//! Index of next vertex sharing the same coordinate index forming a linked ring list,
	//! or MESH_INVALID_VERTEX if none
	unsigned int adjacent;
};

struct mesh_triangle_t {
	//! Index of the three corner vertex
	unsigned int vertex[3];
	//! Index of the three adjacent triangles along the corresponding edge from vertex n -> n+1
	unsigned int adjacent[3];
	//! Index of original face in the case of higher order polygon triangulation
	unsigned int face;
};

struct mesh_face_t {
	//! Number of corners
	unsigned int corner_count;
	//! Index of corner vertex
	unsigned int vertex[FOUNDATION_FLEXIBLE_ARRAY];
};
