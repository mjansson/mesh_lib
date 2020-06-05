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
#include <vector/types.h>

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

#if FOUNDATION_COMPILER_MSVC
#pragma warning(push)
#pragma warning(disable : 4201)
#endif

typedef struct mesh_config_t mesh_config_t;

typedef enum mesh_attribute_type mesh_attribute_type;

typedef union mesh_vector_t mesh_vector_t;

typedef mesh_vector_t mesh_coordinate_t;
typedef mesh_vector_t mesh_normal_t;
typedef mesh_vector_t mesh_tangent_t;
typedef mesh_vector_t mesh_bitangent_t;
typedef mesh_vector_t mesh_color_t;
typedef struct mesh_uv_t mesh_uv_t;
typedef struct mesh_vertex_t mesh_vertex_t;
typedef struct mesh_triangle_t mesh_triangle_t;
typedef struct mesh_vertex_triangle_t mesh_vertex_triangle_t;
typedef struct mesh_attribute_t mesh_attribute_t;
typedef struct mesh_partition_t mesh_partition_t;
typedef struct mesh_topology_t mesh_topology_t;
typedef struct mesh_t mesh_t;

#define MESH_INVALID_VERTEX ((unsigned int)-1)
#define MESH_BUCKET_SIZE 1000

enum mesh_attribute_type { MESH_ATTRIBUTE_FLOAT, MESH_ATTRIBUTE_FLOAT4 };

struct mesh_config_t {
	size_t _unused;
};

union mesh_vector_t {
	struct {
		float x, y, z, w;
	};
	float arr[4];
	vector_t v;
};

struct mesh_uv_t {
	real u;
	real v;
};

struct mesh_attribute_t {
	string_t name;
	mesh_attribute_type type;
	//! Attribute values (mesh_vector_t)
	bucketarray_t values;
	//! Index into attribute values, one per associated item, per-vertex or per-triangle (unsigned int)
	bucketarray_t index;
};

struct mesh_vertex_t {
	//! Coordinate index
	unsigned int coordinate;
	//! Normal index
	unsigned int normal;
	//! UV index
	unsigned int uv[2];
	//! Tangent index
	unsigned int tangent;
	//! Binormal index
	unsigned int binormal;
	//! Color index
	unsigned int color;
	//! Index of next vertex sharing the same coordinate index forming a linked ring list,
	//! or MESH_INVALID_VERTEX if none
	unsigned int adjacent;
};

struct mesh_triangle_t {
	//! Index of the three corner vertex
	unsigned int vertex[3];
	//! Index of the three adjacent triangles along the corresponding edge from vertex n -> n+1
	// An adjacent triangle shares both vertex coordinates and is winded the opposite direction along the edge,
	// such that the adjacent triangle edge goes from vertex n+1 -> n
	unsigned int adjacent[3];
	//! Index of original face in the case of higher order polygon triangulation
	unsigned int face;
	//! Material index
	unsigned int material;
};

struct mesh_face_t {
	//! Number of corners
	unsigned int vertex_count;
	//! Index of corner vertex
	unsigned int vertex[FOUNDATION_FLEXIBLE_ARRAY];
};

struct mesh_vertex_triangle_t {
	//! Valence
	unsigned int valence;
	//! Capacity of triangle index storage
	unsigned int capacity;
	//! Offset into triangle index storage
	size_t offset;
};

struct mesh_topology_t {
	//! Vertex to triangle map, with valence and offset in triangle index storage (mesh_vertex_triangle_t)
	bucketarray_t vertex_triangle_map;
	//! Vertex to triangle map storage (unsigned int)
	bucketarray_t vertex_triangle_store;
};

struct mesh_t {
	//! Coordinate bucket array (mesh_coordinate_t)
	bucketarray_t coordinate;
	//! Normal bucket array (mesh_normal_t)
	bucketarray_t normal;
	//! UV coordinate bucket array (mesh_uv_t)
	bucketarray_t uv[2];
	//! Tangent bucket array (mesh_tangent_t)
	bucketarray_t tangent;
	//! Bitangent bucket array (mesh_bitangent_t)
	bucketarray_t bitangent;
	//! Coordinate to vertex map (unsigned int)
	bucketarray_t coordinate_to_vertex;
	//! Attributes per vertex (attribute count equals vertex count)
	mesh_attribute_t* attribute_vertex;
	//! Attributes per triangle (attribute count equals triangle count)
	mesh_attribute_t* attribute_triangle;
	//! Vertex bucket array (mesh_vertex_t)
	bucketarray_t vertex;
	//! Triangle bucket array (mesh_triangle_t)
	bucketarray_t triangle;
	//! Bounding box minimum
	vector_t bounds_min;
	//! Bounding box maximum
	vector_t bounds_max;
	//! Partitioning of the mesh
	mesh_partition_t* partition;
	//! Topology of the mesh
	mesh_topology_t* topology;
};

#if FOUNDATION_COMPILER_MSVC
#pragma warning(pop)
#endif
