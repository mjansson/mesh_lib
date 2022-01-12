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

typedef vector_t mesh_vector_t;
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
typedef struct mesh_partition_element_t mesh_partition_element_t;
typedef struct mesh_partition_node_t mesh_partition_node_t;
typedef struct mesh_partition_leaf_t mesh_partition_leaf_t;
typedef struct mesh_partition_lookup_t mesh_partition_lookup_t;
typedef struct mesh_topology_t mesh_topology_t;
typedef struct mesh_t mesh_t;

#define MESH_INVALID_INDEX ((unsigned int)-1)

#define MESH_PARTITION_NODE_MAX_ITEMS 16

typedef enum mesh_partition_flag { MESH_PARTITION_FLAG_LEAF = 1 } mesh_partition_flag;

typedef enum mesh_attribute_type { MESH_ATTRIBUTE_FLOAT, MESH_ATTRIBUTE_FLOAT4 } mesh_attribute_type;

struct mesh_config_t {
	size_t unused;
};

struct mesh_uv_t {
	real u;
	real v;
};

FOUNDATION_STATIC_ASSERT(sizeof(mesh_coordinate_t) == 16, "Invalid coordinate size");
FOUNDATION_STATIC_ASSERT(sizeof(mesh_uv_t) == 8, "Invalid uv size");

struct mesh_attribute_t {
	//! Attribute name
	string_t name;
	//! Attribute type
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
	//! Bitangent index
	unsigned int bitangent;
	//! Color index
	unsigned int color;
	//! Index of next vertex sharing the same coordinate index forming a linked ring list,
	//! or MESH_INVALID_VERTEX if none
	unsigned int adjacent;
};

struct mesh_triangle_t {
	//! Index of the three corner vertex
	unsigned int vertex[3];
	//! Index of the three adjacent triangles along the corresponding edge from corner n -> n+1
	// An adjacent triangle shares both vertex coordinates and is winded the opposite direction along the edge,
	// such that the adjacent triangle edge goes from this triangle corner n+1 -> n. Topological
	// adjacency is based on coordinates, other vertex attributes such as normal or uv might differ
	// between corners of the adjacent triangles ("hard" edge for normals, seam for textures etc).
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

struct mesh_partition_node_t {
	//! Split point
	vector_t split;
	//! Child nodes
	unsigned int child[8];
};

struct mesh_partition_leaf_t {
	//! Min point
	vector_t min;
	//! Max point
	vector_t max;
	//! Coordinate index
	unsigned int coordinate[MESH_PARTITION_NODE_MAX_ITEMS];
};

struct mesh_partition_element_t {
	unsigned int flags;
	union {
		mesh_partition_node_t node;
		mesh_partition_leaf_t leaf;
	} data;
};

struct mesh_partition_t {
	//! Node/leaf storage (mesh_partition_element_t)
	bucketarray_t element;
	//! Root node
	mesh_partition_element_t root;
	//! Coordinate to vertex map (unsigned int)
	bucketarray_t coordinate_to_vertex;
};

struct mesh_partition_lookup_t {
	//! Coordinate index
	unsigned int coordinate;
};

struct mesh_t {
	//! ID
	uuid_t id;
	//! Name
	string_t name;
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
	//! Color bucket array (mesh_color_t)
	bucketarray_t color;
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
	//! Next LOD
	mesh_t* next_lod;
};

#if FOUNDATION_COMPILER_MSVC
#pragma warning(pop)
#endif
