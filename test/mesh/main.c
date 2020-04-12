/* mesh.c  -  mesh library  -  Public Domain  -  2020 Mattias Jansson
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

#include <mesh/mesh.h>

#include <foundation/foundation.h>
#include <test/test.h>

static application_t
test_mesh_application(void) {
	application_t app;
	memset(&app, 0, sizeof(app));
	app.name = string_const(STRING_CONST("Mesh tests"));
	app.short_name = string_const(STRING_CONST("test_mesh"));
	app.company = string_const(STRING_CONST(""));
	app.flags = APPLICATION_UTILITY;
	app.exception_handler = test_exception_handler;
	return app;
}

static memory_system_t
test_mesh_memory_system(void) {
	return memory_system_malloc();
}

static foundation_config_t
test_mesh_foundation_config(void) {
	foundation_config_t config;
	memset(&config, 0, sizeof(config));
	return config;
}

static int
test_mesh_initialize(void) {
	mesh_config_t config;
	memset(&config, 0, sizeof(config));
	log_set_suppress(HASH_MESH, ERRORLEVEL_INFO);
	return mesh_module_initialize(config);
}

static void
test_mesh_finalize(void) {
	mesh_module_finalize();
}

static void
test_mesh_declare(void) {
}

static test_suite_t test_mesh_suite = {test_mesh_application,
                                       test_mesh_memory_system,
                                       test_mesh_foundation_config,
                                       test_mesh_declare,
                                       test_mesh_initialize,
                                       test_mesh_finalize,
                                       0};

#if BUILD_MONOLITHIC

int
test_mesh_run(void);

int
test_mesh_run(void) {
	test_suite = test_mesh_suite;
	return test_run_all();
}

#else

test_suite_t
test_suite_define(void);

test_suite_t
test_suite_define(void) {
	return test_mesh_suite;
}

#endif
