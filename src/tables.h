#include <stddef.h> // NULL

#ifndef STATS_TABLES_H
#define STATS_TABLES_H
 
/*
CONSIDER: implement it as a hash table global variable 
once I have the hash table working on arrays, I can hash on (n,k), which will 
be significantly easier, but will require a significant modification

CONSIDER: moving Pascal's triangle to a separate header
*/

//extern size_t * pascals_triangle;
//extern size_t pascals_triangle_size;
extern size_t (*nCr) (size_t, size_t); // alias for pascals_triangle_get

#define PASCALS_TRIANGLE_INIT_SIZE 36

void tables_init(); // module initialization
void table_del(); // cleanup

void pascals_triangle_init(size_t size);
size_t pascals_triangle_get(size_t n, size_t k);
void pascals_triangle_del();

#endif // STATS_TABLES_H
