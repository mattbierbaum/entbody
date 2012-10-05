
// these create the global variable names for the arguments to be stored
#define ARGVNAME(type,name)     type##_##name 

// these create the converter string
#define ARGV_CONVERTER(type)    ARGVNAME(ARGV_CONVERTER, type)

// these are the converters
#define ARGV_CONVERTER_int      atoi
#define ARGV_CONVERTER_double   atof
#define ARGV_CONVERTER_float    (float)atof
