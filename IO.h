#define TRACE(x)
#define DEBUG(x)

struct point * read_etsp(const char *tsp_file_name);
extern FILE *report;
extern FILE *best_report;
extern FILE *try_result;
double round_distance (long int i, long int j);
void write_params( void ); 