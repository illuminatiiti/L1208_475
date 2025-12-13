#ifndef L1208_475_HPP
#define L1208_475_HPP

/* точка x,y */ 
typedef struct {
    int x, y;    
} point;

/* тип для набора точек на котором будем решать tsp (traveling salesman problem) */ 
typedef struct {
    int n;         /* сколько точек в проблеме */          
    point p[1002]; /* массив точек */   
} tsp_instance;

/* структура для представления решения */ 
typedef struct {
    double cost;   /* стоимость решения */
    int n;         /* количество перестановок */
    int p[1002];   /* массив индексов */
} tsp_solution;

void solution_count_update(tsp_solution* s, tsp_instance* t);
void anneal(tsp_instance* t, tsp_solution* s);
void repeated_annealing(tsp_instance* t, int nsamples, tsp_solution* bestsol);

#endif