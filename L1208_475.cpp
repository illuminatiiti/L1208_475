#include <iostream>
#include <cmath>
using namespace std;

int solution_count;        
int const PRINT_FREQUENCY = 10000;
int const COOLING_STEPS = 1000;
int const STEPS_PER_TEMP = 1000;
double const K = 0.01;

void solution_count_update(tsp_solution* s, tsp_instance* t) {
    solution_count = solution_count + 1;
    if ((solution_count == 1) || ((solution_count % PRINT_FREQUENCY) == 0)) {
        printf("%d %7.1f\n", solution_count, solution_cost(s, t));
    }
}

void anneal(tsp_instance* t, tsp_solution* s) {
    int x, y;                       /* пара точек для перестановки */
    int i, j;                       /* счетчики */
    bool accept_win, accept_loss;   /* условия принятия перехода */
    double temperature;             /* текущая температура системы */
    double current_value;           /* текущее значение */
    double start_value;             /* значение на момент начала цикла */
    double delta;                   /* значение после обмена */
    double exponent;                /* показатель степени для функции энерг. состояния */

    temperature = INITIAL_TEMPERATURE;

    initialize_solution(t->n, s);
    current_value = solution_cost(s, t);

    for (i = 1; i <= COOLING_STEPS; i++) {
        temperature *= COOLING_FRACTION;

        start_value = current_value;

        for (j = 1; j <= STEPS_PER_TEMP; j++) {
            /* pick indices of elements to swap */
            x = random_int(1, t->n);
            y = random_int(1, t->n);

            delta = transition(s, t, x, y);
            accept_win = (delta < 0);       /* did swap reduce cost? */

            exponent = (-delta / current_value) / (K * temperature);
            accept_loss = (exp(exponent) > random_float(0, 1));

            if (accept_win || accept_loss) {
                current_value += delta;
            }
            else {
                transition(s, t, x, y);     /* reverse transition */
            }
            solution_count_update(s, t);
        }

        if (current_value < start_value) {  /* rerun at this temp */
            temperature /= COOLING_FRACTION;
        }
    }
}

void repeated_annealing(tsp_instance* t, int nsamples, tsp_solution* bestsol) {
    tsp_solution s;                 /* current tsp solution */
    double best_cost;               /* best cost so far */
    double cost_now;                /* current cost */
    int i;                          /* counter */

    initialize_solution(t->n, &s);
    best_cost = solution_cost(&s, t);
    copy_solution(&s, bestsol);

    for (i = 1; i <= nsamples; i++) {
        anneal(t, &s);
        cost_now = solution_cost(&s, t);
        if (cost_now < best_cost) {
            best_cost = cost_now;
            copy_solution(&s, bestsol);
        }
    }
}
