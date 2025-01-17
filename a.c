#include <libkahypar.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    kahypar_context_t *ctx = kahypar_context_new();
    kahypar_configure_context_from_file(ctx, "config/evo.ini");

    kahypar_hypernode_id_t n = 4;
    kahypar_hyperedge_id_t m = 3;
    double eps = 0.03;
    kahypar_partition_id_t k = 2;
    size_t indices[] = {0, 2, 4, 6};
    kahypar_hyperedge_id_t edges[] = {0, 3, 3, 1, 1, 2};

    kahypar_hyperedge_weight_t obj;
    kahypar_partition_id_t part[] = {0, 0, 0, 0};

    kahypar_partition(n, m, eps, k, NULL, NULL, indices, edges, &obj, ctx, part);

    printf("-> %d\n", obj);
    printf("%d %d %d %d\n", part[0], part[1], part[2], part[3]);

    kahypar_context_free(ctx);
}
