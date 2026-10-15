/* SPDX-License-Identifier: BSD-3-Clause */
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kiss_fftndr.h"

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "Allocation check failed at line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

static int check_shape(const int *dims, int ndims)
{
    int inverse, i;
    for (inverse = 0; inverse < 2; ++inverse) {
        size_t needed = 0, available;
        unsigned char *buffer;
        kiss_fftndr_cfg plan;
        CHECK(kiss_fftndr_alloc(dims, ndims, inverse, NULL, &needed) == NULL);
        CHECK(needed != 0 && needed != (size_t)-1);
        buffer = (unsigned char *)malloc(needed + 32);
        CHECK(buffer != NULL);
        memset(buffer, 0xa5, needed + 32);
        available = needed - 1;
        CHECK(kiss_fftndr_alloc(dims, ndims, inverse, buffer, &available) == NULL);
        CHECK(available == needed && buffer[0] == 0xa5);
        available = needed;
        CHECK(kiss_fftndr_alloc(dims, ndims, inverse, buffer, &available)
              == (kiss_fftndr_cfg)buffer);
        CHECK(available == needed);
        for (i = 0; i < 32; ++i)
            CHECK(buffer[needed + i] == 0xa5);
        free(buffer);
        plan = kiss_fftndr_alloc(dims, ndims, inverse, NULL, NULL);
        CHECK(plan != NULL);
        kiss_fftndr_free(plan);
    }
    return 0;
}

int main(void)
{
    const int valid[][3] = {{54, 54, 0}, {64, 64, 0}, {128, 128, 0},
                            {3, 5, 64}, {1, 2, 0}, {52, 52, 0}};
    const int invalid[][2] = {{0, 64}, {-1, 64}, {64, 0}, {64, 63},
                              {INT_MAX, 2}, {46341, 46342}};
    size_t i;
    for (i = 0; i < sizeof(valid) / sizeof(valid[0]); ++i)
        CHECK(check_shape(valid[i], valid[i][2] == 0 ? 2 : 3) == 0);
    for (i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
        size_t needed = 0;
        CHECK(kiss_fftndr_alloc(invalid[i], 2, 0, NULL, &needed) == NULL);
        CHECK(needed == (size_t)-1);
    }
    {
        /* The index count fits int, but the plan cannot fit a 32-bit size_t.
         * Query only: this must never attempt the multi-gigabyte allocation. */
        const int dims[] = {1 << 28, 2};
        size_t needed = 0;
        CHECK(kiss_fftndr_alloc(dims, 2, 0, NULL, &needed) == NULL);
        CHECK(needed != 0);
        if (sizeof(size_t) == 4)
            CHECK(needed == (size_t)-1);
    }
    puts("Real multidimensional allocation checks passed.");
    return 0;
}
