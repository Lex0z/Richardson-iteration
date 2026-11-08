//
// Created by Алексей on 08.11.2026.
//

#include "Richardson.h"

int check_order(const int order[3])
{
    int used[3] = {0, 0, 0};

    for (int i = 0; i < 3; i++)
    {
        if (order[i] < 0 || order[i] >= 3)
            return 0;

        if (used[order[i]])
            return 0;

        used[order[i]] = 1;
    }

    return 1;
}