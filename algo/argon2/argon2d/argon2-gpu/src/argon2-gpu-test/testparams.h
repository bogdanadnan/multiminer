#ifndef TESTPARAMS_H
#define TESTPARAMS_H

#include "argon2-gpu-common/argon2params.h"

static const argon2::Argon2Params TEST_PARAMS[] = {
    { 1024, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 8,  1 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 8,  1 },
    {   31, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 8,  1 },
    {   16, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 8,  1 },
    {    9, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 8,  1 },
    {    8, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 8,  1 },

    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,   8,  1 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  16,  2 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  24,  3 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  32,  4 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 120, 15 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2, 128, 16 },

    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,   9,  1 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  10,  1 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  11,  1 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  65,  8 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  66,  8 },
    {   32, "somesalt", 8, nullptr, 0, nullptr, 0, 2,  67,  8 },
};

#endif // TESTPARAMS_H
